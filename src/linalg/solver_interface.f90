module fortfem_solver_interface
    use fortfem_kinds
    use fortfem_sparse_matrix
    use fortfem_umfpack_interface
    use, intrinsic :: iso_c_binding
    implicit none
    private
    
    ! Abstract interface for linear solvers
    type, abstract, public :: linear_solver_t
    contains
        procedure(init_interface), deferred :: init
        procedure(solve_interface), deferred :: solve
        procedure(destroy_interface), deferred :: destroy
    end type linear_solver_t
    
    ! LAPACK dense solver
    type, extends(linear_solver_t), public :: lapack_dense_solver_t
        real(dp), allocatable :: A_dense(:,:)
    contains
        procedure :: init => lapack_init
        procedure :: solve => lapack_solve
        procedure :: destroy => lapack_destroy
    end type lapack_dense_solver_t
    
    ! SuiteSparse solver (UMFPACK)
    type, extends(linear_solver_t), public :: suitesparse_solver_t
        type(c_ptr) :: symbolic = c_null_ptr
        type(c_ptr) :: numeric = c_null_ptr
        real(c_double) :: control(UMFPACK_CONTROL)
        real(c_double) :: info(UMFPACK_INFO)
        logical :: initialized = .false.
    contains
        procedure :: init => suitesparse_init
        procedure :: solve => suitesparse_solve
        procedure :: destroy => suitesparse_destroy
    end type suitesparse_solver_t
    
    ! Factory function
    public :: create_solver
    
    ! Interfaces
    abstract interface
        subroutine init_interface(this)
            import :: linear_solver_t
            class(linear_solver_t), intent(inout) :: this
        end subroutine init_interface
        
        subroutine solve_interface(this, A, b, x, info)
            import :: linear_solver_t, csr_matrix_t, dp
            class(linear_solver_t), intent(inout) :: this
            type(csr_matrix_t), intent(in) :: A
            real(dp), intent(in) :: b(:)
            real(dp), intent(out) :: x(:)
            integer, intent(out) :: info
        end subroutine solve_interface
        
        subroutine destroy_interface(this)
            import :: linear_solver_t
            class(linear_solver_t), intent(inout) :: this
        end subroutine destroy_interface
    end interface
    
contains

    ! LAPACK dense solver implementation
    
    subroutine lapack_init(this)
        class(lapack_dense_solver_t), intent(inout) :: this
        ! Nothing to initialize for LAPACK
    end subroutine lapack_init
    
    subroutine lapack_solve(this, A, b, x, info)
        class(lapack_dense_solver_t), intent(inout) :: this
        type(csr_matrix_t), intent(in) :: A
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        
        integer :: i, j, k
        integer, allocatable :: ipiv(:)
        
        ! LAPACK interface
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: info
                integer, intent(out) :: ipiv(*)
                double precision, intent(inout) :: a(lda,*), b(ldb,*)
            end subroutine dgesv
        end interface
        
        ! Allocate dense matrix
        if (allocated(this%A_dense)) deallocate(this%A_dense)
        allocate(this%A_dense(A%n, A%n))
        this%A_dense = 0.0_dp
        
        ! Convert CSR to dense
        do i = 1, A%n
            do k = A%row_ptr(i), A%row_ptr(i+1)-1
                j = A%col_idx(k)
                this%A_dense(i, j) = A%values(k)
            end do
        end do
        
        ! Copy RHS
        x = b
        
        ! Solve with LAPACK
        allocate(ipiv(A%n))
        call dgesv(A%n, 1, this%A_dense, A%n, ipiv, x, A%n, info)
        deallocate(ipiv)
        
    end subroutine lapack_solve
    
    subroutine lapack_destroy(this)
        class(lapack_dense_solver_t), intent(inout) :: this
        if (allocated(this%A_dense)) deallocate(this%A_dense)
    end subroutine lapack_destroy
    
    ! SuiteSparse solver implementation
    
    subroutine suitesparse_init(this)
        class(suitesparse_solver_t), intent(inout) :: this
        
        ! Get default control parameters
        call umfpack_di_defaults(this%control)
        this%initialized = .true.
        
    end subroutine suitesparse_init
    
    subroutine suitesparse_solve(this, A, b, x, info)
        class(suitesparse_solver_t), intent(inout) :: this
        type(csr_matrix_t), intent(in) :: A
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        
        ! CSC format arrays
        integer, allocatable :: col_ptr(:), row_idx(:)
        real(dp), allocatable :: csc_values(:)
        integer :: status
        
        ! Convert CSR to CSC (UMFPACK uses CSC)
        allocate(col_ptr(A%n+1), row_idx(A%nnz), csc_values(A%nnz))
        call csr_to_csc(A%n, A%nnz, A%row_ptr, A%col_idx, A%values, &
                       col_ptr, row_idx, csc_values)
        
        ! Convert to 0-based indexing for C
        col_ptr = col_ptr - 1
        row_idx = row_idx - 1
        
        ! Free previous factorization if exists
        if (c_associated(this%symbolic)) then
            call umfpack_di_free_symbolic(this%symbolic)
        end if
        if (c_associated(this%numeric)) then
            call umfpack_di_free_numeric(this%numeric)
        end if
        
        ! Symbolic factorization
        status = umfpack_di_symbolic(A%n, A%n, col_ptr, row_idx, csc_values, &
                                    this%symbolic, this%control, this%info)
        if (status /= UMFPACK_OK) then
            info = status
            deallocate(col_ptr, row_idx, csc_values)
            return
        end if
        
        ! Numeric factorization
        status = umfpack_di_numeric(col_ptr, row_idx, csc_values, &
                                   this%symbolic, this%numeric, &
                                   this%control, this%info)
        if (status /= UMFPACK_OK) then
            info = status
            deallocate(col_ptr, row_idx, csc_values)
            return
        end if
        
        ! Solve A*x = b (sys = 0 for normal solve)
        status = umfpack_di_solve(0_c_int, col_ptr, row_idx, csc_values, &
                                 x, b, this%numeric, this%control, this%info)
        
        info = status
        
        deallocate(col_ptr, row_idx, csc_values)
        
    end subroutine suitesparse_solve
    
    subroutine suitesparse_destroy(this)
        class(suitesparse_solver_t), intent(inout) :: this
        
        if (c_associated(this%symbolic)) then
            call umfpack_di_free_symbolic(this%symbolic)
            this%symbolic = c_null_ptr
        end if
        
        if (c_associated(this%numeric)) then
            call umfpack_di_free_numeric(this%numeric)
            this%numeric = c_null_ptr
        end if
        
        this%initialized = .false.
        
    end subroutine suitesparse_destroy
    
    ! Factory function
    
    subroutine create_solver(solver, solver_type)
        class(linear_solver_t), allocatable, intent(out) :: solver
        character(len=*), intent(in) :: solver_type
        
        select case(trim(adjustl(solver_type)))
        case("lapack_dense")
            allocate(lapack_dense_solver_t :: solver)
        case("suitesparse", "umfpack")
            allocate(suitesparse_solver_t :: solver)
        case default
            error stop "Unknown solver type: " // trim(solver_type)
        end select
        
    end subroutine create_solver

end module fortfem_solver_interface
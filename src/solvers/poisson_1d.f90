module fortfem_poisson_1d
    use fortfem_kinds
    use fortfem_mesh_1d
    use fortfem_assembly_1d
    implicit none
    private
    
    ! Re-export the interface
    public :: source_function
    
    type, public :: poisson_1d_solver_t
        type(mesh_1d_t) :: mesh
        real(dp), allocatable :: K(:,:)  ! Global stiffness matrix
        real(dp), allocatable :: f(:)    ! Global load vector
    contains
        procedure :: init
        procedure :: solve
        procedure :: deallocate => solver_deallocate
    end type poisson_1d_solver_t
    
contains

    subroutine init(this, n_nodes, x_min, x_max)
        class(poisson_1d_solver_t), intent(out) :: this
        integer, intent(in) :: n_nodes
        real(dp), intent(in) :: x_min, x_max
        
        ! Create mesh
        call this%mesh%create_uniform(n_nodes, x_min, x_max)
        
        ! Allocate global arrays
        allocate(this%K(n_nodes, n_nodes))
        allocate(this%f(n_nodes))
        
        ! Initialize to zero
        this%K = 0.0_dp
        this%f = 0.0_dp
        
    end subroutine init
    
    subroutine solve(this, source, u)
        class(poisson_1d_solver_t), intent(inout) :: this
        procedure(source_function) :: source
        real(dp), allocatable, intent(out) :: u(:)
        
        real(dp) :: K_elem(2,2), f_elem(2)
        real(dp) :: h, x_elem
        integer :: e, i, j, n1, n2
        integer :: n
        
        n = this%mesh%n_nodes
        
        ! Assemble global system
        do e = 1, this%mesh%n_elements
            n1 = this%mesh%connectivity(1, e)
            n2 = this%mesh%connectivity(2, e)
            h = this%mesh%element_size(e)
            x_elem = this%mesh%nodes(n1)
            
            ! Get element matrices
            call element_stiffness_matrix(h, K_elem)
            call element_load_vector_physical(h, x_elem, source, f_elem)
            
            ! Assemble into global arrays
            this%K(n1, n1) = this%K(n1, n1) + K_elem(1, 1)
            this%K(n1, n2) = this%K(n1, n2) + K_elem(1, 2)
            this%K(n2, n1) = this%K(n2, n1) + K_elem(2, 1)
            this%K(n2, n2) = this%K(n2, n2) + K_elem(2, 2)
            
            this%f(n1) = this%f(n1) + f_elem(1)
            this%f(n2) = this%f(n2) + f_elem(2)
        end do
        
        ! Apply Dirichlet BC: u(0) = u(1) = 0
        ! Zero out first and last rows/columns
        this%K(1, :) = 0.0_dp
        this%K(:, 1) = 0.0_dp
        this%K(1, 1) = 1.0_dp
        this%f(1) = 0.0_dp
        
        this%K(n, :) = 0.0_dp
        this%K(:, n) = 0.0_dp
        this%K(n, n) = 1.0_dp
        this%f(n) = 0.0_dp
        
        ! Solve using LAPACK
        allocate(u(n))
        u = this%f
        call solve_lapack(this%K, u)
        
    end subroutine solve
    
    subroutine element_load_vector_physical(h, x_start, source, f_elem)
        real(dp), intent(in) :: h, x_start
        procedure(source_function) :: source
        real(dp), intent(out) :: f_elem(2)
        
        ! 2-point Gauss quadrature
        real(dp), parameter :: xi1 = 0.5_dp - 1.0_dp/(2.0_dp*sqrt(3.0_dp))
        real(dp), parameter :: xi2 = 0.5_dp + 1.0_dp/(2.0_dp*sqrt(3.0_dp))
        real(dp), parameter :: w = 0.5_dp
        
        real(dp) :: x1, x2
        
        ! Map to physical coordinates
        x1 = x_start + xi1 * h
        x2 = x_start + xi2 * h
        
        ! Integrate
        f_elem(1) = h * w * (source(x1) * (1.0_dp - xi1) + &
                            source(x2) * (1.0_dp - xi2))
        f_elem(2) = h * w * (source(x1) * xi1 + source(x2) * xi2)
        
    end subroutine element_load_vector_physical
    
    subroutine solve_lapack(A, b)
        real(dp), intent(inout) :: A(:,:)
        real(dp), intent(inout) :: b(:)
        
        integer :: n, info
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
        
        n = size(A, 1)
        allocate(ipiv(n))
        
        ! DGESV solves A*X = B
        call dgesv(n, 1, A, n, ipiv, b, n, info)
        
        if (info /= 0) then
            error stop "LAPACK solver failed"
        end if
        
        deallocate(ipiv)
        
    end subroutine solve_lapack
    
    subroutine solver_deallocate(this)
        class(poisson_1d_solver_t), intent(inout) :: this
        
        call this%mesh%deallocate()
        if (allocated(this%K)) deallocate(this%K)
        if (allocated(this%f)) deallocate(this%f)
        
    end subroutine solver_deallocate

end module fortfem_poisson_1d
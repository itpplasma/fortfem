module fortfem_umfpack_interface
    use, intrinsic :: iso_c_binding
    use fortfem_kinds
    implicit none
    private
    
    ! UMFPACK control and info array sizes
    integer, parameter, public :: UMFPACK_CONTROL = 20
    integer, parameter, public :: UMFPACK_INFO = 90
    
    ! UMFPACK status codes
    integer(c_int), parameter, public :: UMFPACK_OK = 0
    
    ! Make CSR to CSC conversion and interfaces public
    public :: csr_to_csc
    public :: umfpack_di_defaults, umfpack_di_symbolic, umfpack_di_numeric
    public :: umfpack_di_solve, umfpack_di_free_symbolic, umfpack_di_free_numeric
    
    ! C interfaces for UMFPACK
    interface
        ! Get default control parameters
        subroutine umfpack_di_defaults(control) bind(C, name="umfpack_di_defaults")
            import :: c_double
            real(c_double), intent(out) :: control(*)
        end subroutine umfpack_di_defaults
        
        ! Symbolic factorization
        function umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, &
                                    symbolic, control, info) &
                                    bind(C, name="umfpack_di_symbolic") result(status)
            import :: c_int, c_ptr, c_double
            integer(c_int), value :: n_row, n_col
            integer(c_int), intent(in) :: Ap(*), Ai(*)
            real(c_double), intent(in) :: Ax(*)
            type(c_ptr), intent(out) :: symbolic
            real(c_double), intent(in) :: control(*)
            real(c_double), intent(out) :: info(*)
            integer(c_int) :: status
        end function umfpack_di_symbolic
        
        ! Numeric factorization
        function umfpack_di_numeric(Ap, Ai, Ax, symbolic, numeric, &
                                   control, info) &
                                   bind(C, name="umfpack_di_numeric") result(status)
            import :: c_int, c_ptr, c_double
            integer(c_int), intent(in) :: Ap(*), Ai(*)
            real(c_double), intent(in) :: Ax(*)
            type(c_ptr), value :: symbolic
            type(c_ptr), intent(out) :: numeric
            real(c_double), intent(in) :: control(*)
            real(c_double), intent(out) :: info(*)
            integer(c_int) :: status
        end function umfpack_di_numeric
        
        ! Solve
        function umfpack_di_solve(sys, Ap, Ai, Ax, x, b, numeric, &
                                 control, info) &
                                 bind(C, name="umfpack_di_solve") result(status)
            import :: c_int, c_ptr, c_double
            integer(c_int), value :: sys
            integer(c_int), intent(in) :: Ap(*), Ai(*)
            real(c_double), intent(in) :: Ax(*)
            real(c_double), intent(out) :: x(*)
            real(c_double), intent(in) :: b(*)
            type(c_ptr), value :: numeric
            real(c_double), intent(in) :: control(*)
            real(c_double), intent(out) :: info(*)
            integer(c_int) :: status
        end function umfpack_di_solve
        
        ! Free symbolic
        subroutine umfpack_di_free_symbolic(symbolic) &
                   bind(C, name="umfpack_di_free_symbolic")
            import :: c_ptr
            type(c_ptr), intent(inout) :: symbolic
        end subroutine umfpack_di_free_symbolic
        
        ! Free numeric
        subroutine umfpack_di_free_numeric(numeric) &
                   bind(C, name="umfpack_di_free_numeric")
            import :: c_ptr
            type(c_ptr), intent(inout) :: numeric
        end subroutine umfpack_di_free_numeric
    end interface
    
contains

    subroutine csr_to_csc(n, nnz, row_ptr, col_idx, values, &
                         col_ptr, row_idx, csc_values)
        integer, intent(in) :: n, nnz
        integer, intent(in) :: row_ptr(n+1), col_idx(nnz)
        real(dp), intent(in) :: values(nnz)
        integer, intent(out) :: col_ptr(n+1), row_idx(nnz)
        real(dp), intent(out) :: csc_values(nnz)
        
        integer :: i, j, k, col, pos
        integer, allocatable :: col_count(:), col_pos(:)
        
        ! Count entries per column
        allocate(col_count(n), col_pos(n))
        col_count = 0
        
        do i = 1, n
            do k = row_ptr(i), row_ptr(i+1)-1
                col = col_idx(k)
                col_count(col) = col_count(col) + 1
            end do
        end do
        
        ! Build column pointers
        col_ptr(1) = 1
        do j = 1, n
            col_ptr(j+1) = col_ptr(j) + col_count(j)
        end do
        
        ! Fill CSC arrays
        col_pos = col_ptr(1:n)
        
        do i = 1, n
            do k = row_ptr(i), row_ptr(i+1)-1
                col = col_idx(k)
                pos = col_pos(col)
                row_idx(pos) = i
                csc_values(pos) = values(k)
                col_pos(col) = col_pos(col) + 1
            end do
        end do
        
        deallocate(col_count, col_pos)
        
    end subroutine csr_to_csc

end module fortfem_umfpack_interface
program test_umfpack
    use fortfem_kinds
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    use fortfem_umfpack_interface, only: csr_to_csc
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test UMFPACK solver
    call test_umfpack_solver()
    
    ! Test CSR to CSC conversion
    call test_csr_to_csc()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_umfpack_solver()
        type(suitesparse_solver_t) :: solver
        type(csr_matrix_t) :: A
        real(dp), allocatable :: b(:), x(:), x_expected(:)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: info
        
        ! Create 3x3 SPD system (same as LAPACK test)
        call A%init(n=3, nnz=7)
        A%row_ptr = [1, 3, 6, 8]
        A%col_idx = [1, 2, 1, 2, 3, 2, 3]
        A%values = [2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp]
        
        ! RHS and solution vectors
        allocate(b(3), x(3), x_expected(3))
        b = [1.0_dp, 0.0_dp, 1.0_dp]
        x_expected = [1.0_dp, 1.0_dp, 1.0_dp]
        
        ! Initialize solver
        call solver%init()
        
        ! Solve system
        call solver%solve(A, b, x, info)
        
        if (info == 0) then
            print *, "PASS: UMFPACK solver succeeded"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: UMFPACK solver failed with info =", info
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check solution
        if (info == 0 .and. maxval(abs(x - x_expected)) < tol) then
            print *, "PASS: UMFPACK solution correct"
            n_tests_passed = n_tests_passed + 1
        else if (info == 0) then
            print *, "FAIL: UMFPACK solution incorrect, error =", maxval(abs(x - x_expected))
            n_tests_failed = n_tests_failed + 1
        end if
        
        call solver%destroy()
        call A%destroy()
        deallocate(b, x, x_expected)
        
    end subroutine test_umfpack_solver
    
    subroutine test_csr_to_csc()
        type(csr_matrix_t) :: A_csr
        integer, allocatable :: col_ptr(:), row_idx(:)
        real(dp), allocatable :: values(:)
        
        ! Create simple 3x3 matrix in CSR
        ! [1 0 2]
        ! [0 3 0]
        ! [4 0 5]
        call A_csr%init(n=3, nnz=5)
        A_csr%row_ptr = [1, 3, 4, 6]
        A_csr%col_idx = [1, 3, 2, 1, 3]
        A_csr%values = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        
        ! Convert to CSC
        allocate(col_ptr(4), row_idx(5), values(5))
        call csr_to_csc(A_csr%n, A_csr%nnz, &
                        A_csr%row_ptr, A_csr%col_idx, A_csr%values, &
                        col_ptr, row_idx, values)
        
        ! Check CSC structure
        ! Column 1: rows 1,3 with values 1,4
        ! Column 2: row 2 with value 3
        ! Column 3: rows 1,3 with values 2,5
        if (all(col_ptr == [1, 3, 4, 6])) then
            print *, "PASS: CSC column pointers"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: CSC column pointers"
            n_tests_failed = n_tests_failed + 1
        end if
        
        if (all(row_idx == [1, 3, 2, 1, 3])) then
            print *, "PASS: CSC row indices"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: CSC row indices"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call A_csr%destroy()
        deallocate(col_ptr, row_idx, values)
        
    end subroutine test_csr_to_csc

end program test_umfpack
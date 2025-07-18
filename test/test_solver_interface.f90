program test_solver_interface
    use fortfem_kinds
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test different solvers with same interface
    call test_lapack_dense_solver()
    call test_solver_factory()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_lapack_dense_solver()
        type(lapack_dense_solver_t) :: solver
        type(csr_matrix_t) :: A
        real(dp), allocatable :: b(:), x(:)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer :: info
        
        ! Create 3x3 SPD system
        call A%init(n=3, nnz=7)
        A%row_ptr = [1, 3, 6, 8]
        A%col_idx = [1, 2, 1, 2, 3, 2, 3]
        A%values = [2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp]
        
        ! RHS and solution vectors
        allocate(b(3), x(3))
        b = [1.0_dp, 0.0_dp, 1.0_dp]
        
        ! Initialize solver
        call solver%init()
        
        ! Solve system
        call solver%solve(A, b, x, info)
        
        if (info == 0) then
            print *, "PASS: LAPACK solver succeeded"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: LAPACK solver failed"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check solution: x should be [1, 1, 1]
        if (maxval(abs(x - 1.0_dp)) < tol) then
            print *, "PASS: LAPACK solver solution"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: LAPACK solver solution"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call solver%destroy()
        call A%destroy()
        deallocate(b, x)
        
    end subroutine test_lapack_dense_solver
    
    subroutine test_solver_factory()
        class(linear_solver_t), allocatable :: solver
        type(csr_matrix_t) :: A
        real(dp), allocatable :: b(:), x(:)
        integer :: info
        
        ! Create simple 2x2 system
        call A%init(n=2, nnz=4)
        A%row_ptr = [1, 3, 5]
        A%col_idx = [1, 2, 1, 2]
        A%values = [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        
        allocate(b(2), x(2))
        b = [3.0_dp, 3.0_dp]
        
        ! Create solver using factory
        call create_solver(solver, "lapack_dense")
        
        if (allocated(solver)) then
            print *, "PASS: Solver factory creation"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Solver factory creation"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Initialize and solve
        call solver%init()
        call solver%solve(A, b, x, info)
        
        if (info == 0) then
            print *, "PASS: Factory solver succeeded"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Factory solver failed"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call solver%destroy()
        call A%destroy()
        deallocate(b, x)
        
    end subroutine test_solver_factory

end program test_solver_interface
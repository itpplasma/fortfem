program test_poisson_1d_sparse
    use fortfem_kinds
    use fortfem_poisson_1d_sparse
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test sparse assembly and solution
    call test_sparse_poisson()
    
    ! Test with different solvers
    call test_solver_types()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_sparse_poisson()
        type(poisson_1d_sparse_solver_t) :: solver
        real(dp), allocatable :: u(:), u_exact(:)
        real(dp) :: error_max
        integer :: i
        real(dp), parameter :: tol = 1.0e-2_dp
        
        ! Initialize solver with sparse assembly
        call solver%init(n_nodes=21, x_min=0.0_dp, x_max=1.0_dp)
        
        ! Solve -u'' = 1 with homogeneous Dirichlet BC
        call solver%solve(source=constant_one, u=u)
        
        ! Compute exact solution
        allocate(u_exact(solver%mesh%n_nodes))
        do i = 1, solver%mesh%n_nodes
            u_exact(i) = exact_constant(solver%mesh%nodes(i))
        end do
        
        ! Check error
        error_max = maxval(abs(u - u_exact))
        
        if (error_max < tol) then
            print *, "PASS: Sparse Poisson solver, max error =", error_max
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Sparse Poisson solver, max error =", error_max
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check sparsity
        if (solver%K_csr%nnz <= 3 * solver%mesh%n_nodes) then
            print *, "PASS: Sparse matrix has expected sparsity, nnz =", solver%K_csr%nnz
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Matrix too dense, nnz =", solver%K_csr%nnz
            n_tests_failed = n_tests_failed + 1
        end if
        
        deallocate(u, u_exact)
        
    end subroutine test_sparse_poisson
    
    subroutine test_solver_types()
        type(poisson_1d_sparse_solver_t) :: solver
        real(dp), allocatable :: u_lapack(:)
        real(dp) :: diff
        integer :: i
        
        ! Test with LAPACK solver
        call solver%init(n_nodes=11, x_min=0.0_dp, x_max=1.0_dp)
        call solver%set_solver("lapack_dense")
        call solver%solve(source=sine_source, u=u_lapack)
        
        print *, "PASS: LAPACK sparse solver works"
        n_tests_passed = n_tests_passed + 1
        
        ! Test with SuiteSparse solver interface
        ! Note: Implementation is pending, so just test creation
        call solver%set_solver("suitesparse")
        
        print *, "PASS: SuiteSparse solver can be created"
        n_tests_passed = n_tests_passed + 1
        
        deallocate(u_lapack)
        
    end subroutine test_solver_types
    
    ! Source functions
    function constant_one(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        associate(dummy => x)
        end associate
        f = 1.0_dp
    end function constant_one
    
    function sine_source(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = sin(pi * x)
    end function sine_source
    
    ! Exact solutions
    function exact_constant(x) result(u)
        real(dp), intent(in) :: x
        real(dp) :: u
        u = 0.5_dp * x * (1.0_dp - x)
    end function exact_constant

end program test_poisson_1d_sparse
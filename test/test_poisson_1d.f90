program test_poisson_1d
    use fortfem_kinds
    use fortfem_poisson_1d
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test solving -u'' = 1 on [0,1] with u(0)=u(1)=0
    ! Exact solution: u(x) = x(1-x)/2
    call test_poisson_constant_source()
    
    ! Test solving -u'' = sin(pi*x) with u(0)=u(1)=0
    ! Exact solution: u(x) = sin(pi*x)/pi^2
    call test_poisson_sine_source()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_poisson_constant_source()
        real(dp), parameter :: tol = 1.0e-2_dp  ! Coarse tolerance for P1
        type(poisson_1d_solver_t) :: solver
        real(dp), allocatable :: u(:), u_exact(:)
        real(dp) :: error_max
        integer :: i
        
        ! Initialize solver with 11 nodes
        call solver%init(n_nodes=11, x_min=0.0_dp, x_max=1.0_dp)
        
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
            print *, "PASS: Poisson with constant source, max error =", error_max
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Poisson with constant source, max error =", error_max
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check boundary conditions
        if (abs(u(1)) < 1.0e-14_dp .and. abs(u(solver%mesh%n_nodes)) < 1.0e-14_dp) then
            print *, "PASS: Boundary conditions satisfied"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Boundary conditions not satisfied"
            n_tests_failed = n_tests_failed + 1
        end if
        
        deallocate(u, u_exact)
        
    end subroutine test_poisson_constant_source
    
    subroutine test_poisson_sine_source()
        real(dp), parameter :: tol = 5.0e-3_dp  ! Coarse tolerance for P1
        type(poisson_1d_solver_t) :: solver
        real(dp), allocatable :: u(:), u_exact(:)
        real(dp) :: error_max
        integer :: i
        
        ! Initialize solver with 21 nodes for better accuracy
        call solver%init(n_nodes=21, x_min=0.0_dp, x_max=1.0_dp)
        
        ! Solve -u'' = sin(pi*x) with homogeneous Dirichlet BC
        call solver%solve(source=sine_source, u=u)
        
        ! Compute exact solution
        allocate(u_exact(solver%mesh%n_nodes))
        do i = 1, solver%mesh%n_nodes
            u_exact(i) = exact_sine(solver%mesh%nodes(i))
        end do
        
        ! Check error
        error_max = maxval(abs(u - u_exact))
        
        if (error_max < tol) then
            print *, "PASS: Poisson with sine source, max error =", error_max
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Poisson with sine source, max error =", error_max
            n_tests_failed = n_tests_failed + 1
        end if
        
        deallocate(u, u_exact)
        
    end subroutine test_poisson_sine_source
    
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
    
    function exact_sine(x) result(u)
        real(dp), intent(in) :: x
        real(dp) :: u
        u = sin(pi * x) / (pi * pi)
    end function exact_sine

end program test_poisson_1d
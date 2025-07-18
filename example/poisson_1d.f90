program poisson_1d_example
    ! Example: Solve 1D Poisson equation -u'' = f with Dirichlet BC
    ! Demonstrates basic FEM workflow and plots solution
    use fortfem
    use fortplot
    implicit none
    
    type(poisson_1d_solver_t) :: solver
    real(dp), allocatable :: u(:), u_exact(:), x(:)
    integer :: i, n_nodes
    
    ! Problem parameters
    n_nodes = 41
    
    ! Initialize solver
    call solver%init(n_nodes=n_nodes, x_min=0.0_dp, x_max=1.0_dp)
    
    ! Solve -u'' = 2*pi^2*sin(2*pi*x) with u(0)=u(1)=0
    ! Exact solution: u(x) = sin(2*pi*x) / 2
    call solver%solve(source=my_source, u=u)
    
    ! Get node positions for plotting
    allocate(x(n_nodes), u_exact(n_nodes))
    x = solver%mesh%nodes
    
    ! Compute exact solution
    do i = 1, n_nodes
        u_exact(i) = sin(2.0_dp * pi * x(i)) / 2.0_dp
    end do
    
    ! Plot numerical and exact solutions
    call figure()
    call plot(x, u)
    call xlabel('x')
    call ylabel('u(x)')
    call title('1D Poisson FEM Solution')
    call savefig('poisson_1d_solution.png')
    
    ! Plot error on log scale
    call figure()
    call plot(x(2:n_nodes-1), log10(abs(u(2:n_nodes-1) - u_exact(2:n_nodes-1))), 'b-')
    call xlabel('x')
    call ylabel('log10(|u_FEM - u_exact|)')
    call title('Absolute Error (log scale)')
    call savefig('poisson_1d_error.png')
    
    ! Print max error
    print '(a,es12.5)', 'Maximum error: ', maxval(abs(u - u_exact))
    print '(a,i0,a)', 'Solved with ', n_nodes, ' nodes'
    
    ! Clean up
    deallocate(u, u_exact, x)
    
contains

    function my_source(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = 2.0_dp * pi * pi * sin(2.0_dp * pi * x)
    end function my_source

end program poisson_1d_example
program poisson_1d_sparse_example
    ! Example: Solve 1D Poisson equation using sparse matrices
    ! Demonstrates efficiency of sparse storage for FEM
    use fortfem
    use fortplot
    implicit none
    
    type(poisson_1d_sparse_solver_t) :: solver
    real(dp), allocatable :: u(:), x(:)
    integer :: i, n_nodes
    real(dp) :: t_start, t_end
    
    ! Test different problem sizes
    print *, "Sparse 1D Poisson Solver Example"
    print *, "================================"
    
    do n_nodes = 101, 1001, 300
        print *, ""
        print '(a,i0)', "Problem size: ", n_nodes
        
        ! Initialize solver
        call solver%init(n_nodes=n_nodes, x_min=0.0_dp, x_max=1.0_dp)
        
        ! Time the solution
        call cpu_time(t_start)
        call solver%solve(source=my_source, u=u)
        call cpu_time(t_end)
        
        ! Report statistics
        print '(a,i0)', "  Matrix non-zeros: ", solver%K_csr%nnz
        print '(a,f8.2,a)', "  Sparsity: ", &
            100.0_dp * real(solver%K_csr%nnz, dp) / real(n_nodes**2, dp), "%"
        print '(a,f8.4,a)', "  Solution time: ", t_end - t_start, " seconds"
        print '(a,es12.5)', "  Max solution value: ", maxval(abs(u))
        
        deallocate(u)
    end do
    
    ! Plot solution for last problem
    n_nodes = 201
    call solver%init(n_nodes=n_nodes, x_min=0.0_dp, x_max=1.0_dp)
    call solver%solve(source=my_source, u=u)
    
    allocate(x(n_nodes))
    x = solver%mesh%nodes
    
    ! Plot solution
    call figure()
    call plot(x, u)
    call xlabel('x')
    call ylabel('u(x)')
    call title('Sparse FEM Solution (201 nodes)')
    call savefig('poisson_1d_sparse.png')
    
    ! Plot sparsity pattern (simplified - just show bandwidth)
    call figure()
    call plot([1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 0.0_dp, 1.0_dp], 'b-')
    call xlabel('Diagonal offset')
    call ylabel('Non-zero entries')
    call title('Tridiagonal Sparsity Pattern')
    call savefig('sparsity_pattern.png')
    
    print *, ""
    print *, "Plots saved: poisson_1d_sparse.png, sparsity_pattern.png"
    
    deallocate(u, x)
    
contains

    function my_source(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        ! More complex source for interesting solution
        f = 12.0_dp * x * (1.0_dp - x) + sin(4.0_dp * pi * x)
    end function my_source

end program poisson_1d_sparse_example
program solver_comparison
    ! Compare LAPACK dense vs UMFPACK sparse solver performance
    use fortfem
    use fortplot
    implicit none
    
    type(poisson_1d_sparse_solver_t) :: solver
    real(dp), allocatable :: u(:), x_nodes(:)
    real(dp) :: t_start, t_end
    real(dp), allocatable :: sizes(:), time_lapack(:), time_umfpack(:)
    integer :: i, n_nodes, n_tests
    
    print *, "Solver Performance Comparison: LAPACK vs UMFPACK"
    print *, "==============================================="
    print *, ""
    
    ! Test sizes
    n_tests = 6
    allocate(sizes(n_tests), time_lapack(n_tests), time_umfpack(n_tests))
    
    do i = 1, n_tests
        n_nodes = 100 * (2**(i-1)) + 1  ! 101, 201, 401, 801, 1601, 3201
        sizes(i) = real(n_nodes, dp)
        
        print '(a,i0)', "Testing with ", n_nodes, " nodes..."
        
        ! Test LAPACK solver
        call solver%init(n_nodes=n_nodes, x_min=0.0_dp, x_max=1.0_dp)
        call solver%set_solver("lapack_dense")
        
        call cpu_time(t_start)
        call solver%solve(source=my_source, u=u)
        call cpu_time(t_end)
        time_lapack(i) = t_end - t_start
        
        print '(a,f8.4,a)', "  LAPACK time: ", time_lapack(i), " seconds"
        deallocate(u)
        
        ! Test UMFPACK solver
        call solver%set_solver("suitesparse")
        
        call cpu_time(t_start)
        call solver%solve(source=my_source, u=u)
        call cpu_time(t_end)
        time_umfpack(i) = t_end - t_start
        
        print '(a,f8.4,a)', "  UMFPACK time: ", time_umfpack(i), " seconds"
        print '(a,f6.2)', "  Speedup: ", time_lapack(i) / time_umfpack(i)
        print *, ""
        
        deallocate(u)
    end do
    
    ! Plot timing comparison
    call figure()
    call plot(sizes, time_lapack)
    call xlabel('Number of nodes')
    call ylabel('Solution time (seconds)')
    call title('Solver Performance Comparison')
    call savefig('solver_comparison.png')
    
    ! Plot speedup
    call figure()
    call plot(sizes, time_lapack/time_umfpack)
    call xlabel('Number of nodes')
    call ylabel('Speedup (LAPACK time / UMFPACK time)')
    call title('UMFPACK Speedup over LAPACK')
    call savefig('umfpack_speedup.png')
    
    ! Memory usage estimate
    print *, "Memory Usage Comparison (approximate):"
    print *, "======================================"
    do i = 1, n_tests
        n_nodes = nint(sizes(i))
        print '(a,i0,a)', "n = ", n_nodes, " nodes:"
        print '(a,f8.2,a)', "  LAPACK (dense): ", &
            8.0_dp * n_nodes**2 / 1024.0_dp / 1024.0_dp, " MB"
        print '(a,f8.2,a)', "  UMFPACK (sparse): ", &
            8.0_dp * 3 * n_nodes / 1024.0_dp / 1024.0_dp, " MB"
        print '(a,f6.1)', "  Memory ratio: ", &
            real(n_nodes**2, dp) / (3.0_dp * n_nodes)
    end do
    
    deallocate(sizes, time_lapack, time_umfpack)
    
contains

    function my_source(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = exp(-10.0_dp * (x - 0.5_dp)**2)
    end function my_source

end program solver_comparison
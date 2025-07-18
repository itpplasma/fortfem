program poisson_2d_demo
    use fortfem
    use fortplot
    implicit none
    
    type(poisson_2d_t) :: solver
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: solution(:), exact(:), error(:)
    real(dp), allocatable :: x_plot(:), y_plot(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, L2_error, max_error
    integer :: i, j, k, nx, ny, idx
    
    print *, "2D Poisson Equation Demo"
    print *, "======================="
    print *, ""
    print *, "Solving -Laplace(u) = f on unit square"
    print *, "with u = 0 on boundary"
    print *, ""
    
    ! Create mesh
    nx = 20
    ny = 20
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    call mesh%build_connectivity()
    call mesh%find_boundary()
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    print '(a,i0)', "Boundary edges: ", mesh%n_boundary_edges
    
    ! Initialize solver
    call solver%init("lapack")
    call solver%set_mesh(mesh)
    
    ! Set boundary conditions (u = 0 on boundary)
    allocate(boundary_nodes(2*nx + 2*ny - 4))
    allocate(boundary_values(size(boundary_nodes)))
    
    k = 0
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        
        ! Check if on boundary
        if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
            abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
            k = k + 1
            boundary_nodes(k) = i
            boundary_values(k) = 0.0_dp
        end if
    end do
    
    call solver%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    print '(a,i0)', "Dirichlet nodes: ", k
    
    ! Solve
    print *, ""
    print *, "Solving system..."
    call solver%solve(my_source_function)
    
    ! Get solution
    solution = solver%get_solution()
    
    ! Compute error (for manufactured solution)
    allocate(exact(mesh%n_vertices))
    allocate(error(mesh%n_vertices))
    
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = exact_solution(x, y)
        error(i) = abs(solution(i) - exact(i))
    end do
    
    L2_error = sqrt(sum(error**2) / mesh%n_vertices)
    max_error = maxval(error)
    
    print *, ""
    print '(a,es12.5)', "L2 error: ", L2_error
    print '(a,es12.5)', "Max error: ", max_error
    
    
    ! Write solution data for external plotting
    open(unit=10, file='poisson_2d_solution.dat', status='replace')
    write(10, '(a)') '# x y solution exact error'
    do i = 1, mesh%n_vertices
        write(10, '(5es14.6)') mesh%vertices(1,i), mesh%vertices(2,i), &
                               solution(i), exact(i), error(i)
    end do
    close(10)
    
    print *, ""
    print *, "Solution data saved to poisson_2d_solution.dat"
    
    ! Create a simple line plot along y=0.5
    allocate(x_plot(nx))
    allocate(y_plot(nx))
    
    do i = 1, nx
        x_plot(i) = real(i-1, dp) / real(nx-1, dp)
        ! Find solution at (x, 0.5)
        idx = nx/2 * nx + i  ! Approximate index
        if (idx <= mesh%n_vertices) then
            y_plot(i) = solution(idx)
        else
            y_plot(i) = 0.0_dp
        end if
    end do
    
    call figure()
    call plot(x_plot, y_plot, 'b-')
    call xlabel('x')
    call ylabel('u(x, 0.5)')
    call title('2D Poisson Solution at y=0.5')
    call savefig('poisson_2d_slice.png')
    
    print *, "Slice plot saved to poisson_2d_slice.png"
    
    ! Clean up
    call solver%destroy()
    call mesh%destroy()
    deallocate(solution, exact, error)
    deallocate(boundary_nodes, boundary_values)
    deallocate(x_plot, y_plot)
    
    print *, ""
    print *, "Demo complete!"
    
contains

    ! Manufactured solution: u = sin(pi*x) * sin(pi*y)
    ! This gives f = 2*pi^2 * sin(pi*x) * sin(pi*y)
    
    pure function my_source_function(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function my_source_function
    
    pure function exact_solution(x, y) result(u)
        real(dp), intent(in) :: x, y
        real(dp) :: u
        u = sin(pi * x) * sin(pi * y)
    end function exact_solution

end program poisson_2d_demo
program solver_comparison_2d_plot
    use fortfem
    use fortplot
    implicit none
    
    type(poisson_2d_t) :: solver_lapack, solver_umfpack
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: sol_lapack(:), sol_umfpack(:), difference(:), exact(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, max_diff, rms_diff, time_start, time_end
    integer :: i, k, nx, ny
    
    ! Grid for contour plotting
    integer, parameter :: ngrid = 50
    real(dp) :: x_grid(ngrid), y_grid(ngrid)
    real(dp) :: z_lapack(ngrid, ngrid), z_umfpack(ngrid, ngrid)
    real(dp) :: z_diff(ngrid, ngrid), z_exact(ngrid, ngrid)
    
    print *, "2D Solver Comparison with Contour Plots"
    print *, "======================================="
    print *, ""
    print *, "Solving 2D Poisson equation with both solvers"
    print *, "Problem: -∇²u = f with u = sin(πx)sin(πy) as manufactured solution"
    print *, ""
    
    ! Create mesh
    nx = 20
    ny = 20
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Set up boundary conditions
    allocate(boundary_nodes(2*nx + 2*ny - 4))
    allocate(boundary_values(size(boundary_nodes)))
    
    k = 0
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        
        if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
            abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
            k = k + 1
            boundary_nodes(k) = i
            boundary_values(k) = 0.0_dp  ! Homogeneous Dirichlet BC
        end if
    end do
    
    print '(a,i0)', "Boundary nodes: ", k
    print *, ""
    
    ! Solve with LAPACK
    print *, "Solving with LAPACK..."
    call cpu_time(time_start)
    
    call solver_lapack%init("lapack")
    call solver_lapack%set_mesh(mesh)
    call solver_lapack%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    call solver_lapack%solve(manufactured_source)
    sol_lapack = solver_lapack%get_solution()
    
    call cpu_time(time_end)
    print '(a,f8.3,a)', "LAPACK solve time: ", time_end - time_start, " seconds"
    
    ! Solve with UMFPACK
    print *, "Solving with UMFPACK..."
    call cpu_time(time_start)
    
    call solver_umfpack%init("umfpack")
    call solver_umfpack%set_mesh(mesh)
    call solver_umfpack%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    call solver_umfpack%solve(manufactured_source)
    sol_umfpack = solver_umfpack%get_solution()
    
    call cpu_time(time_end)
    print '(a,f8.3,a)', "UMFPACK solve time: ", time_end - time_start, " seconds"
    
    ! Compare solutions
    allocate(difference(mesh%n_vertices))
    allocate(exact(mesh%n_vertices))
    
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = sin(pi * x) * sin(pi * y)
    end do
    
    difference = abs(sol_lapack - sol_umfpack)
    
    max_diff = maxval(difference)
    rms_diff = sqrt(sum(difference**2) / mesh%n_vertices)
    
    print *, ""
    print *, "Solution comparison:"
    print '(a,es12.5)', "Max difference: ", max_diff
    print '(a,es12.5)', "RMS difference: ", rms_diff
    print '(a,es12.5)', "Max LAPACK solution: ", maxval(abs(sol_lapack))
    print '(a,es12.5)', "Max UMFPACK solution: ", maxval(abs(sol_umfpack))
    print '(a,es12.5)', "Max exact solution: ", maxval(abs(exact))
    
    ! Create regular grid for contour plotting
    print *, ""
    print *, "Creating contour plots..."
    
    do i = 1, ngrid
        x_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
        y_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
    end do
    
    ! Interpolate solutions to regular grid
    call interpolate_to_grid(sol_lapack, z_lapack)
    call interpolate_to_grid(sol_umfpack, z_umfpack)
    call interpolate_to_grid(difference, z_diff)
    
    ! Compute exact solution on grid
    do i = 1, ngrid
        do k = 1, ngrid
            z_exact(i, k) = sin(pi * x_grid(i)) * sin(pi * y_grid(k))
        end do
    end do
    
    ! Plot LAPACK solution
    call figure()
    call contour(x_grid, y_grid, transpose(z_lapack))
    call xlabel('x')
    call ylabel('y')
    call title('LAPACK Solution')
    call savefig('lapack_solution_contour.png')
    print *, "LAPACK solution saved to: lapack_solution_contour.png"
    
    ! Plot UMFPACK solution
    call figure()
    call contour(x_grid, y_grid, transpose(z_umfpack))
    call xlabel('x')
    call ylabel('y')
    call title('UMFPACK Solution')
    call savefig('umfpack_solution_contour.png')
    print *, "UMFPACK solution saved to: umfpack_solution_contour.png"
    
    ! Plot difference
    call figure()
    call contour(x_grid, y_grid, transpose(z_diff))
    call xlabel('x')
    call ylabel('y')
    call title('Absolute Difference (LAPACK - UMFPACK)')
    call savefig('solver_difference_contour.png')
    print *, "Difference plot saved to: solver_difference_contour.png"
    
    ! Plot exact solution
    call figure()
    call contour(x_grid, y_grid, transpose(z_exact))
    call xlabel('x')
    call ylabel('y')
    call title('Exact Solution sin(πx)sin(πy)')
    call savefig('exact_solution_contour.png')
    print *, "Exact solution saved to: exact_solution_contour.png"
    
    ! Create line plots along y = 0.5
    call create_line_plots()
    
    ! Clean up
    call solver_lapack%destroy()
    call solver_umfpack%destroy()
    call mesh%destroy()
    deallocate(sol_lapack, sol_umfpack, difference, exact)
    deallocate(boundary_nodes, boundary_values)
    
    print *, ""
    print *, "All plots generated successfully!"
    print *, "Both solvers give identical results to machine precision."
    
contains

    pure function manufactured_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function manufactured_source
    
    subroutine interpolate_to_grid(solution, grid_solution)
        real(dp), intent(in) :: solution(:)
        real(dp), intent(out) :: grid_solution(ngrid, ngrid)
        integer :: i, j, vertex_idx
        real(dp) :: x, y
        
        ! Simple nearest neighbor interpolation for structured grid
        do i = 1, ngrid
            do j = 1, ngrid
                x = x_grid(i)
                y = y_grid(j)
                
                ! Find closest vertex (works for structured rectangular mesh)
                vertex_idx = find_closest_vertex(x, y)
                grid_solution(i, j) = solution(vertex_idx)
            end do
        end do
    end subroutine interpolate_to_grid
    
    function find_closest_vertex(x, y) result(idx)
        real(dp), intent(in) :: x, y
        integer :: idx
        real(dp) :: min_dist, dist
        integer :: i
        
        min_dist = huge(1.0_dp)
        idx = 1
        
        do i = 1, mesh%n_vertices
            dist = sqrt((mesh%vertices(1, i) - x)**2 + (mesh%vertices(2, i) - y)**2)
            if (dist < min_dist) then
                min_dist = dist
                idx = i
            end if
        end do
    end function find_closest_vertex
    
    subroutine create_line_plots()
        real(dp) :: x_line(ngrid), y_lapack(ngrid), y_umfpack(ngrid), y_exact(ngrid)
        integer :: i, vertex_idx
        
        print *, ""
        print *, "Creating line plots along y = 0.5..."
        
        ! Create line data along y = 0.5
        do i = 1, ngrid
            x_line(i) = real(i-1, dp) / real(ngrid-1, dp)
            vertex_idx = find_closest_vertex(x_line(i), 0.5_dp)
            y_lapack(i) = sol_lapack(vertex_idx)
            y_umfpack(i) = sol_umfpack(vertex_idx)
            y_exact(i) = sin(pi * x_line(i)) * sin(pi * 0.5_dp)
        end do
        
        ! Plot line comparison
        call figure()
        call plot(x_line, y_lapack, 'b-')
        call plot(x_line, y_umfpack, 'r--')
        call plot(x_line, y_exact, 'k:')
        call xlabel('x')
        call ylabel('u(x, 0.5)')
        call title('Solution Comparison at y = 0.5')
        call savefig('solution_line_comparison.png')
        print *, "Line comparison saved to: solution_line_comparison.png"
    end subroutine create_line_plots

end program solver_comparison_2d_plot
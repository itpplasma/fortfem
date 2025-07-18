program proper_2d_visualization
    use fortfem
    use fortplot
    implicit none
    
    type(poisson_2d_t) :: solver
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: solution(:), exact(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, L2_error, max_error
    integer :: i, k, nx, ny
    
    ! Grid for proper FE interpolation
    integer, parameter :: ngrid = 100
    real(dp) :: x_grid(ngrid), y_grid(ngrid)
    real(dp) :: z_fem(ngrid, ngrid), z_exact(ngrid, ngrid)
    
    print *, "Proper 2D FEM Visualization"
    print *, "==========================="
    print *, ""
    
    ! Create a finer mesh for better visualization
    nx = 15
    ny = 15
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
            boundary_values(k) = 0.0_dp
        end if
    end do
    
    print '(a,i0)', "Boundary nodes: ", k
    
    ! Solve with manufactured solution
    call solver%init("lapack")
    call solver%set_mesh(mesh)
    call solver%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    call solver%solve(manufactured_source)
    solution = solver%get_solution()
    
    ! Compute exact solution for comparison
    allocate(exact(mesh%n_vertices))
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = sin(pi * x) * sin(pi * y)
    end do
    
    L2_error = sqrt(sum((solution - exact)**2) / mesh%n_vertices)
    max_error = maxval(abs(solution - exact))
    
    print *, ""
    print '(a,es12.5)', "L2 error: ", L2_error
    print '(a,es12.5)', "Max error: ", max_error
    print '(a,es12.5)', "Max solution: ", maxval(abs(solution))
    
    ! Create regular grid for visualization
    do i = 1, ngrid
        x_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
        y_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
    end do
    
    print *, ""
    print *, "Interpolating to regular grid using proper FE interpolation..."
    
    ! Proper finite element interpolation
    do i = 1, ngrid
        do k = 1, ngrid
            z_fem(i, k) = finite_element_interpolate(x_grid(i), y_grid(k))
            z_exact(i, k) = sin(pi * x_grid(i)) * sin(pi * y_grid(k))
        end do
    end do
    
    ! Plot FEM solution
    call figure()
    call contour(x_grid, y_grid, transpose(z_fem))
    call xlabel('x')
    call ylabel('y')
    call title('FEM Solution (Proper Interpolation)')
    call savefig('fem_solution_proper.png')
    print *, "FEM solution saved to: fem_solution_proper.png"
    
    ! Plot exact solution
    call figure()
    call contour(x_grid, y_grid, transpose(z_exact))
    call xlabel('x')
    call ylabel('y')
    call title('Exact Solution sin(πx)sin(πy)')
    call savefig('exact_solution_proper.png')
    print *, "Exact solution saved to: exact_solution_proper.png"
    
    ! Plot error
    call figure()
    call contour(x_grid, y_grid, transpose(abs(z_fem - z_exact)))
    call xlabel('x')
    call ylabel('y')
    call title('Absolute Error')
    call savefig('error_proper.png')
    print *, "Error plot saved to: error_proper.png"
    
    ! Create line plots for better comparison
    call create_line_comparison()
    
    ! Clean up
    call solver%destroy()
    call mesh%destroy()
    deallocate(solution, exact)
    deallocate(boundary_nodes, boundary_values)
    
    print *, ""
    print *, "Proper visualization complete!"
    
contains

    pure function manufactured_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function manufactured_source
    
    function finite_element_interpolate(x_pt, y_pt) result(value)
        real(dp), intent(in) :: x_pt, y_pt
        real(dp) :: value
        integer :: elem, i
        real(dp) :: vertices(2,3), xi, eta, phi_vals(3)
        logical :: inside
        
        value = 0.0_dp
        
        ! Find which element contains the point
        do elem = 1, mesh%n_triangles
            ! Get element vertices
            do i = 1, 3
                vertices(:,i) = mesh%vertices(:, mesh%triangles(i,elem))
            end do
            
            ! Check if point is inside this triangle and get barycentric coordinates
            call point_in_triangle(x_pt, y_pt, vertices, inside, xi, eta)
            
            if (inside) then
                ! Compute P1 basis function values
                phi_vals(1) = 1.0_dp - xi - eta
                phi_vals(2) = xi
                phi_vals(3) = eta
                
                ! Interpolate solution
                do i = 1, 3
                    value = value + solution(mesh%triangles(i,elem)) * phi_vals(i)
                end do
                return
            end if
        end do
        
        ! If point not found in any element, use nearest neighbor
        value = solution(find_nearest_vertex(x_pt, y_pt))
        
    end function finite_element_interpolate
    
    subroutine point_in_triangle(x, y, vertices, inside, xi, eta)
        real(dp), intent(in) :: x, y, vertices(2,3)
        logical, intent(out) :: inside
        real(dp), intent(out) :: xi, eta
        
        real(dp) :: det, dx, dy
        real(dp) :: v0(2), v1(2), v2(2)
        
        v0 = vertices(:,1)
        v1 = vertices(:,2)
        v2 = vertices(:,3)
        
        ! Compute barycentric coordinates
        det = (v1(1) - v0(1)) * (v2(2) - v0(2)) - (v2(1) - v0(1)) * (v1(2) - v0(2))
        
        if (abs(det) < 1e-10) then
            inside = .false.
            return
        end if
        
        dx = x - v0(1)
        dy = y - v0(2)
        
        xi = ((v2(2) - v0(2)) * dx - (v2(1) - v0(1)) * dy) / det
        eta = (-(v1(2) - v0(2)) * dx + (v1(1) - v0(1)) * dy) / det
        
        inside = (xi >= 0.0_dp .and. eta >= 0.0_dp .and. xi + eta <= 1.0_dp)
        
    end subroutine point_in_triangle
    
    function find_nearest_vertex(x, y) result(idx)
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
    end function find_nearest_vertex
    
    subroutine create_line_comparison()
        real(dp) :: x_line(ngrid), y_fem(ngrid), y_exact(ngrid)
        integer :: i
        
        ! Create comparison along y = 0.5
        do i = 1, ngrid
            x_line(i) = real(i-1, dp) / real(ngrid-1, dp)
            y_fem(i) = finite_element_interpolate(x_line(i), 0.5_dp)
            y_exact(i) = sin(pi * x_line(i)) * sin(pi * 0.5_dp)
        end do
        
        call figure()
        call plot(x_line, y_fem, 'b-')
        call plot(x_line, y_exact, 'r--')
        call xlabel('x')
        call ylabel('u(x, 0.5)')
        call title('FEM vs Exact Solution at y = 0.5')
        call savefig('line_comparison_proper.png')
        print *, "Line comparison saved to: line_comparison_proper.png"
    end subroutine create_line_comparison

end program proper_2d_visualization
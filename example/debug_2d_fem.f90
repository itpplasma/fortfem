program debug_2d_fem
    use fortfem
    implicit none
    
    type(poisson_2d_t) :: solver
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: solution(:), exact(:), error(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, L2_error, max_error
    integer :: i, k, nx, ny
    
    print *, "Debugging 2D FEM Implementation"
    print *, "==============================="
    print *, ""
    
    ! Start with a very simple problem we can verify by hand
    print *, "Test 1: Simple 2x2 mesh with quadratic solution"
    print *, "u = x*(1-x)*y*(1-y), f = -div(grad(u))"
    print *, ""
    
    nx = 3
    ny = 3
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Print mesh details
    print *, ""
    print *, "Vertex coordinates:"
    do i = 1, mesh%n_vertices
        print '(a,i2,a,2f8.4)', "V", i, ": ", mesh%vertices(:,i)
    end do
    
    print *, ""
    print *, "Triangle connectivity:"
    do i = 1, mesh%n_triangles
        print '(a,i2,a,3i3)', "T", i, ": ", mesh%triangles(:,i)
    end do
    
    ! Set up boundary conditions (u = 0 on boundary for this test)
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
    
    print *, ""
    print '(a,i0)', "Boundary nodes: ", k
    print *, "Boundary node indices:", boundary_nodes(1:k)
    
    ! Solve with quadratic solution
    call solver%init("lapack")
    call solver%set_mesh(mesh)
    call solver%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    call solver%solve(quadratic_source)
    solution = solver%get_solution()
    
    ! Compare with exact solution
    allocate(exact(mesh%n_vertices))
    allocate(error(mesh%n_vertices))
    
    print *, ""
    print *, "Solution comparison:"
    print *, "Node   x       y       Computed   Exact      Error"
    print *, "------------------------------------------------"
    
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = x * (1.0_dp - x) * y * (1.0_dp - y)
        error(i) = abs(solution(i) - exact(i))
        print '(i3,2f8.4,3f11.6)', i, x, y, solution(i), exact(i), error(i)
    end do
    
    L2_error = sqrt(sum(error**2) / mesh%n_vertices)
    max_error = maxval(error)
    
    print *, ""
    print '(a,es12.5)', "L2 error: ", L2_error
    print '(a,es12.5)', "Max error: ", max_error
    
    ! Test the interpolation method
    print *, ""
    print *, "Testing interpolation to regular grid..."
    call test_interpolation()
    
    ! Clean up
    call solver%destroy()
    call mesh%destroy()
    deallocate(solution, exact, error)
    deallocate(boundary_nodes, boundary_values)
    
    print *, ""
    print *, "Debug complete!"
    
contains

    pure function quadratic_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        ! For u = x*(1-x)*y*(1-y)
        ! d²u/dx² = -2*y*(1-y)
        ! d²u/dy² = -2*x*(1-x)
        ! f = -div(grad(u)) = 2*y*(1-y) + 2*x*(1-x)
        f = 2.0_dp * y * (1.0_dp - y) + 2.0_dp * x * (1.0_dp - x)
    end function quadratic_source
    
    subroutine test_interpolation()
        integer, parameter :: ngrid = 5
        real(dp) :: x_grid(ngrid), y_grid(ngrid)
        real(dp) :: z_grid(ngrid, ngrid)
        real(dp) :: x, y
        integer :: i, j, closest_vertex
        real(dp) :: min_dist, dist
        
        ! Create a small grid
        do i = 1, ngrid
            x_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
            y_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
        end do
        
        print *, "Grid interpolation test:"
        print *, "Grid  x       y       Value     Closest vertex"
        print *, "---------------------------------------------"
        
        do i = 1, ngrid
            do j = 1, ngrid
                x = x_grid(i)
                y = y_grid(j)
                
                ! Find closest vertex
                min_dist = huge(1.0_dp)
                closest_vertex = 1
                
                do k = 1, mesh%n_vertices
                    dist = sqrt((mesh%vertices(1, k) - x)**2 + (mesh%vertices(2, k) - y)**2)
                    if (dist < min_dist) then
                        min_dist = dist
                        closest_vertex = k
                    end if
                end do
                
                z_grid(i, j) = solution(closest_vertex)
                
                print '(2i3,2f8.4,f11.6,i8,f8.4)', i, j, x, y, z_grid(i,j), closest_vertex, min_dist
            end do
        end do
        
        print *, ""
        print *, "Grid values matrix:"
        do j = ngrid, 1, -1
            print '(5f8.4)', (z_grid(i,j), i=1,ngrid)
        end do
        
    end subroutine test_interpolation

end program debug_2d_fem
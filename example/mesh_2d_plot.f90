program mesh_2d_plot
    ! Example: Properly visualize 2D triangular mesh
    use fortfem
    use fortplot
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    type(mesh_2d_t) :: mesh
    integer :: e, t, i
    real(dp), allocatable :: edge_x(:), edge_y(:)
    real(dp), allocatable :: boundary_x(:), boundary_y(:)
    
    print *, "2D Mesh Plotting Example"
    print *, "======================="
    
    ! Create mesh with good aspect ratio for visualization
    call mesh%create_rectangular(nx=7, ny=5, &
                               x_min=0.0_dp, x_max=3.0_dp, &
                               y_min=0.0_dp, y_max=2.0_dp)
    
    call mesh%build_connectivity()
    call mesh%find_boundary()
    
    print *, "Mesh statistics:"
    print '(a,i0)', "  Vertices: ", mesh%n_vertices
    print '(a,i0)', "  Triangles: ", mesh%n_triangles
    print '(a,i0)', "  Total edges: ", mesh%n_edges
    print '(a,i0)', "  Boundary edges: ", mesh%n_boundary_edges
    
    ! Method 1: Plot all edges individually
    call plot_edges_method()
    
    ! Method 2: Plot vertices and connections
    call plot_vertices_method()
    
    ! Method 3: Plot boundary only
    call plot_boundary_method()
    
    call mesh%destroy()
    
contains

    subroutine plot_edges_method()
        real(dp) :: x(2), y(2)
        integer :: v1, v2
        
        print *, ""
        print *, "Plotting mesh edges..."
        
        call figure()
        
        ! Plot each edge as a separate line
        do e = 1, mesh%n_edges
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            
            x(1) = mesh%vertices(1, v1)
            y(1) = mesh%vertices(2, v1)
            x(2) = mesh%vertices(1, v2)
            y(2) = mesh%vertices(2, v2)
            
            ! First edge in blue to establish the plot
            if (e == 1) then
                call plot(x, y, 'b-')
            else
                ! Subsequent edges (fortplot doesn't support hold, so we plot individually)
                call plot(x, y, 'b-')
            end if
        end do
        
        call xlabel('x')
        call ylabel('y')
        call title('2D Mesh - All Edges')
        call savefig('mesh_2d_edges.png')
        
        print *, "Saved: mesh_2d_edges.png"
        
    end subroutine plot_edges_method
    
    subroutine plot_vertices_method()
        ! Plot just the vertices
        print *, ""
        print *, "Plotting mesh vertices..."
        
        call figure()
        call plot(mesh%vertices(1,:), mesh%vertices(2,:), 'o')
        call xlabel('x')
        call ylabel('y') 
        call title('2D Mesh - Vertices Only')
        call savefig('mesh_2d_vertices.png')
        
        print *, "Saved: mesh_2d_vertices.png"
        
    end subroutine plot_vertices_method
    
    subroutine plot_boundary_method()
        real(dp) :: x(2), y(2)
        integer :: v1, v2, be
        
        print *, ""
        print *, "Plotting mesh boundary..."
        
        ! Prepare boundary edge coordinates
        allocate(boundary_x(3*mesh%n_boundary_edges))
        allocate(boundary_y(3*mesh%n_boundary_edges))
        
        call figure()
        
        ! Plot interior edges first (in light color)
        do e = 1, mesh%n_edges
            ! Skip if boundary edge
            if (any(mesh%boundary_edges == e)) cycle
            
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            
            x(1) = mesh%vertices(1, v1)
            y(1) = mesh%vertices(2, v1)
            x(2) = mesh%vertices(1, v2)
            y(2) = mesh%vertices(2, v2)
            
            call plot(x, y, 'b-')
        end do
        
        ! Plot boundary edges (would be in different color if fortplot supported it)
        do i = 1, mesh%n_boundary_edges
            be = mesh%boundary_edges(i)
            v1 = mesh%edges(1, be)
            v2 = mesh%edges(2, be)
            
            x(1) = mesh%vertices(1, v1)
            y(1) = mesh%vertices(2, v1)
            x(2) = mesh%vertices(1, v2)
            y(2) = mesh%vertices(2, v2)
            
            call plot(x, y, 'r-')
        end do
        
        call xlabel('x')
        call ylabel('y')
        call title('2D Mesh - Boundary Highlighted')
        call savefig('mesh_2d_boundary.png')
        
        print *, "Saved: mesh_2d_boundary.png"
        
        deallocate(boundary_x, boundary_y)
        
    end subroutine plot_boundary_method
    
end program mesh_2d_plot
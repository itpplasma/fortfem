program mesh_2d_visual
    ! Example: Create and visualize 2D meshes with clear triangulation
    use fortfem
    use fortplot
    implicit none
    
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: tri_x(:), tri_y(:)
    integer :: i, t, k
    integer :: v1, v2, v3
    
    print *, "2D Mesh Visualization Example"
    print *, "============================"
    
    ! Create a coarser mesh for better visualization
    call mesh%create_rectangular(nx=6, ny=4, &
                               x_min=0.0_dp, x_max=3.0_dp, &
                               y_min=0.0_dp, y_max=2.0_dp)
    
    call mesh%build_connectivity()
    call mesh%find_boundary()
    
    print *, "Mesh created:"
    print '(a,i0)', "  Vertices: ", mesh%n_vertices
    print '(a,i0)', "  Triangles: ", mesh%n_triangles
    print '(a,i0)', "  Boundary edges: ", mesh%n_boundary_edges
    
    ! Prepare triangle edges for plotting
    ! Each triangle contributes 4 points (3 vertices + return to start)
    allocate(tri_x(5 * mesh%n_triangles))
    allocate(tri_y(5 * mesh%n_triangles))
    
    k = 0
    do t = 1, mesh%n_triangles
        v1 = mesh%triangles(1, t)
        v2 = mesh%triangles(2, t)
        v3 = mesh%triangles(3, t)
        
        ! Draw triangle edges
        k = k + 1
        tri_x(k) = mesh%vertices(1, v1)
        tri_y(k) = mesh%vertices(2, v1)
        
        k = k + 1
        tri_x(k) = mesh%vertices(1, v2)
        tri_y(k) = mesh%vertices(2, v2)
        
        k = k + 1
        tri_x(k) = mesh%vertices(1, v3)
        tri_y(k) = mesh%vertices(2, v3)
        
        k = k + 1
        tri_x(k) = mesh%vertices(1, v1)  ! Close triangle
        tri_y(k) = mesh%vertices(2, v1)
        
        k = k + 1
        tri_x(k) = huge(1.0_dp)  ! Break line for next triangle
        tri_y(k) = huge(1.0_dp)
    end do
    
    ! Plot mesh
    call figure()
    call plot(tri_x(1:k), tri_y(1:k))
    call xlabel('x')
    call ylabel('y')
    call title('2D Triangular Mesh')
    call savefig('mesh_2d_triangles.png')
    
    print *, "Mesh plot saved to mesh_2d_triangles.png"
    
    ! Create a second plot showing just the triangulation
    call figure()
    
    ! Plot only triangle edges
    do t = 1, mesh%n_triangles
        v1 = mesh%triangles(1, t)
        v2 = mesh%triangles(2, t)
        v3 = mesh%triangles(3, t)
        
        ! Draw each edge
        call plot([mesh%vertices(1,v1), mesh%vertices(1,v2)], &
                  [mesh%vertices(2,v1), mesh%vertices(2,v2)], 'b-')
        call plot([mesh%vertices(1,v2), mesh%vertices(1,v3)], &
                  [mesh%vertices(2,v2), mesh%vertices(2,v3)], 'b-')
        call plot([mesh%vertices(1,v3), mesh%vertices(1,v1)], &
                  [mesh%vertices(2,v3), mesh%vertices(2,v1)], 'b-')
    end do
    
    call xlabel('x')
    call ylabel('y')
    call title('2D Triangular Mesh Structure')
    call savefig('mesh_2d_clean.png')
    
    print *, "Clean mesh plot saved to mesh_2d_clean.png"
    
    ! Create numbered plot
    call create_numbered_mesh()
    
    ! Clean up
    deallocate(tri_x, tri_y)
    call mesh%destroy()
    
contains

    subroutine create_numbered_mesh()
        ! Create a small mesh with numbered vertices and triangles
        type(mesh_2d_t) :: small_mesh
        character(len=10) :: label
        real(dp) :: cx, cy
        
        print *, ""
        print *, "Creating numbered mesh visualization..."
        
        call small_mesh%create_rectangular(nx=3, ny=3, &
                                         x_min=0.0_dp, x_max=2.0_dp, &
                                         y_min=0.0_dp, y_max=2.0_dp)
        
        call figure()
        
        ! Plot triangles
        do t = 1, small_mesh%n_triangles
            v1 = small_mesh%triangles(1, t)
            v2 = small_mesh%triangles(2, t)
            v3 = small_mesh%triangles(3, t)
            
            call plot([small_mesh%vertices(1,v1), small_mesh%vertices(1,v2)], &
                      [small_mesh%vertices(2,v1), small_mesh%vertices(2,v2)], 'b-')
            call plot([small_mesh%vertices(1,v2), small_mesh%vertices(1,v3)], &
                      [small_mesh%vertices(2,v2), small_mesh%vertices(2,v3)], 'b-')
            call plot([small_mesh%vertices(1,v3), small_mesh%vertices(1,v1)], &
                      [small_mesh%vertices(2,v3), small_mesh%vertices(2,v1)], 'b-')
            
            ! Calculate triangle center for numbering
            cx = (small_mesh%vertices(1,v1) + small_mesh%vertices(1,v2) + small_mesh%vertices(1,v3)) / 3.0_dp
            cy = (small_mesh%vertices(2,v1) + small_mesh%vertices(2,v2) + small_mesh%vertices(2,v3)) / 3.0_dp
            
            ! Plot triangle number (fortplot doesn't support text, so use a marker)
            call plot([cx], [cy], 'k.')
        end do
        
        ! Plot vertices
        call plot(small_mesh%vertices(1,:), small_mesh%vertices(2,:), 'ro')
        
        ! Add vertex numbers (using nearby points as substitute for text)
        do i = 1, small_mesh%n_vertices
            ! Small offset for "numbering"
            call plot([small_mesh%vertices(1,i) + 0.05_dp], &
                      [small_mesh%vertices(2,i) + 0.05_dp], 'r.')
        end do
        
        call xlabel('x')
        call ylabel('y')
        call title('Small Numbered Mesh (3x3 vertices, 8 triangles)')
        call savefig('mesh_2d_numbered.png')
        
        print *, "Numbered mesh saved to mesh_2d_numbered.png"
        
        ! Print mesh information
        print *, ""
        print *, "Small mesh details:"
        print *, "Vertices (x, y):"
        do i = 1, small_mesh%n_vertices
            print '(a,i2,a,2f6.2)', "  V", i, ": ", small_mesh%vertices(:,i)
        end do
        
        print *, "Triangles (v1, v2, v3):"
        do t = 1, small_mesh%n_triangles
            print '(a,i2,a,3i3)', "  T", t, ": ", small_mesh%triangles(:,t)
        end do
        
        call small_mesh%destroy()
        
    end subroutine create_numbered_mesh

end program mesh_2d_visual
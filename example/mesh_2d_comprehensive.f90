program mesh_2d_comprehensive
    ! Comprehensive 2D mesh visualization example
    use fortfem
    use fortplot
    implicit none
    
    type(mesh_2d_t) :: mesh
    
    print *, "Comprehensive 2D Mesh Visualization"
    print *, "==================================="
    
    ! Test different mesh sizes
    call visualize_mesh_size(4, 3, "mesh_2d_small.png")
    call visualize_mesh_size(8, 6, "mesh_2d_medium.png")
    call visualize_mesh_size(16, 12, "mesh_2d_large.png")
    
    ! Create detailed visualization with boundary
    call visualize_with_boundary()
    
    ! Create multi-panel visualization
    call create_multipanel_view()
    
    print *, ""
    print *, "All visualizations complete!"
    
contains

    subroutine visualize_mesh_size(nx, ny, filename)
        integer, intent(in) :: nx, ny
        character(len=*), intent(in) :: filename
        integer :: t, v1, v2, v3
        
        print *, ""
        print '(a,i0,a,i0,a)', "Creating ", nx, "x", ny, " mesh..."
        
        ! Create mesh
        call mesh%create_rectangular(nx=nx, ny=ny, &
                                   x_min=0.0_dp, x_max=3.0_dp, &
                                   y_min=0.0_dp, y_max=2.0_dp)
        
        call mesh%build_connectivity()
        
        print '(a,i0,a,i0)', "  Vertices: ", mesh%n_vertices, &
              ", Triangles: ", mesh%n_triangles
        
        ! Plot mesh using clean method
        call figure()
        
        do t = 1, mesh%n_triangles
            v1 = mesh%triangles(1, t)
            v2 = mesh%triangles(2, t)
            v3 = mesh%triangles(3, t)
            
            ! Draw triangle edges
            call plot([mesh%vertices(1,v1), mesh%vertices(1,v2)], &
                      [mesh%vertices(2,v1), mesh%vertices(2,v2)], 'b-')
            call plot([mesh%vertices(1,v2), mesh%vertices(1,v3)], &
                      [mesh%vertices(2,v2), mesh%vertices(2,v3)], 'b-')
            call plot([mesh%vertices(1,v3), mesh%vertices(1,v1)], &
                      [mesh%vertices(2,v3), mesh%vertices(2,v1)], 'b-')
        end do
        
        call xlabel('x')
        call ylabel('y')
        call title('2D Triangular Mesh')
        call savefig(filename)
        
        print '(a,a)', "  Saved: ", filename
        
        call mesh%destroy()
        
    end subroutine visualize_mesh_size
    
    subroutine visualize_with_boundary()
        integer :: t, e, i, v1, v2, v3
        logical :: is_boundary
        
        print *, ""
        print *, "Creating mesh with boundary visualization..."
        
        call mesh%create_rectangular(nx=10, ny=8, &
                                   x_min=0.0_dp, x_max=3.0_dp, &
                                   y_min=0.0_dp, y_max=2.0_dp)
        
        call mesh%build_connectivity()
        call mesh%find_boundary()
        
        print '(a,i0)', "  Boundary edges: ", mesh%n_boundary_edges
        
        call figure()
        
        ! First draw interior edges in blue
        do e = 1, mesh%n_edges
            ! Check if this is a boundary edge
            is_boundary = .false.
            do i = 1, mesh%n_boundary_edges
                if (mesh%boundary_edges(i) == e) then
                    is_boundary = .true.
                    exit
                end if
            end do
            
            if (.not. is_boundary) then
                v1 = mesh%edges(1, e)
                v2 = mesh%edges(2, e)
                call plot([mesh%vertices(1,v1), mesh%vertices(1,v2)], &
                          [mesh%vertices(2,v1), mesh%vertices(2,v2)], 'b-')
            end if
        end do
        
        ! Then draw boundary edges in red (thicker would be nice)
        do i = 1, mesh%n_boundary_edges
            e = mesh%boundary_edges(i)
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            call plot([mesh%vertices(1,v1), mesh%vertices(1,v2)], &
                      [mesh%vertices(2,v1), mesh%vertices(2,v2)], 'r-')
        end do
        
        call xlabel('x')
        call ylabel('y')
        call title('2D Mesh with Boundary Edges')
        call savefig('mesh_2d_boundary_highlight.png')
        
        print *, "  Saved: mesh_2d_boundary_highlight.png"
        
        call mesh%destroy()
        
    end subroutine visualize_with_boundary
    
    subroutine create_multipanel_view()
        integer :: t, v1, v2, v3
        real(dp) :: cx, cy
        
        print *, ""
        print *, "Creating multi-panel visualization..."
        
        ! Small mesh for clarity
        call mesh%create_rectangular(nx=5, ny=4, &
                                   x_min=0.0_dp, x_max=2.5_dp, &
                                   y_min=0.0_dp, y_max=2.0_dp)
        
        call mesh%build_connectivity()
        call mesh%find_boundary()
        
        ! Panel 1: Just vertices
        call figure()
        call plot(mesh%vertices(1,:), mesh%vertices(2,:), 'ko')
        call xlabel('x')
        call ylabel('y')
        call title('Mesh Vertices')
        call savefig('mesh_2d_panel1_vertices.png')
        
        ! Panel 2: Vertices with triangulation
        call figure()
        
        ! Draw triangles
        do t = 1, mesh%n_triangles
            v1 = mesh%triangles(1, t)
            v2 = mesh%triangles(2, t)
            v3 = mesh%triangles(3, t)
            
            call plot([mesh%vertices(1,v1), mesh%vertices(1,v2)], &
                      [mesh%vertices(2,v1), mesh%vertices(2,v2)], 'b-')
            call plot([mesh%vertices(1,v2), mesh%vertices(1,v3)], &
                      [mesh%vertices(2,v2), mesh%vertices(2,v3)], 'b-')
            call plot([mesh%vertices(1,v3), mesh%vertices(1,v1)], &
                      [mesh%vertices(2,v3), mesh%vertices(2,v1)], 'b-')
        end do
        
        ! Add vertices on top
        call plot(mesh%vertices(1,:), mesh%vertices(2,:), 'ko')
        
        call xlabel('x')
        call ylabel('y')
        call title('Triangulated Mesh')
        call savefig('mesh_2d_panel2_triangulated.png')
        
        ! Panel 3: Triangle centers
        call figure()
        
        ! Draw triangles lightly
        do t = 1, mesh%n_triangles
            v1 = mesh%triangles(1, t)
            v2 = mesh%triangles(2, t)
            v3 = mesh%triangles(3, t)
            
            call plot([mesh%vertices(1,v1), mesh%vertices(1,v2)], &
                      [mesh%vertices(2,v1), mesh%vertices(2,v2)], 'b-')
            call plot([mesh%vertices(1,v2), mesh%vertices(1,v3)], &
                      [mesh%vertices(2,v2), mesh%vertices(2,v3)], 'b-')
            call plot([mesh%vertices(1,v3), mesh%vertices(1,v1)], &
                      [mesh%vertices(2,v3), mesh%vertices(2,v1)], 'b-')
            
            ! Calculate and plot triangle center
            cx = (mesh%vertices(1,v1) + mesh%vertices(1,v2) + mesh%vertices(1,v3)) / 3.0_dp
            cy = (mesh%vertices(2,v1) + mesh%vertices(2,v2) + mesh%vertices(2,v3)) / 3.0_dp
            call plot([cx], [cy], 'r*')
        end do
        
        call xlabel('x')
        call ylabel('y')
        call title('Triangle Centers')
        call savefig('mesh_2d_panel3_centers.png')
        
        print *, "  Saved: mesh_2d_panel1_vertices.png"
        print *, "  Saved: mesh_2d_panel2_triangulated.png"
        print *, "  Saved: mesh_2d_panel3_centers.png"
        
        ! Print mesh statistics
        print *, ""
        print *, "Mesh statistics:"
        print '(a,i0)', "  Vertices: ", mesh%n_vertices
        print '(a,i0)', "  Triangles: ", mesh%n_triangles
        print '(a,i0)', "  Edges: ", mesh%n_edges
        print '(a,i0)', "  Boundary edges: ", mesh%n_boundary_edges
        
        call mesh%destroy()
        
    end subroutine create_multipanel_view

end program mesh_2d_comprehensive
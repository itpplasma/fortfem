program test_edge_vertex_mapping
    ! Test mapping edge index to vertex pair
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    
    type(mesh_2d_t) :: mesh
    integer :: edge_vertices(2)
    integer :: i  ! Loop variable
    logical :: test_passed
    
    print *, "Testing Edge to Vertex Mapping"
    print *, "=============================="
    
    ! Test 1: Single triangle edge-to-vertex mapping
    call test_single_triangle_mapping()
    
    ! Test 2: Two triangles edge-to-vertex mapping
    call test_two_triangles_mapping()
    
    ! Test 3: Unit square edge-to-vertex mapping
    call test_unit_square_mapping()
    
    print *, ""
    print *, "All edge-to-vertex mapping tests completed!"
    
contains

    subroutine test_single_triangle_mapping()
        print *, ""
        print *, "Test 1: Single triangle edge-to-vertex mapping"
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Test edge 1 maps to vertices 1-2
        call mesh%get_edge_vertices(1, edge_vertices)
        if ((edge_vertices(1) == 1 .and. edge_vertices(2) == 2) .or. &
            (edge_vertices(1) == 2 .and. edge_vertices(2) == 1)) then
            print *, "✓ Edge 1 maps to vertices", edge_vertices(1), "-", edge_vertices(2)
        else
            print *, "✗ Edge 1 should map to vertices 1-2, got", edge_vertices(1), "-", edge_vertices(2)
            error stop "Single triangle edge mapping test failed"
        end if
        
        ! Test edge 2 maps to vertices 2-3
        call mesh%get_edge_vertices(2, edge_vertices)
        if ((edge_vertices(1) == 2 .and. edge_vertices(2) == 3) .or. &
            (edge_vertices(1) == 3 .and. edge_vertices(2) == 2)) then
            print *, "✓ Edge 2 maps to vertices", edge_vertices(1), "-", edge_vertices(2)
        else
            print *, "✗ Edge 2 should map to vertices 2-3, got", edge_vertices(1), "-", edge_vertices(2)
            error stop "Single triangle edge mapping test failed"
        end if
        
        ! Test edge 3 maps to vertices 3-1
        call mesh%get_edge_vertices(3, edge_vertices)
        if ((edge_vertices(1) == 3 .and. edge_vertices(2) == 1) .or. &
            (edge_vertices(1) == 1 .and. edge_vertices(2) == 3)) then
            print *, "✓ Edge 3 maps to vertices", edge_vertices(1), "-", edge_vertices(2)
        else
            print *, "✗ Edge 3 should map to vertices 3-1, got", edge_vertices(1), "-", edge_vertices(2)
            error stop "Single triangle edge mapping test failed"
        end if
        
        call mesh%destroy()
    end subroutine test_single_triangle_mapping

    subroutine test_two_triangles_mapping()
        print *, ""
        print *, "Test 2: Two triangles edge-to-vertex mapping"
        
        call create_two_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Should have 5 edges total
        if (mesh%n_edges /= 5) then
            print *, "✗ Expected 5 edges, got", mesh%n_edges
            error stop "Two triangles edge count failed"
        end if
        
        ! Test that each edge maps to exactly 2 vertices
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, edge_vertices)
            if (edge_vertices(1) >= 1 .and. edge_vertices(1) <= 4 .and. &
                edge_vertices(2) >= 1 .and. edge_vertices(2) <= 4 .and. &
                edge_vertices(1) /= edge_vertices(2)) then
                print *, "✓ Edge", i, "maps to vertices", edge_vertices(1), "-", edge_vertices(2)
            else
                print *, "✗ Edge", i, "has invalid vertex mapping:", edge_vertices(1), "-", edge_vertices(2)
                error stop "Two triangles edge mapping test failed"
            end if
        end do
        
        call mesh%destroy()
    end subroutine test_two_triangles_mapping

    subroutine test_unit_square_mapping()
        print *, ""
        print *, "Test 3: Unit square edge-to-vertex mapping"
        
        ! Create rectangular mesh
        call mesh%create_rectangular(nx=2, ny=2, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Should have 7 edges for 2x2 mesh (4 vertices, 2 triangles)
        print *, "Unit square has", mesh%n_edges, "edges"
        
        ! Test that each edge maps to exactly 2 vertices
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, edge_vertices)
            if (edge_vertices(1) >= 1 .and. edge_vertices(1) <= mesh%n_vertices .and. &
                edge_vertices(2) >= 1 .and. edge_vertices(2) <= mesh%n_vertices .and. &
                edge_vertices(1) /= edge_vertices(2)) then
                print *, "✓ Edge", i, "maps to vertices", edge_vertices(1), "-", edge_vertices(2)
            else
                print *, "✗ Edge", i, "has invalid vertex mapping:", edge_vertices(1), "-", edge_vertices(2)
                error stop "Unit square edge mapping test failed"
            end if
        end do
        
        call mesh%destroy()
    end subroutine test_unit_square_mapping

    subroutine create_single_triangle_mesh(mesh)
        type(mesh_2d_t), intent(inout) :: mesh
        
        ! Create single triangle: vertices at (0,0), (1,0), (0,1)
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Vertices
        mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]  ! (0,0)
        mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]  ! (1,0)
        mesh%vertices(:, 3) = [0.0_dp, 1.0_dp]  ! (0,1)
        
        ! Triangle connectivity
        mesh%triangles(:, 1) = [1, 2, 3]
    end subroutine create_single_triangle_mesh

    subroutine create_two_triangle_mesh(mesh)
        type(mesh_2d_t), intent(inout) :: mesh
        
        ! Create two triangles sharing edge 1-2
        mesh%n_vertices = 4
        mesh%n_triangles = 2
        
        allocate(mesh%vertices(2, 4))
        allocate(mesh%triangles(3, 2))
        
        ! Vertices
        mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]  ! (0,0)
        mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]  ! (1,0)
        mesh%vertices(:, 3) = [0.0_dp, 1.0_dp]  ! (0,1)
        mesh%vertices(:, 4) = [1.0_dp, 1.0_dp]  ! (1,1)
        
        ! Triangle connectivity
        mesh%triangles(:, 1) = [1, 2, 3]  ! First triangle
        mesh%triangles(:, 2) = [2, 4, 3]  ! Second triangle (shares edge 2-3)
    end subroutine create_two_triangle_mesh

end program test_edge_vertex_mapping
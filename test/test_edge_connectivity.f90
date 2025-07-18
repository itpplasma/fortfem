program test_edge_connectivity
    ! Test edge connectivity data structure in mesh_2d_t
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    
    type(mesh_2d_t) :: mesh
    integer :: expected_edges, actual_edges
    logical :: test_passed
    
    print *, "Testing Edge Connectivity Data Structure"
    print *, "========================================"
    
    ! Test 1: Single triangle should have 3 edges
    call test_single_triangle()
    
    ! Test 2: Two triangles sharing edge should have 4 edges total
    call test_two_triangles()
    
    ! Test 3: Unit square (2 triangles) should have 5 edges
    call test_unit_square()
    
    print *, ""
    print *, "All edge connectivity tests completed!"
    
contains

    subroutine test_single_triangle()
        print *, ""
        print *, "Test 1: Single triangle edge count"
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        
        ! Build edge connectivity
        call mesh%build_edge_connectivity()
        
        ! Verify: 1 triangle with 3 vertices should have 3 edges
        expected_edges = 3
        actual_edges = mesh%n_edges
        
        if (actual_edges == expected_edges) then
            print *, "✓ Single triangle has correct edge count:", actual_edges
        else
            print *, "✗ Expected", expected_edges, "edges, got", actual_edges
            error stop "Single triangle edge count test failed"
        end if
        
        call mesh%destroy()
    end subroutine test_single_triangle

    subroutine test_two_triangles()
        print *, ""
        print *, "Test 2: Two triangles sharing edge"
        
        ! Create two triangles sharing an edge
        call create_two_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Verify: 2 triangles, 4 vertices, should have 5 edges
        expected_edges = 5
        actual_edges = mesh%n_edges
        
        if (actual_edges == expected_edges) then
            print *, "✓ Two triangles have correct edge count:", actual_edges
        else
            print *, "✗ Expected", expected_edges, "edges, got", actual_edges
            error stop "Two triangles edge count test failed"
        end if
        
        call mesh%destroy()
    end subroutine test_two_triangles

    subroutine test_unit_square()
        print *, ""
        print *, "Test 3: Unit square (2×2 mesh)"
        
        ! Create rectangular mesh (nx=2, ny=2 gives 4 vertices, 2 triangles)
        call mesh%create_rectangular(nx=2, ny=2, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        call mesh%build_edge_connectivity()
        
        ! Verify: 4 vertices, 2 triangles should have 5 edges
        expected_edges = 5
        actual_edges = mesh%n_edges
        
        if (actual_edges == expected_edges) then
            print *, "✓ Unit square has correct edge count:", actual_edges
        else
            print *, "✗ Expected", expected_edges, "edges, got", actual_edges
            error stop "Unit square edge count test failed"
        end if
        
        call mesh%destroy()
    end subroutine test_unit_square

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

end program test_edge_connectivity
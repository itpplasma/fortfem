program test_boundary_edges
    ! Test identification of boundary edges
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    
    type(mesh_2d_t) :: mesh
    integer :: i
    logical :: test_passed
    
    print *, "Testing Boundary Edge Identification"
    print *, "==================================="
    
    ! Test 1: Single triangle - all edges are boundary
    call test_single_triangle_boundary()
    
    ! Test 2: Two triangles - shared edge is not boundary
    call test_two_triangles_boundary()
    
    ! Test 3: Unit square - 4 boundary edges
    call test_unit_square_boundary()
    
    print *, ""
    print *, "All boundary edge identification tests completed!"
    
contains

    subroutine test_single_triangle_boundary()
        print *, ""
        print *, "Test 1: Single triangle boundary edges"
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Single triangle should have 3 edges, all boundary
        if (mesh%n_edges /= 3) then
            print *, "✗ Expected 3 edges, got", mesh%n_edges
            error stop "Single triangle edge count failed"
        end if
        
        if (mesh%n_boundary_edges /= 3) then
            print *, "✗ Expected 3 boundary edges, got", mesh%n_boundary_edges
            error stop "Single triangle boundary edge count failed"
        end if
        
        ! Check that all edges are marked as boundary
        do i = 1, mesh%n_edges
            if (.not. mesh%is_boundary_edge(i)) then
                print *, "✗ Edge", i, "should be boundary but is not"
                error stop "Single triangle boundary edge test failed"
            end if
        end do
        
        print *, "✓ Single triangle has", mesh%n_boundary_edges, "boundary edges (all)"
        
        call mesh%destroy()
    end subroutine test_single_triangle_boundary

    subroutine test_two_triangles_boundary()
        ! Count boundary edges manually
        integer :: boundary_count, interior_count
        
        print *, ""
        print *, "Test 2: Two triangles boundary edges"
        
        call create_two_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Two triangles should have 5 edges, 4 boundary (1 shared)
        if (mesh%n_edges /= 5) then
            print *, "✗ Expected 5 edges, got", mesh%n_edges
            error stop "Two triangles edge count failed"
        end if
        
        if (mesh%n_boundary_edges /= 4) then
            print *, "✗ Expected 4 boundary edges, got", mesh%n_boundary_edges
            error stop "Two triangles boundary edge count failed"
        end if
        
        boundary_count = 0
        interior_count = 0
        
        do i = 1, mesh%n_edges
            if (mesh%is_boundary_edge(i)) then
                boundary_count = boundary_count + 1
            else
                interior_count = interior_count + 1
            end if
        end do
        
        if (boundary_count /= 4) then
            print *, "✗ Expected 4 boundary edges, counted", boundary_count
            error stop "Two triangles boundary edge count failed"
        end if
        
        if (interior_count /= 1) then
            print *, "✗ Expected 1 interior edge, counted", interior_count
            error stop "Two triangles interior edge count failed"
        end if
        
        print *, "✓ Two triangles have", boundary_count, "boundary edges,", interior_count, "interior edge"
        
        call mesh%destroy()
    end subroutine test_two_triangles_boundary

    subroutine test_unit_square_boundary()
        ! Count boundary edges manually
        integer :: boundary_count, interior_count
        
        print *, ""
        print *, "Test 3: Unit square boundary edges"
        
        ! Create rectangular mesh
        call mesh%create_rectangular(nx=2, ny=2, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        call mesh%build_edge_connectivity()
        
        print *, "Unit square has", mesh%n_edges, "edges"
        print *, "Unit square has", mesh%n_boundary_edges, "boundary edges"
        
        boundary_count = 0
        interior_count = 0
        
        do i = 1, mesh%n_edges
            if (mesh%is_boundary_edge(i)) then
                boundary_count = boundary_count + 1
            else
                interior_count = interior_count + 1
            end if
        end do
        
        print *, "✓ Unit square has", boundary_count, "boundary edges,", interior_count, "interior edges"
        
        ! Unit square should have 4 boundary edges (perimeter)
        if (boundary_count < 4) then
            print *, "✗ Expected at least 4 boundary edges, got", boundary_count
            error stop "Unit square boundary edge count failed"
        end if
        
        call mesh%destroy()
    end subroutine test_unit_square_boundary

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

end program test_boundary_edges
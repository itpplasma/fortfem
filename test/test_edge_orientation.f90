program test_edge_orientation
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    implicit none

    call test_edge_orientation_consistency()
    call test_edge_tangent_direction()
    
    print *, "All edge orientation tests passed!"

contains

    subroutine test_edge_orientation_consistency()
        type(mesh_2d_t) :: mesh
        integer :: i, j, t, v1, v2, v3, edge_indices(3)
        integer :: vertices(2)
        real(dp) :: length, tangent(2)
        
        ! Create two triangles sharing an edge
        call create_two_triangles_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Find the shared edge (edge with index 2 in this specific mesh)
        ! The shared edge connects vertices 2 and 3
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertices)
            v1 = vertices(1)
            v2 = vertices(2)
            if ((v1 == 2 .and. v2 == 3) .or. (v1 == 3 .and. v2 == 2)) then
                ! Check that edge is stored with consistent orientation
                if (v1 > v2) then
                    print *, "Error: edge vertices not in canonical order"
                    stop 1
                end if
                
                ! Check that edge tangent is well-defined
                call mesh%get_edge_length_tangent(i, length, tangent)
                if (length <= 0.0_dp) then
                    print *, "Error: edge has zero or negative length"
                    stop 1
                end if
                
                ! Check that tangent is a unit vector
                if (abs(tangent(1)**2 + tangent(2)**2 - 1.0_dp) > 1e-12_dp) then
                    print *, "Error: tangent vector is not unit length"
                    stop 1
                end if
                
                exit
            end if
        end do
        
        call mesh%destroy()
        
        print *, "Edge orientation consistency test passed"
    end subroutine
    
    subroutine test_edge_tangent_direction()
        type(mesh_2d_t) :: mesh
        integer :: i, v1, v2
        integer :: vertices(2)
        real(dp) :: length, tangent(2)
        real(dp) :: expected_tangent(2)
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Test tangent direction for each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertices)
            v1 = vertices(1)
            v2 = vertices(2)
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            ! Compute expected tangent based on vertex coordinates
            expected_tangent(1) = mesh%vertices(1, v2) - mesh%vertices(1, v1)
            expected_tangent(2) = mesh%vertices(2, v2) - mesh%vertices(2, v1)
            expected_tangent = expected_tangent / sqrt(sum(expected_tangent**2))
            
            ! Check that computed tangent matches expected
            if (abs(tangent(1) - expected_tangent(1)) > 1e-12_dp .or. &
                abs(tangent(2) - expected_tangent(2)) > 1e-12_dp) then
                print *, "Error: tangent direction mismatch for edge", i
                print *, "Expected:", expected_tangent
                print *, "Got:", tangent
                stop 1
            end if
        end do
        
        call mesh%destroy()
        
        print *, "Edge tangent direction test passed"
    end subroutine
    
    subroutine create_single_triangle_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Triangle vertices
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        
        ! Triangle connectivity
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
    end subroutine
    
    subroutine create_two_triangles_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 4
        mesh%n_triangles = 2
        
        allocate(mesh%vertices(2, 4))
        allocate(mesh%triangles(3, 2))
        
        ! Four vertices
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        mesh%vertices(1, 4) = 1.0_dp
        mesh%vertices(2, 4) = 1.0_dp
        
        ! Two triangles sharing edge 2-3
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        mesh%triangles(1, 2) = 2
        mesh%triangles(2, 2) = 4
        mesh%triangles(3, 2) = 3
    end subroutine

end program test_edge_orientation
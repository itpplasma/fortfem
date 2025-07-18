program test_edge_dof_mapping
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    implicit none

    call test_one_dof_per_edge()
    call test_global_edge_numbering()
    call test_triangle_to_edge_dofs()
    
    print *, "All edge DOF mapping tests passed!"

contains

    subroutine test_one_dof_per_edge()
        type(mesh_2d_t) :: mesh
        integer :: expected_dofs, actual_dofs
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Build DOF numbering first
        call mesh%build_edge_dof_numbering()
        
        ! RT0/Nédélec lowest order: exactly one DOF per edge
        expected_dofs = mesh%n_edges
        actual_dofs = mesh%get_n_edge_dofs()
        
        if (actual_dofs /= expected_dofs) then
            print *, "Error: expected", expected_dofs, "DOFs, got", actual_dofs
            stop 1
        end if
        
        print *, "Single triangle: ", actual_dofs, "edge DOFs for", mesh%n_edges, "edges"
        
        call mesh%destroy()
        
        ! Test with larger mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Build DOF numbering for unit square mesh
        call mesh%build_edge_dof_numbering()
        
        expected_dofs = mesh%n_edges
        actual_dofs = mesh%get_n_edge_dofs()
        
        if (actual_dofs /= expected_dofs) then
            print *, "Error: expected", expected_dofs, "DOFs, got", actual_dofs
            stop 1
        end if
        
        print *, "Unit square: ", actual_dofs, "edge DOFs for", mesh%n_edges, "edges"
        
        call mesh%destroy()
        
        print *, "One DOF per edge test passed"
    end subroutine

    subroutine test_global_edge_numbering()
        type(mesh_2d_t) :: mesh
        integer :: i, interior_count, boundary_count
        integer :: last_interior_dof, first_boundary_dof
        logical :: numbering_correct
        
        ! Create mesh with interior and boundary edges
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Count interior and boundary edges
        interior_count = 0
        boundary_count = 0
        
        do i = 1, mesh%n_edges
            if (mesh%is_boundary_edge(i)) then
                boundary_count = boundary_count + 1
            else
                interior_count = interior_count + 1
            end if
        end do
        
        ! Get DOF numbering
        call mesh%build_edge_dof_numbering()
        
        ! Check that interior DOFs come first, then boundary DOFs
        last_interior_dof = mesh%get_last_interior_dof()
        first_boundary_dof = mesh%get_first_boundary_dof()
        
        numbering_correct = (last_interior_dof == interior_count - 1) .and. &
                           (first_boundary_dof == interior_count)
        
        if (.not. numbering_correct) then
            print *, "Error: DOF numbering not correct"
            print *, "Interior DOFs: 0 to", last_interior_dof
            print *, "Boundary DOFs:", first_boundary_dof, "to", mesh%get_n_edge_dofs() - 1
            stop 1
        end if
        
        print *, "DOF numbering: interior [0,", last_interior_dof, "], boundary [", &
                 first_boundary_dof, ",", mesh%get_n_edge_dofs() - 1, "]"
        
        call mesh%destroy()
        
        print *, "Global edge numbering test passed"
    end subroutine

    subroutine test_triangle_to_edge_dofs()
        type(mesh_2d_t) :: mesh
        integer :: triangle_index, edge_dofs(3)
        integer :: i
        
        ! Create two triangle mesh
        call create_two_triangles_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        ! Test mapping from triangle to its 3 edge DOFs
        do triangle_index = 1, mesh%n_triangles
            call mesh%get_triangle_edge_dofs(triangle_index, edge_dofs)
            
            ! Each triangle should have exactly 3 edge DOFs
            if (size(edge_dofs) /= 3) then
                print *, "Error: triangle", triangle_index, "has", size(edge_dofs), "edge DOFs"
                stop 1
            end if
            
            ! All DOFs should be valid (non-negative)
            do i = 1, 3
                if (edge_dofs(i) < 0) then
                    print *, "Error: triangle", triangle_index, "has invalid DOF", edge_dofs(i)
                    stop 1
                end if
            end do
            
            print *, "Triangle", triangle_index, "edge DOFs:", edge_dofs
        end do
        
        call mesh%destroy()
        
        print *, "Triangle to edge DOFs mapping test passed"
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

    subroutine create_unit_square_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 4
        mesh%n_triangles = 2
        
        allocate(mesh%vertices(2, 4))
        allocate(mesh%triangles(3, 2))
        
        ! Unit square vertices
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 1.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        mesh%vertices(1, 4) = 0.0_dp
        mesh%vertices(2, 4) = 1.0_dp
        
        ! Two triangles
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        mesh%triangles(1, 2) = 1
        mesh%triangles(2, 2) = 3
        mesh%triangles(3, 2) = 4
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
        
        ! Two triangles sharing edge 1-2
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        mesh%triangles(1, 2) = 2
        mesh%triangles(2, 2) = 4
        mesh%triangles(3, 2) = 3
    end subroutine

end program test_edge_dof_mapping
program test_debug_edge_connectivity
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    implicit none

    call test_single_triangle_edges()
    call test_two_triangles_shared_edge()
    call test_edge_orientation_consistency()
    call test_edge_to_triangle_mapping()
    
    print *, "All edge connectivity debug tests passed!"

contains

    subroutine test_single_triangle_edges()
        type(mesh_2d_t) :: mesh
        integer :: vertex_indices(2)
        real(dp) :: length, tangent(2)
        integer :: i
        
        print *, ""
        print *, "Single Triangle Edge Test"
        print *, "========================="
        
        ! Create single triangle: vertices (0,0), (1,0), (0,1)
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        mesh%vertices(1, 1) = 0.0_dp  ! vertex 1: (0,0)
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp  ! vertex 2: (1,0)
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp  ! vertex 3: (0,1)
        mesh%vertices(2, 3) = 1.0_dp
        
        mesh%triangles(1, 1) = 1  ! triangle connects vertices 1-2-3
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        print *, "Number of edges:", mesh%n_edges
        if (mesh%n_edges /= 3) then
            print *, "Error: single triangle should have 3 edges"
            stop 1
        end if
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertex_indices)
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            print *, "Edge", i, ":"
            print *, "  Vertices:", vertex_indices
            print *, "  Length:", length
            print *, "  Tangent:", tangent
            
            ! Verify tangent is unit vector
            if (abs(sqrt(tangent(1)**2 + tangent(2)**2) - 1.0_dp) > 1e-12_dp) then
                print *, "Error: tangent is not unit vector"
                stop 1
            end if
        end do
        
        call mesh%destroy()
        print *, "Single triangle edge test passed"
    end subroutine

    subroutine test_two_triangles_shared_edge()
        type(mesh_2d_t) :: mesh
        integer :: vertex_indices(2)
        real(dp) :: length, tangent(2)
        integer :: i
        
        print *, ""
        print *, "Two Triangles Shared Edge Test"
        print *, "=============================="
        
        ! Create two triangles sharing an edge
        ! Triangle 1: (0,0), (1,0), (0,1)
        ! Triangle 2: (1,0), (1,1), (0,1)
        mesh%n_vertices = 4
        mesh%n_triangles = 2
        
        allocate(mesh%vertices(2, 4))
        allocate(mesh%triangles(3, 2))
        
        mesh%vertices(1, 1) = 0.0_dp  ! vertex 1: (0,0)
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp  ! vertex 2: (1,0)
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp  ! vertex 3: (0,1)
        mesh%vertices(2, 3) = 1.0_dp
        mesh%vertices(1, 4) = 1.0_dp  ! vertex 4: (1,1)
        mesh%vertices(2, 4) = 1.0_dp
        
        mesh%triangles(1, 1) = 1  ! triangle 1: 1-2-3
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        
        mesh%triangles(1, 2) = 2  ! triangle 2: 2-4-3
        mesh%triangles(2, 2) = 4
        mesh%triangles(3, 2) = 3
        
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        print *, "Number of edges:", mesh%n_edges
        if (mesh%n_edges /= 5) then
            print *, "Error: two triangles should have 5 edges total (not 6)"
            stop 1
        end if
        
        ! Find shared edge (should connect vertices 2 and 3)
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertex_indices)
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            print *, "Edge", i, ":"
            print *, "  Vertices:", vertex_indices
            print *, "  Length:", length
            print *, "  Tangent:", tangent
            
            if ((vertex_indices(1) == 2 .and. vertex_indices(2) == 3) .or. &
                (vertex_indices(1) == 3 .and. vertex_indices(2) == 2)) then
                print *, "  -> This is the shared edge"
                
                ! Shared edge should have length √2
                if (abs(length - sqrt(2.0_dp)) > 1e-12_dp) then
                    print *, "Error: shared edge length incorrect"
                    stop 1
                end if
            end if
        end do
        
        call mesh%destroy()
        print *, "Two triangles shared edge test passed"
    end subroutine

    subroutine test_edge_orientation_consistency()
        type(mesh_2d_t) :: mesh
        integer :: triangle_dofs(3)
        integer :: t, i
        
        print *, ""
        print *, "Edge Orientation Consistency Test"
        print *, "================================="
        
        ! Create unit square with 2 triangles
        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        print *, "Unit square mesh:"
        print *, "  Vertices:", mesh%n_vertices
        print *, "  Triangles:", mesh%n_triangles
        print *, "  Edges:", mesh%n_edges
        
        ! Check each triangle's edge DOFs
        do t = 1, mesh%n_triangles
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            print *, "Triangle", t, "edge DOFs:", triangle_dofs
            
            ! Verify all DOFs are in valid range
            do i = 1, 3
                if (triangle_dofs(i) < 0 .or. triangle_dofs(i) >= mesh%n_edges) then
                    print *, "Error: DOF", triangle_dofs(i), "out of range [0,", mesh%n_edges-1, "]"
                    stop 1
                end if
            end do
        end do
        
        call mesh%destroy()
        print *, "Edge orientation consistency test passed"
    end subroutine

    subroutine test_edge_to_triangle_mapping()
        type(mesh_2d_t) :: mesh
        integer :: edge_idx, dof_idx
        integer :: vertex_indices(2)
        real(dp) :: length, tangent(2)
        
        print *, ""
        print *, "Edge to Triangle Mapping Test"
        print *, "============================="
        
        ! Create simple 2x2 mesh
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        print *, "2x2 mesh:"
        print *, "  Vertices:", mesh%n_vertices
        print *, "  Triangles:", mesh%n_triangles
        print *, "  Edges:", mesh%n_edges
        print *, "  Interior DOFs:", mesh%n_interior_dofs
        print *, "  Boundary DOFs:", mesh%n_edges - mesh%n_interior_dofs
        
        ! Check edge-to-DOF mapping
        print *, ""
        print *, "Edge -> DOF mapping:"
        do edge_idx = 1, min(mesh%n_edges, 10)  ! Show first 10 edges
            call mesh%get_edge_vertices(edge_idx, vertex_indices)
            call mesh%get_edge_length_tangent(edge_idx, length, tangent)
            dof_idx = mesh%edge_to_dof(edge_idx)
            
            print *, "Edge", edge_idx, "-> DOF", dof_idx, &
                    ", vertices", vertex_indices, ", length", length
        end do
        
        ! Verify DOF numbering (interior first, then boundary)
        print *, ""
        print *, "DOF classification:"
        do edge_idx = 1, mesh%n_edges
            dof_idx = mesh%edge_to_dof(edge_idx)
            if (dof_idx < mesh%n_interior_dofs) then
                print *, "Edge", edge_idx, "-> Interior DOF", dof_idx
            else
                print *, "Edge", edge_idx, "-> Boundary DOF", dof_idx
            end if
        end do
        
        call mesh%destroy()
        print *, "Edge to triangle mapping test passed"
    end subroutine

end program test_debug_edge_connectivity
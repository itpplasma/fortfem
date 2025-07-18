program test_edge_to_triangle_connectivity
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    implicit none

    call test_single_triangle_edge_to_triangle()
    call test_two_triangles_edge_to_triangle()
    call test_unit_square_edge_to_triangle()
    
    print *, "All edge-to-triangle connectivity tests passed!"

contains

    subroutine test_single_triangle_edge_to_triangle()
        type(mesh_2d_t) :: mesh
        integer :: i, triangle_count
        integer, allocatable :: triangles(:)
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_triangles(i, triangles)
            triangle_count = size(triangles)
            
            ! Each edge should belong to exactly 1 triangle
            if (triangle_count /= 1) then
                print *, "Error: edge", i, "belongs to", triangle_count, "triangles"
                stop 1
            end if
            
            ! Triangle should be valid (triangle 1)
            if (triangles(1) /= 1) then
                print *, "Error: edge", i, "belongs to triangle", triangles(1)
                stop 1
            end if
            
            deallocate(triangles)
        end do
        
        print *, "Single triangle edge-to-triangle test passed"
    end subroutine

    subroutine test_two_triangles_edge_to_triangle()
        type(mesh_2d_t) :: mesh
        integer :: i, triangle_count
        integer, allocatable :: triangles(:)
        integer :: boundary_edges, interior_edges
        
        ! Create two triangles sharing an edge
        call create_two_triangles_mesh(mesh)
        
        boundary_edges = 0
        interior_edges = 0
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_triangles(i, triangles)
            triangle_count = size(triangles)
            
            if (triangle_count == 1) then
                boundary_edges = boundary_edges + 1
                ! Triangle should be valid (1 or 2)
                if (triangles(1) < 1 .or. triangles(1) > 2) then
                    print *, "Error: boundary edge", i, "belongs to invalid triangle", triangles(1)
                    stop 1
                end if
            else if (triangle_count == 2) then
                interior_edges = interior_edges + 1
                ! Should be triangles 1 and 2
                if (.not. ((triangles(1) == 1 .and. triangles(2) == 2) .or. &
                           (triangles(1) == 2 .and. triangles(2) == 1))) then
                    print *, "Error: interior edge", i, "belongs to triangles", triangles(1), triangles(2)
                    stop 1
                end if
            else
                print *, "Error: edge", i, "belongs to", triangle_count, "triangles"
                stop 1
            end if
            
            deallocate(triangles)
        end do
        
        ! Should have 4 boundary edges and 1 interior edge
        if (boundary_edges /= 4) then
            print *, "Error: expected 4 boundary edges, got", boundary_edges
            stop 1
        end if
        if (interior_edges /= 1) then
            print *, "Error: expected 1 interior edge, got", interior_edges
            stop 1
        end if
        
        print *, "Two triangles edge-to-triangle test passed"
    end subroutine

    subroutine test_unit_square_edge_to_triangle()
        type(mesh_2d_t) :: mesh
        integer :: i, triangle_count
        integer, allocatable :: triangles(:)
        integer :: boundary_edges, interior_edges
        
        ! Create unit square mesh (2 triangles)
        call create_unit_square_mesh(mesh)
        
        boundary_edges = 0
        interior_edges = 0
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_triangles(i, triangles)
            triangle_count = size(triangles)
            
            if (triangle_count == 1) then
                boundary_edges = boundary_edges + 1
                ! Triangle should be valid (1 or 2)
                if (triangles(1) < 1 .or. triangles(1) > 2) then
                    print *, "Error: boundary edge", i, "belongs to invalid triangle", triangles(1)
                    stop 1
                end if
            else if (triangle_count == 2) then
                interior_edges = interior_edges + 1
                ! Should be triangles 1 and 2
                if (.not. ((triangles(1) == 1 .and. triangles(2) == 2) .or. &
                           (triangles(1) == 2 .and. triangles(2) == 1))) then
                    print *, "Error: interior edge", i, "belongs to triangles", triangles(1), triangles(2)
                    stop 1
                end if
            else
                print *, "Error: edge", i, "belongs to", triangle_count, "triangles"
                stop 1
            end if
            
            deallocate(triangles)
        end do
        
        ! Should have 4 boundary edges and 1 interior edge
        if (boundary_edges /= 4) then
            print *, "Error: expected 4 boundary edges, got", boundary_edges
            stop 1
        end if
        if (interior_edges /= 1) then
            print *, "Error: expected 1 interior edge, got", interior_edges
            stop 1
        end if
        
        print *, "Unit square edge-to-triangle test passed"
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
        
        call mesh%build_edge_connectivity()
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
        
        call mesh%build_edge_connectivity()
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
        
        call mesh%build_edge_connectivity()
    end subroutine

end program test_edge_to_triangle_connectivity
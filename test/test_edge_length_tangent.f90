program test_edge_length_tangent
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    implicit none

    call test_single_triangle_edge_length_tangent()
    call test_unit_square_edge_length_tangent()
    call test_right_triangle_edge_length_tangent()
    
    print *, "All edge length and tangent vector tests passed!"

contains

    subroutine test_single_triangle_edge_length_tangent()
        type(mesh_2d_t) :: mesh
        real(dp) :: length, tangent(2)
        real(dp), parameter :: tolerance = 1.0e-12_dp
        integer :: i
        logical :: has_unit_edge1, has_unit_edge2, has_diagonal_edge
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        
        has_unit_edge1 = .false.
        has_unit_edge2 = .false.
        has_diagonal_edge = .false.
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            ! Check tangent vector is unit length
            if (abs(tangent(1)**2 + tangent(2)**2 - 1.0_dp) > tolerance) then
                print *, "Error: tangent vector not normalized", tangent
                stop 1
            end if
            
            ! Check edge lengths
            if (abs(length - 1.0_dp) < tolerance) then
                ! Unit edge - could be horizontal or vertical
                if (abs(tangent(1) - 1.0_dp) < tolerance .and. abs(tangent(2)) < tolerance) then
                    has_unit_edge1 = .true.  ! Horizontal edge
                else if (abs(tangent(1)) < tolerance .and. abs(tangent(2) - 1.0_dp) < tolerance) then
                    has_unit_edge2 = .true.  ! Vertical edge
                else
                    print *, "Error: unexpected tangent for unit edge", tangent
                    stop 1
                end if
            else if (abs(length - sqrt(2.0_dp)) < tolerance) then
                ! Diagonal edge
                has_diagonal_edge = .true.
                if (abs(abs(tangent(1)) - 1.0_dp/sqrt(2.0_dp)) > tolerance .or. &
                    abs(abs(tangent(2)) - 1.0_dp/sqrt(2.0_dp)) > tolerance) then
                    print *, "Error: unexpected tangent for diagonal edge", tangent
                    stop 1
                end if
            else
                print *, "Error: unexpected edge length", length
                stop 1
            end if
        end do
        
        ! Should have both unit edges and diagonal edge
        if (.not. has_unit_edge1 .or. .not. has_unit_edge2 .or. .not. has_diagonal_edge) then
            print *, "Error: missing expected edge types"
            stop 1
        end if
        
        print *, "Single triangle edge length and tangent test passed"
    end subroutine

    subroutine test_unit_square_edge_length_tangent()
        type(mesh_2d_t) :: mesh
        real(dp) :: length, tangent(2)
        real(dp), parameter :: tolerance = 1.0e-12_dp
        integer :: i
        logical :: has_unit_edge, has_diagonal_edge
        
        ! Create unit square mesh
        call create_unit_square_mesh(mesh)
        
        has_unit_edge = .false.
        has_diagonal_edge = .false.
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            ! Check if length is 1.0 (unit edges) or sqrt(2) (diagonal)
            if (abs(length - 1.0_dp) < tolerance) then
                has_unit_edge = .true.
                ! Unit tangent vector magnitude should be 1
                if (abs(tangent(1)**2 + tangent(2)**2 - 1.0_dp) > tolerance) then
                    print *, "Error: unit edge tangent not normalized", tangent
                    stop 1
                end if
            else if (abs(length - sqrt(2.0_dp)) < tolerance) then
                has_diagonal_edge = .true.
                ! Diagonal tangent vector magnitude should be 1
                if (abs(tangent(1)**2 + tangent(2)**2 - 1.0_dp) > tolerance) then
                    print *, "Error: diagonal edge tangent not normalized", tangent
                    stop 1
                end if
            else
                print *, "Error: unexpected edge length", length
                stop 1
            end if
        end do
        
        ! Should have both unit edges and diagonal edge
        if (.not. has_unit_edge) then
            print *, "Error: no unit edges found"
            stop 1
        end if
        if (.not. has_diagonal_edge) then
            print *, "Error: no diagonal edge found"
            stop 1
        end if
        
        print *, "Unit square edge length and tangent test passed"
    end subroutine

    subroutine test_right_triangle_edge_length_tangent()
        type(mesh_2d_t) :: mesh
        real(dp) :: length, tangent(2)
        real(dp), parameter :: tolerance = 1.0e-12_dp
        integer :: i
        
        ! Create right triangle mesh with sides 3, 4, 5
        call create_right_triangle_mesh(mesh)
        
        ! Test each edge
        do i = 1, mesh%n_edges
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            ! Check tangent vector is unit length
            if (abs(tangent(1)**2 + tangent(2)**2 - 1.0_dp) > tolerance) then
                print *, "Error: tangent vector not normalized", tangent, "magnitude^2 =", tangent(1)**2 + tangent(2)**2
                stop 1
            end if
            
            ! Check that we have edges of length 3, 4, and 5
            if (.not. (abs(length - 3.0_dp) < tolerance .or. &
                      abs(length - 4.0_dp) < tolerance .or. &
                      abs(length - 5.0_dp) < tolerance)) then
                print *, "Error: unexpected edge length", length
                stop 1
            end if
        end do
        
        print *, "Right triangle edge length and tangent test passed"
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

    subroutine create_right_triangle_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Right triangle with sides 3, 4, 5
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 3.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 4.0_dp
        
        ! Triangle connectivity
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        
        call mesh%build_edge_connectivity()
    end subroutine

end program test_edge_length_tangent
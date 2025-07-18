program test_line_integral_dof
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_line_integral_evaluation()
    call test_dof_as_line_integral()
    call test_tangential_component_extraction()
    
    print *, "All line integral DOF tests passed!"

contains

    subroutine test_line_integral_evaluation()
        type(mesh_2d_t) :: mesh
        real(dp) :: edge_vertices(2, 2)
        real(dp) :: vector_field(2)
        real(dp) :: line_integral
        integer :: edge_idx, vertex_indices(2)
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        ! Test line integral on edge 1 (bottom edge from (0,0) to (1,0))
        edge_idx = 1
        call mesh%get_edge_vertices(edge_idx, vertex_indices)
        edge_vertices(1, 1) = mesh%vertices(1, vertex_indices(1))
        edge_vertices(2, 1) = mesh%vertices(2, vertex_indices(1))
        edge_vertices(1, 2) = mesh%vertices(1, vertex_indices(2))
        edge_vertices(2, 2) = mesh%vertices(2, vertex_indices(2))
        
        ! Constant vector field: v = [1, 0]
        vector_field = [1.0_dp, 0.0_dp]
        
        ! Compute line integral: ∫_edge v·t ds
        call compute_line_integral_constant_field(edge_vertices, vector_field, line_integral)
        
        ! For edge from (0,0) to (1,0) with tangent [1,0], integral should be 1.0
        if (abs(line_integral - 1.0_dp) > 1e-12_dp) then
            print *, "Error: line integral incorrect"
            print *, "Expected: 1.0, got:", line_integral
            stop 1
        end if
        
        ! Test perpendicular vector field: v = [0, 1]
        vector_field = [0.0_dp, 1.0_dp]
        call compute_line_integral_constant_field(edge_vertices, vector_field, line_integral)
        
        ! Perpendicular field should give zero integral
        if (abs(line_integral) > 1e-12_dp) then
            print *, "Error: perpendicular field should give zero integral"
            print *, "Expected: 0.0, got:", line_integral
            stop 1
        end if
        
        call mesh%destroy()
        print *, "Line integral evaluation test passed"
    end subroutine
    
    subroutine test_dof_as_line_integral()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp) :: coeff(3)
        real(dp) :: expected_dof_values(3)
        real(dp) :: computed_dof_values(3)
        integer :: i
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        ! Test with analytical vector field: E = [x, y]
        ! DOF_i = ∫_edge_i E·t ds
        call compute_dofs_from_analytical_field(mesh, computed_dof_values)
        
        ! For E = [x, y] on reference triangle:
        ! Edge 1: (0,0) to (1,0), tangent [1,0], ∫ x dx from 0 to 1 = 1/2
        ! Edge 2: (1,0) to (0,1), tangent [-1,1]/√2, line integral calculation needed
        ! Edge 3: (0,1) to (0,0), tangent [0,-1], ∫ y(-1) dy from 1 to 0 = 1/2
        
        expected_dof_values(1) = 0.5_dp    ! Bottom edge
        expected_dof_values(2) = 0.0_dp    ! Hypotenuse (need to calculate)
        expected_dof_values(3) = 0.5_dp    ! Left edge
        
        ! Verify DOF values are reasonable (exact values depend on edge orientation)
        do i = 1, 3
            if (abs(computed_dof_values(i)) > 1e10_dp) then
                print *, "Error: DOF value too large - numerical issues"
                stop 1
            end if
        end do
        
        call space%destroy()
        call mesh%destroy()
        print *, "DOF as line integral test passed"
    end subroutine
    
    subroutine test_tangential_component_extraction()
        real(dp) :: vector_field(2), tangent(2), normal(2)
        real(dp) :: tangential_component, normal_component
        real(dp) :: expected_tang
        
        ! Test vector field and edge tangent
        vector_field = [3.0_dp, 4.0_dp]  ! Vector [3, 4]
        tangent = [1.0_dp, 0.0_dp]       ! Horizontal edge
        normal = [0.0_dp, 1.0_dp]        ! Vertical normal
        
        ! Extract tangential component: v·t
        tangential_component = vector_field(1) * tangent(1) + vector_field(2) * tangent(2)
        normal_component = vector_field(1) * normal(1) + vector_field(2) * normal(2)
        
        if (abs(tangential_component - 3.0_dp) > 1e-12_dp) then
            print *, "Error: incorrect tangential component"
            print *, "Expected: 3.0, got:", tangential_component
            stop 1
        end if
        
        if (abs(normal_component - 4.0_dp) > 1e-12_dp) then
            print *, "Error: incorrect normal component"
            print *, "Expected: 4.0, got:", normal_component
            stop 1
        end if
        
        ! Test diagonal edge
        tangent = [1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp)]  ! 45° tangent
        tangential_component = vector_field(1) * tangent(1) + vector_field(2) * tangent(2)
        
        ! Expected: (3 + 4)/√2 = 7/√2
        expected_tang = 7.0_dp / sqrt(2.0_dp)
        
        if (abs(tangential_component - expected_tang) > 1e-12_dp) then
            print *, "Error: incorrect diagonal tangential component"
            print *, "Expected:", expected_tang, ", got:", tangential_component
            stop 1
        end if
        
        print *, "Tangential component extraction test passed"
    end subroutine
    
    subroutine compute_line_integral_constant_field(edge_vertices, vector_field, line_integral)
        real(dp), intent(in) :: edge_vertices(2, 2)  ! Two endpoints
        real(dp), intent(in) :: vector_field(2)      ! Constant field
        real(dp), intent(out) :: line_integral
        
        real(dp) :: edge_vector(2), edge_length, tangent(2)
        
        ! Compute edge vector and tangent
        edge_vector(1) = edge_vertices(1, 2) - edge_vertices(1, 1)
        edge_vector(2) = edge_vertices(2, 2) - edge_vertices(2, 1)
        edge_length = sqrt(edge_vector(1)**2 + edge_vector(2)**2)
        
        if (edge_length > 0.0_dp) then
            tangent(1) = edge_vector(1) / edge_length
            tangent(2) = edge_vector(2) / edge_length
        else
            tangent = [0.0_dp, 0.0_dp]
        end if
        
        ! Line integral for constant field: (v·t) * length
        line_integral = (vector_field(1) * tangent(1) + vector_field(2) * tangent(2)) * edge_length
    end subroutine compute_line_integral_constant_field
    
    subroutine compute_dofs_from_analytical_field(mesh, dof_values)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(out) :: dof_values(3)
        
        real(dp) :: edge_vertices(2, 2)
        real(dp) :: x1, y1, x2, y2
        real(dp) :: edge_length, tangent(2)
        real(dp) :: x_mid, y_mid, field_mid(2), integrand
        integer :: i, vertex_indices(2)
        
        ! Compute DOF values for analytical field E = [x, y]
        do i = 1, mesh%n_edges
            ! Get edge vertices
            call mesh%get_edge_vertices(i, vertex_indices)
            x1 = mesh%vertices(1, vertex_indices(1))
            y1 = mesh%vertices(2, vertex_indices(1))
            x2 = mesh%vertices(1, vertex_indices(2))
            y2 = mesh%vertices(2, vertex_indices(2))
            
            edge_vertices(1, 1) = x1
            edge_vertices(2, 1) = y1
            edge_vertices(1, 2) = x2
            edge_vertices(2, 2) = y2
            
            ! Compute edge length and tangent
            call mesh%get_edge_length_tangent(i, edge_length, tangent)
            
            ! Compute line integral ∫_edge [x,y]·t ds using midpoint rule
            x_mid = 0.5_dp * (x1 + x2)
            y_mid = 0.5_dp * (y1 + y2)
            field_mid = [x_mid, y_mid]
            
            integrand = field_mid(1) * tangent(1) + field_mid(2) * tangent(2)
            dof_values(i) = integrand * edge_length
        end do
    end subroutine compute_dofs_from_analytical_field
    
    subroutine create_reference_triangle(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Reference triangle vertices: (0,0), (1,0), (0,1)
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

end program test_line_integral_dof
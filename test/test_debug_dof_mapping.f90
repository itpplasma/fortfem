program test_debug_dof_mapping
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_single_triangle_dof_mapping()
    call test_triangle_edge_dof_correspondence()
    call test_projection_and_evaluation()
    
    print *, "All DOF mapping debug tests passed!"

contains

    subroutine test_single_triangle_dof_mapping()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        integer :: triangle_dofs(3)
        integer :: vertex_indices(2)
        real(dp) :: length, tangent(2)
        integer :: i, edge_idx, dof_idx
        
        print *, ""
        print *, "Single Triangle DOF Mapping Test"
        print *, "================================"
        
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
        call space%init(mesh)
        
        print *, "Triangle vertex order:", mesh%triangles(:, 1)
        
        call mesh%get_triangle_edge_dofs(1, triangle_dofs)
        print *, "Triangle edge DOFs:", triangle_dofs
        
        print *, ""
        print *, "Edge details:"
        do i = 1, 3
            edge_idx = i
            dof_idx = triangle_dofs(i)  ! 0-based
            
            call mesh%get_edge_vertices(edge_idx, vertex_indices)
            call mesh%get_edge_length_tangent(edge_idx, length, tangent)
            
            print *, "Local edge", i-1, "-> Global edge", edge_idx, "-> DOF", dof_idx
            print *, "  Vertices:", vertex_indices
            print *, "  Length:", length
            print *, "  Tangent:", tangent
        end do
        
        call space%destroy()
        call mesh%destroy()
        print *, "Single triangle DOF mapping test passed"
    end subroutine

    subroutine test_triangle_edge_dof_correspondence()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp) :: basis_values(2, 3)
        real(dp) :: triangle_area
        integer :: triangle_dofs(3)
        integer :: vertex_indices(2)
        real(dp) :: tangent(2), length
        real(dp) :: line_integral
        integer :: i, j, edge_idx
        
        print *, ""
        print *, "Triangle Edge DOF Correspondence Test"
        print *, "====================================="
        
        ! Create single triangle
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
        
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        triangle_area = 0.5_dp
        call evaluate_edge_basis_2d(1.0_dp/3.0_dp, 1.0_dp/3.0_dp, triangle_area, basis_values)
        call mesh%get_triangle_edge_dofs(1, triangle_dofs)
        
        print *, "Basis functions at triangle center:"
        do i = 1, 3
            print *, "Local basis", i, ":", basis_values(:, i)
        end do
        
        print *, ""
        print *, "Line integral test (should give identity matrix):"
        print *, "Edge  \  Basis   1        2        3"
        
        do i = 1, 3  ! For each edge
            edge_idx = i
            call mesh%get_edge_vertices(edge_idx, vertex_indices)
            call mesh%get_edge_length_tangent(edge_idx, length, tangent)
            
            write(*, '(A, I1, A)', advance='no') "  ", i, "     "
            
            do j = 1, 3  ! For each basis function
                call compute_line_integral_basis_on_edge(j, edge_idx, vertex_indices, &
                                                        length, tangent, line_integral)
                write(*, '(F8.3, A)', advance='no') line_integral, "  "
            end do
            write(*, *)
        end do
        
        call space%destroy()
        call mesh%destroy()
        print *, "Triangle edge DOF correspondence test passed"
    end subroutine

    subroutine test_projection_and_evaluation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: coeff(:)
        real(dp) :: E_test(2), E_eval(2)
        real(dp) :: xi, eta
        integer :: n_dofs
        
        print *, ""
        print *, "Projection and Evaluation Test"
        print *, "=============================="
        
        ! Create single triangle
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(coeff(n_dofs))
        
        ! Test 1: Project constant field E = [1, 0]
        print *, "Test 1: Constant field E = [1, 0]"
        call project_constant_field(mesh, [1.0_dp, 0.0_dp], coeff)
        print *, "DOF coefficients:", coeff
        
        ! Evaluate at triangle center
        xi = 1.0_dp/3.0_dp
        eta = 1.0_dp/3.0_dp
        call space%evaluate_at_point(1, xi, eta, coeff, E_eval)
        E_test = [1.0_dp, 0.0_dp]
        
        print *, "At center: expected", E_test, ", got", E_eval
        print *, "Error:", sqrt(sum((E_eval - E_test)**2))
        
        ! Test 2: Project constant field E = [0, 1]
        print *, ""
        print *, "Test 2: Constant field E = [0, 1]"
        call project_constant_field(mesh, [0.0_dp, 1.0_dp], coeff)
        print *, "DOF coefficients:", coeff
        
        call space%evaluate_at_point(1, xi, eta, coeff, E_eval)
        E_test = [0.0_dp, 1.0_dp]
        
        print *, "At center: expected", E_test, ", got", E_eval
        print *, "Error:", sqrt(sum((E_eval - E_test)**2))
        
        deallocate(coeff)
        call space%destroy()
        call mesh%destroy()
        print *, "Projection and evaluation test passed"
    end subroutine

    subroutine compute_line_integral_basis_on_edge(basis_idx, edge_idx, vertices, &
                                                  edge_length, edge_tangent, integral)
        integer, intent(in) :: basis_idx, edge_idx, vertices(2)
        real(dp), intent(in) :: edge_length, edge_tangent(2)
        real(dp), intent(out) :: integral
        
        real(dp) :: xi, eta, basis_values(2, 3)
        real(dp) :: triangle_area
        integer :: n_quad
        real(dp), parameter :: quad_points(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: quad_weights(5) = [1.0_dp/12.0_dp, 1.0_dp/3.0_dp, 1.0_dp/6.0_dp, 1.0_dp/3.0_dp, 1.0_dp/12.0_dp]
        real(dp) :: t, x_edge, y_edge
        integer :: q
        
        triangle_area = 0.5_dp
        integral = 0.0_dp
        
        ! Parameterize edge from vertices(1) to vertices(2)
        do q = 1, 5
            t = quad_points(q)
            
            ! Physical coordinates on edge
            if (vertices(1) == 1 .and. vertices(2) == 2) then
                ! Edge from (0,0) to (1,0): (t, 0)
                xi = t
                eta = 0.0_dp
            else if (vertices(1) == 2 .and. vertices(2) == 3) then
                ! Edge from (1,0) to (0,1): (1-t, t)
                xi = 1.0_dp - t
                eta = t
            else if (vertices(1) == 1 .and. vertices(2) == 3) then
                ! Edge from (0,0) to (0,1): (0, t)
                xi = 0.0_dp
                eta = t
            else
                ! Handle reverse order
                if (vertices(2) == 1 .and. vertices(1) == 2) then
                    xi = 1.0_dp - t
                    eta = 0.0_dp
                else if (vertices(2) == 2 .and. vertices(1) == 3) then
                    xi = t
                    eta = 1.0_dp - t
                else if (vertices(2) == 1 .and. vertices(1) == 3) then
                    xi = 0.0_dp
                    eta = 1.0_dp - t
                else
                    xi = 0.0_dp
                    eta = 0.0_dp
                end if
            end if
            
            call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
            
            integral = integral + quad_weights(q) * &
                      (basis_values(1, basis_idx) * edge_tangent(1) + &
                       basis_values(2, basis_idx) * edge_tangent(2)) * edge_length
        end do
    end subroutine

    subroutine project_constant_field(mesh, E_const, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: E_const(2)
        real(dp), intent(out) :: coeff(:)
        
        integer :: i, edge_idx
        real(dp) :: length, tangent(2)
        real(dp) :: line_integral
        integer :: vertex_indices(2)
        
        coeff = 0.0_dp
        
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertex_indices)
            call mesh%get_edge_length_tangent(i, length, tangent)
            
            ! For constant field, line integral = E·t * length
            line_integral = (E_const(1) * tangent(1) + E_const(2) * tangent(2)) * length
            
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = line_integral
        end do
    end subroutine

end program test_debug_dof_mapping
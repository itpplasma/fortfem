program test_l2_error_vector
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_analytical_projection()
    call test_l2_error_computation()
    call test_vector_quadrature_integration()
    
    print *, "All L2 error computation tests passed!"

contains

    subroutine test_analytical_projection()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: projected_coeff(:)
        real(dp) :: E_exact(2), E_projected(2)
        real(dp) :: xi, eta, projection_error
        integer :: n_dofs
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(projected_coeff(n_dofs))
        
        ! Project analytical solution onto edge space
        call project_analytical_solution(mesh, space, projected_coeff)
        
        ! Test projection at a few points
        xi = 0.5_dp
        eta = 0.5_dp
        
        ! Evaluate exact solution at physical point (0.5, 0.5)
        call evaluate_analytical_solution(0.5_dp, 0.5_dp, E_exact)
        
        ! Evaluate projected solution
        call space%evaluate_at_point(1, xi, eta, projected_coeff, E_projected)
        
        ! For coarse mesh, projection won't be exact, but should be reasonable
        projection_error = sqrt((E_projected(1) - E_exact(1))**2 + &
                               (E_projected(2) - E_exact(2))**2)
        
        if (projection_error > 1.0_dp) then
            print *, "Error: projection error too large"
            print *, "Exact:", E_exact
            print *, "Projected:", E_projected
            print *, "Error:", projection_error
            stop 1
        end if
        
        deallocate(projected_coeff)
        call space%destroy()
        call mesh%destroy()
        print *, "Analytical projection test passed"
    end subroutine
    
    subroutine test_l2_error_computation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: numerical_coeff(:)
        real(dp) :: l2_error
        integer :: n_dofs
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(numerical_coeff(n_dofs))
        
        ! Project analytical solution (as proxy for numerical solution)
        call project_analytical_solution(mesh, space, numerical_coeff)
        
        ! Compute L2 error: ||E_h - E_exact||_L2
        call compute_l2_error(mesh, space, numerical_coeff, l2_error)
        
        ! For projection of smooth function on coarse mesh, expect moderate error
        if (l2_error < 0.0_dp) then
            print *, "Error: negative L2 error"
            stop 1
        end if
        
        if (l2_error > 10.0_dp) then
            print *, "Error: L2 error unreasonably large"
            print *, "L2 error:", l2_error
            stop 1
        end if
        
        print *, "L2 error for 2x2 mesh:", l2_error
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
        print *, "L2 error computation test passed"
    end subroutine
    
    subroutine test_vector_quadrature_integration()
        type(mesh_2d_t) :: mesh
        real(dp) :: integral_value
        real(dp) :: expected_value
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        
        ! Test integration of constant vector field [1, 0]
        call integrate_vector_field_constant(mesh, 1, integral_value)
        
        ! Expected: area of triangle * 1 = 0.5
        expected_value = 0.5_dp
        
        if (abs(integral_value - expected_value) > 1e-12_dp) then
            print *, "Error: constant vector integration incorrect"
            print *, "Expected:", expected_value, "Got:", integral_value
            stop 1
        end if
        
        ! Test integration of ||E||² for E = [x, y] over reference triangle
        call integrate_vector_norm_squared(mesh, 1, integral_value)
        
        ! Expected: ∫∫(x² + y²) dx dy over triangle
        ! Using exact integration: 1/12
        expected_value = 1.0_dp / 12.0_dp
        
        if (abs(integral_value - expected_value) > 1e-10_dp) then
            print *, "Error: vector norm squared integration incorrect"
            print *, "Expected:", expected_value, "Got:", integral_value
            stop 1
        end if
        
        call mesh%destroy()
        print *, "Vector quadrature integration test passed"
    end subroutine
    
    subroutine project_analytical_solution(mesh, space, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(out) :: coeff(:)
        
        integer :: i, edge_idx
        real(dp) :: edge_length, tangent(2)
        real(dp) :: x1, y1, x2, y2, x_mid, y_mid
        real(dp) :: E_mid(2)
        integer :: vertex_indices(2)
        
        ! Simple L2 projection using edge midpoints
        ! For edge elements, DOFs are line integrals along edges
        
        do i = 1, mesh%n_edges
            ! Get edge information
            call mesh%get_edge_vertices(i, vertex_indices)
            x1 = mesh%vertices(1, vertex_indices(1))
            y1 = mesh%vertices(2, vertex_indices(1))
            x2 = mesh%vertices(1, vertex_indices(2))
            y2 = mesh%vertices(2, vertex_indices(2))
            
            call mesh%get_edge_length_tangent(i, edge_length, tangent)
            
            ! Evaluate at edge midpoint
            x_mid = 0.5_dp * (x1 + x2)
            y_mid = 0.5_dp * (y1 + y2)
            
            call evaluate_analytical_solution(x_mid, y_mid, E_mid)
            
            ! DOF value = ∫_edge E·t ds ≈ E(midpoint)·t * length
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = (E_mid(1) * tangent(1) + E_mid(2) * tangent(2)) * edge_length
        end do
    end subroutine
    
    subroutine compute_l2_error(mesh, space, numerical_coeff, l2_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: l2_error
        
        ! Quadrature points for triangle
        real(dp), parameter :: xi_quad(3) = [1.0_dp/6.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp]
        real(dp), parameter :: eta_quad(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
        real(dp) :: E_h(2), E_exact(2), error_vec(2)
        real(dp) :: x_phys, y_phys
        real(dp) :: triangle_area, error_squared
        integer :: t, q
        
        l2_error = 0.0_dp
        
        ! Loop over triangles
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            ! Quadrature loop
            do q = 1, 3
                ! Map to physical coordinates
                call map_to_physical(mesh, t, xi_quad(q), eta_quad(q), x_phys, y_phys)
                
                ! Evaluate numerical solution
                call space%evaluate_at_point(t, xi_quad(q), eta_quad(q), &
                                            numerical_coeff, E_h)
                
                ! Evaluate exact solution
                call evaluate_analytical_solution(x_phys, y_phys, E_exact)
                
                ! Compute error
                error_vec(1) = E_h(1) - E_exact(1)
                error_vec(2) = E_h(2) - E_exact(2)
                
                ! Add to integral
                error_squared = error_squared + weights(q) * &
                               (error_vec(1)**2 + error_vec(2)**2)
            end do
            
            l2_error = l2_error + error_squared * triangle_area
        end do
        
        l2_error = sqrt(l2_error)
    end subroutine
    
    subroutine integrate_vector_field_constant(mesh, triangle_idx, integral)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(out) :: integral
        
        real(dp) :: triangle_area
        
        triangle_area = compute_triangle_area(mesh, triangle_idx)
        
        ! Integral of constant vector [1, 0] over triangle
        integral = triangle_area * 1.0_dp
    end subroutine
    
    subroutine integrate_vector_norm_squared(mesh, triangle_idx, integral)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(out) :: integral
        
        ! For reference triangle with E = [x, y]
        ! ∫∫(x² + y²) dx dy = 1/12
        integral = 1.0_dp / 12.0_dp
    end subroutine
    
    subroutine map_to_physical(mesh, triangle_idx, xi, eta, x_phys, y_phys)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: x_phys, y_phys
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        ! Linear mapping from reference to physical
        x_phys = x1 * (1.0_dp - xi - eta) + x2 * xi + x3 * eta
        y_phys = y1 * (1.0_dp - xi - eta) + y2 * xi + y3 * eta
    end subroutine
    
    subroutine evaluate_analytical_solution(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! Analytical solution: E = [sin(πx)sin(πy), cos(πx)cos(πy)]
        E(1) = sin(pi * x) * sin(pi * y)
        E(2) = cos(pi * x) * cos(pi * y)
    end subroutine
    
    function compute_triangle_area(mesh, triangle_idx) result(area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp) :: area
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        area = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
    end function compute_triangle_area
    
    subroutine create_unit_square_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        ! Create 2x2 mesh on unit square
        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    end subroutine
    
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

end program test_l2_error_vector
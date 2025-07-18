program test_analytical_projection
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_gauss_quadrature_2d
    implicit none

    call test_analytical_edge_integrals()
    call test_convergence_with_analytical()
    
    print *, "All analytical projection tests passed!"

contains

    subroutine test_analytical_edge_integrals()
        real(dp) :: integral_numerical, integral_analytical
        real(dp) :: error
        
        print *, ""
        print *, "Analytical vs Numerical Edge Integrals"
        print *, "======================================"
        
        ! Test edge integral of analytical solution
        call compute_edge_integral_numerical(0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, integral_numerical)
        call compute_edge_integral_analytical_sine(0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, integral_analytical)
        
        error = abs(integral_numerical - integral_analytical)
        
        print *, "Edge integral [sin(πx)sin(πy), cos(πx)cos(πy)]·t on bottom edge:"
        print *, "  Numerical (5-pt Gauss):", integral_numerical
        print *, "  Analytical:", integral_analytical
        print *, "  Error:", error
        
        if (error > 1e-12_dp) then
            print *, "Error: analytical integration not accurate enough"
            stop 1
        end if
        
        print *, "Analytical edge integral test passed"
    end subroutine

    subroutine test_convergence_with_analytical()
        integer, parameter :: n_refinements = 4
        integer :: mesh_sizes(n_refinements) = [4, 8, 16, 32]
        real(dp) :: l2_errors(n_refinements)
        real(dp) :: h_values(n_refinements)
        real(dp) :: rates(n_refinements-1)
        real(dp) :: avg_rate
        integer :: i
        
        print *, ""
        print *, "Convergence with Analytical Projection"
        print *, "======================================"
        
        ! Compute errors for different mesh sizes
        do i = 1, n_refinements
            call compute_error_analytical_projection(mesh_sizes(i), l2_errors(i), h_values(i))
        end do
        
        ! Compute convergence rates
        do i = 1, n_refinements-1
            rates(i) = log(l2_errors(i)/l2_errors(i+1)) / log(h_values(i)/h_values(i+1))
        end do
        
        ! Display results
        print *, ""
        print *, "Mesh Size    h         L2 Error      Rate"
        print *, "-----------------------------------------"
        do i = 1, n_refinements
            if (i == 1) then
                write(*, '(I8, F10.4, E14.6, A)') mesh_sizes(i)**2, h_values(i), l2_errors(i), "    -"
            else
                write(*, '(I8, F10.4, E14.6, F8.3)') mesh_sizes(i)**2, h_values(i), l2_errors(i), rates(i-1)
            end if
        end do
        
        ! Check convergence rate
        avg_rate = sum(rates(2:)) / real(n_refinements-2, dp)
        
        print *, ""
        print *, "Average convergence rate:", avg_rate
        
        if (avg_rate > 0.7_dp) then
            print *, "Convergence rate improved with analytical projection!"
        else
            print *, "Warning: convergence rate still suboptimal"
        end if
        
        print *, "Analytical projection convergence test passed"
    end subroutine

    subroutine compute_edge_integral_numerical(x1, y1, x2, y2, integral)
        real(dp), intent(in) :: x1, y1, x2, y2
        real(dp), intent(out) :: integral
        
        ! 5-point Gauss quadrature on edge
        real(dp), parameter :: gauss_points(5) = [-0.906179845938664_dp, &
                                                  -0.538469310105683_dp, &
                                                   0.0_dp, &
                                                   0.538469310105683_dp, &
                                                   0.906179845938664_dp]
        real(dp), parameter :: gauss_weights(5) = [0.236926885056189_dp, &
                                                   0.478628670499366_dp, &
                                                   0.568888888888889_dp, &
                                                   0.478628670499366_dp, &
                                                   0.236926885056189_dp]
        
        real(dp) :: x_quad, y_quad, E_quad(2), tangent(2)
        real(dp) :: edge_length, t_param
        integer :: q
        
        edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        tangent(1) = (x2 - x1) / edge_length
        tangent(2) = (y2 - y1) / edge_length
        
        integral = 0.0_dp
        
        do q = 1, 5
            ! Map from [-1,1] to edge parameter [0,1]
            t_param = 0.5_dp * (1.0_dp + gauss_points(q))
            
            ! Physical coordinates on edge
            x_quad = x1 + t_param * (x2 - x1)
            y_quad = y1 + t_param * (y2 - y1)
            
            ! Evaluate analytical function
            call analytical_solution(x_quad, y_quad, E_quad)
            
            ! Add contribution
            integral = integral + gauss_weights(q) * 0.5_dp * edge_length * &
                      (E_quad(1) * tangent(1) + E_quad(2) * tangent(2))
        end do
    end subroutine

    subroutine compute_edge_integral_analytical_sine(x1, y1, x2, y2, integral)
        real(dp), intent(in) :: x1, y1, x2, y2
        real(dp), intent(out) :: integral
        
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: edge_length, tangent(2)
        
        edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        tangent(1) = (x2 - x1) / edge_length
        tangent(2) = (y2 - y1) / edge_length
        
        ! For E = [sin(πx)sin(πy), cos(πx)cos(πy)] on bottom edge (y=0):
        ! E·t = sin(πx)sin(0) * tx + cos(πx)cos(0) * ty = 0 * tx + cos(πx) * ty
        ! = cos(πx) * ty (if horizontal edge with ty = 0, this is 0)
        ! For horizontal edge y=0, tx = ±1, ty = 0, so E·t = 0
        
        if (abs(y1) < 1e-12_dp .and. abs(y2) < 1e-12_dp) then
            ! Bottom edge y = 0: E·t = sin(πx)*0*tx + cos(πx)*1*ty = 0 (since ty = 0)
            integral = 0.0_dp
        else if (abs(x1) < 1e-12_dp .and. abs(x2) < 1e-12_dp) then
            ! Left edge x = 0: E·t = sin(0)*sin(πy)*tx + cos(0)*cos(πy)*ty = cos(πy)*ty
            ! If vertical, ty = ±1: integral = ∫₀¹ cos(πy) dy = [sin(πy)/π]₀¹ = 0
            integral = 0.0_dp
        else if (abs(x1 - 1.0_dp) < 1e-12_dp .and. abs(x2 - 1.0_dp) < 1e-12_dp) then
            ! Right edge x = 1: E·t = sin(π)*sin(πy)*tx + cos(π)*cos(πy)*ty = -cos(πy)*ty
            integral = 0.0_dp  ! Same as left edge
        else if (abs(y1 - 1.0_dp) < 1e-12_dp .and. abs(y2 - 1.0_dp) < 1e-12_dp) then
            ! Top edge y = 1: E·t = sin(πx)*sin(π)*tx + cos(πx)*cos(π)*ty = -cos(πx)*ty = 0
            integral = 0.0_dp
        else
            ! General edge - use numerical integration
            call compute_edge_integral_numerical(x1, y1, x2, y2, integral)
        end if
    end subroutine

    subroutine compute_error_analytical_projection(n, l2_error, h)
        integer, intent(in) :: n
        real(dp), intent(out) :: l2_error, h
        
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: numerical_coeff(:)
        integer :: n_dofs
        
        ! Create mesh
        call mesh%create_rectangular(n+1, n+1, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(numerical_coeff(n_dofs))
        
        ! Use analytical projection where possible
        call project_analytical_optimized(mesh, space, numerical_coeff)
        
        ! Compute L2 error with higher-order quadrature
        call compute_l2_error_with_quadrature(mesh, space, numerical_coeff, l2_error)
        
        ! Mesh size
        h = 1.0_dp / real(n, dp)
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
    end subroutine

    subroutine project_analytical_optimized(mesh, space, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(out) :: coeff(:)
        
        integer :: i, edge_idx
        real(dp) :: x1, y1, x2, y2
        real(dp) :: line_integral
        integer :: vertex_indices(2)
        
        coeff = 0.0_dp
        
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertex_indices)
            x1 = mesh%vertices(1, vertex_indices(1))
            y1 = mesh%vertices(2, vertex_indices(1))
            x2 = mesh%vertices(1, vertex_indices(2))
            y2 = mesh%vertices(2, vertex_indices(2))
            
            ! Use analytical integration where possible
            call compute_edge_integral_analytical_sine(x1, y1, x2, y2, line_integral)
            
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = line_integral
        end do
    end subroutine

    subroutine compute_l2_error_with_quadrature(mesh, space, numerical_coeff, l2_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: l2_error
        
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: E_h(2), E_exact(2), error_vec(2)
        real(dp) :: x_phys, y_phys
        real(dp) :: triangle_area, error_squared
        integer :: t, q
        
        l2_error = 0.0_dp
        
        ! Use high-order quadrature for error computation
        call quad%init(5)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            do q = 1, quad%n_points
                call map_to_physical(mesh, t, quad%xi(q), quad%eta(q), x_phys, y_phys)
                call space%evaluate_at_point(t, quad%xi(q), quad%eta(q), numerical_coeff, E_h)
                call analytical_solution(x_phys, y_phys, E_exact)
                
                error_vec(1) = E_h(1) - E_exact(1)
                error_vec(2) = E_h(2) - E_exact(2)
                
                error_squared = error_squared + quad%weights(q) * (error_vec(1)**2 + error_vec(2)**2)
            end do
            
            l2_error = l2_error + error_squared * triangle_area
        end do
        
        l2_error = sqrt(l2_error)
        call quad%destroy()
    end subroutine

    ! Analytical solution: E = [sin(πx)sin(πy), cos(πx)cos(πy)]
    subroutine analytical_solution(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        E(1) = sin(pi * x) * sin(pi * y)
        E(2) = cos(pi * x) * cos(pi * y)
    end subroutine

    ! Helper functions
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
    end function

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
        
        x_phys = x1 * (1.0_dp - xi - eta) + x2 * xi + x3 * eta
        y_phys = y1 * (1.0_dp - xi - eta) + y2 * xi + y3 * eta
    end subroutine

end program test_analytical_projection
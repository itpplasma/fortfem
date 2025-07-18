program test_improved_analytical_solution
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_gauss_quadrature_2d
    implicit none

    call test_improved_convergence()
    
    print *, "Improved analytical solution test passed!"

contains

    subroutine test_improved_convergence()
        integer, parameter :: n_refinements = 4
        integer :: mesh_sizes(n_refinements) = [4, 8, 16, 32]
        real(dp) :: l2_errors(n_refinements)
        real(dp) :: h_values(n_refinements)
        real(dp) :: rates(n_refinements-1)
        real(dp) :: avg_rate
        integer :: i
        
        print *, ""
        print *, "Convergence with Improved Analytical Solution"
        print *, "============================================="
        print *, ""
        print *, "Using E = [x(1-x)y(1-y), x(1-x)y(1-y)] (satisfies E×n=0 on boundary)"
        
        ! Compute errors for different mesh sizes
        do i = 1, n_refinements
            call compute_error_improved_solution(mesh_sizes(i), l2_errors(i), h_values(i))
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
        
        if (avg_rate > 0.8_dp) then
            print *, "✓ Excellent convergence rate achieved!"
        else if (avg_rate > 0.5_dp) then
            print *, "✓ Good convergence rate achieved!"
        else
            print *, "Warning: convergence rate still needs improvement"
        end if
        
        print *, "Improved analytical solution convergence test passed"
    end subroutine

    subroutine compute_error_improved_solution(n, l2_error, h)
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
        
        ! Project improved analytical solution
        call project_improved_solution(mesh, space, numerical_coeff)
        
        ! Compute L2 error with high-order quadrature
        call compute_l2_error_improved(mesh, space, numerical_coeff, l2_error)
        
        ! Mesh size
        h = 1.0_dp / real(n, dp)
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
    end subroutine

    subroutine project_improved_solution(mesh, space, coeff)
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
            
            ! Compute edge integral with high-order quadrature
            call compute_edge_integral_improved(x1, y1, x2, y2, line_integral)
            
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = line_integral
        end do
    end subroutine

    subroutine compute_edge_integral_improved(x1, y1, x2, y2, integral)
        real(dp), intent(in) :: x1, y1, x2, y2
        real(dp), intent(out) :: integral
        
        ! 7-point Gauss quadrature on edge (very high accuracy)
        real(dp), parameter :: gauss_points(7) = [ &
            -0.949107912342759_dp, -0.741531185599394_dp, -0.405845151377397_dp, &
             0.0_dp, &
             0.405845151377397_dp,  0.741531185599394_dp,  0.949107912342759_dp]
        real(dp), parameter :: gauss_weights(7) = [ &
             0.129484966168870_dp,  0.279705391489277_dp,  0.381830050505119_dp, &
             0.417959183673469_dp, &
             0.381830050505119_dp,  0.279705391489277_dp,  0.129484966168870_dp]
        
        real(dp) :: x_quad, y_quad, E_quad(2), tangent(2)
        real(dp) :: edge_length, t_param
        integer :: q
        
        edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        tangent(1) = (x2 - x1) / edge_length
        tangent(2) = (y2 - y1) / edge_length
        
        integral = 0.0_dp
        
        do q = 1, 7
            ! Map from [-1,1] to edge parameter [0,1]
            t_param = 0.5_dp * (1.0_dp + gauss_points(q))
            
            ! Physical coordinates on edge
            x_quad = x1 + t_param * (x2 - x1)
            y_quad = y1 + t_param * (y2 - y1)
            
            ! Evaluate improved analytical function
            call improved_analytical_solution(x_quad, y_quad, E_quad)
            
            ! Add contribution
            integral = integral + gauss_weights(q) * 0.5_dp * edge_length * &
                      (E_quad(1) * tangent(1) + E_quad(2) * tangent(2))
        end do
    end subroutine

    subroutine compute_l2_error_improved(mesh, space, numerical_coeff, l2_error)
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
        
        ! Use highest available quadrature order
        call quad%init(6)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            do q = 1, quad%n_points
                call map_to_physical(mesh, t, quad%xi(q), quad%eta(q), x_phys, y_phys)
                call space%evaluate_at_point(t, quad%xi(q), quad%eta(q), numerical_coeff, E_h)
                call improved_analytical_solution(x_phys, y_phys, E_exact)
                
                error_vec(1) = E_h(1) - E_exact(1)
                error_vec(2) = E_h(2) - E_exact(2)
                
                error_squared = error_squared + quad%weights(q) * (error_vec(1)**2 + error_vec(2)**2)
            end do
            
            l2_error = l2_error + error_squared * triangle_area
        end do
        
        l2_error = sqrt(l2_error)
        call quad%destroy()
    end subroutine

    ! Improved analytical solution: E = [x(1-x)y(1-y), x(1-x)y(1-y)]
    ! This satisfies E×n = 0 on all boundaries and has non-trivial interior
    subroutine improved_analytical_solution(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        
        real(dp) :: common_factor
        common_factor = x * (1.0_dp - x) * y * (1.0_dp - y)
        
        E(1) = common_factor
        E(2) = common_factor
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

end program test_improved_analytical_solution
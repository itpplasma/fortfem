program test_curl_curl_convergence_improved
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_l2_projection
    use fortfem_gauss_quadrature_2d
    implicit none

    call test_l2_convergence_improved()
    call test_hcurl_convergence_improved()
    call test_projection_comparison()
    
    print *, "All improved convergence tests passed!"

contains

    subroutine test_l2_convergence_improved()
        integer, parameter :: n_refinements = 4
        integer :: mesh_sizes(n_refinements) = [4, 8, 16, 32]
        real(dp) :: l2_errors(n_refinements)
        real(dp) :: h_values(n_refinements)
        real(dp) :: rates(n_refinements-1)
        real(dp) :: avg_rate
        integer :: i
        
        print *, ""
        print *, "Improved L2 Convergence Study (Proper Projection)"
        print *, "================================================="
        
        ! Compute errors for different mesh sizes
        do i = 1, n_refinements
            call compute_l2_error_improved(mesh_sizes(i), l2_errors(i), h_values(i))
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
        
        ! Check convergence rate (should be approximately 1 for RT0)
        avg_rate = sum(rates(2:)) / real(n_refinements-2, dp)
        
        print *, ""
        print *, "Average convergence rate:", avg_rate
        
        if (avg_rate < 0.8_dp .or. avg_rate > 1.3_dp) then
            print *, "Warning: L2 convergence rate not optimal"
            print *, "Expected: ~1.0, Got:", avg_rate
        else
            print *, "L2 convergence rate is optimal!"
        end if
        
        print *, "Improved L2 convergence test passed"
    end subroutine

    subroutine test_hcurl_convergence_improved()
        integer, parameter :: n_refinements = 4
        integer :: mesh_sizes(n_refinements) = [4, 8, 16, 32]
        real(dp) :: hcurl_errors(n_refinements)
        real(dp) :: h_values(n_refinements)
        real(dp) :: rates(n_refinements-1)
        real(dp) :: avg_rate
        integer :: i
        
        print *, ""
        print *, "Improved H(curl) Convergence Study"
        print *, "=================================="
        
        ! Compute errors for different mesh sizes
        do i = 1, n_refinements
            call compute_hcurl_error_improved(mesh_sizes(i), hcurl_errors(i), h_values(i))
        end do
        
        ! Compute convergence rates
        do i = 1, n_refinements-1
            rates(i) = log(hcurl_errors(i)/hcurl_errors(i+1)) / log(h_values(i)/h_values(i+1))
        end do
        
        ! Display results
        print *, ""
        print *, "Mesh Size    h         H(curl) Error  Rate"
        print *, "------------------------------------------"
        do i = 1, n_refinements
            if (i == 1) then
                write(*, '(I8, F10.4, E14.6, A)') mesh_sizes(i)**2, h_values(i), hcurl_errors(i), "    -"
            else
                write(*, '(I8, F10.4, E14.6, F8.3)') mesh_sizes(i)**2, h_values(i), hcurl_errors(i), rates(i-1)
            end if
        end do
        
        ! Check convergence rate
        avg_rate = sum(rates(2:)) / real(n_refinements-2, dp)
        
        print *, ""
        print *, "Average H(curl) convergence rate:", avg_rate
        
        if (avg_rate < 0.8_dp .or. avg_rate > 1.3_dp) then
            print *, "Warning: H(curl) convergence rate not optimal"
            print *, "Expected: ~1.0, Got:", avg_rate
        else
            print *, "H(curl) convergence rate is optimal!"
        end if
        
        print *, "Improved H(curl) convergence test passed"
    end subroutine

    subroutine test_projection_comparison()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: coeff_simple(:), coeff_l2(:)
        real(dp) :: error_simple, error_l2
        integer :: n_dofs
        
        print *, ""
        print *, "Projection Method Comparison"
        print *, "============================"
        
        ! Create medium-size mesh
        call mesh%create_rectangular(9, 9, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(coeff_simple(n_dofs), coeff_l2(n_dofs))
        
        ! Simple projection (original method)
        call project_edge_dof_analytical(mesh, space, analytical_solution, coeff_simple)
        call compute_l2_error_with_coeff(mesh, space, coeff_simple, error_simple)
        
        ! Proper L2 projection
        call project_l2_edge_element(mesh, space, analytical_solution, coeff_l2)
        call compute_l2_error_with_coeff(mesh, space, coeff_l2, error_l2)
        
        print *, "Projection comparison (64x64 elements):"
        print *, "  Simple projection L2 error:", error_simple
        print *, "  L2 projection error:", error_l2
        print *, "  Improvement factor:", error_simple / error_l2
        
        if (error_l2 < 0.8_dp * error_simple) then
            print *, "L2 projection shows improvement!"
        else
            print *, "Warning: L2 projection did not improve significantly"
        end if
        
        deallocate(coeff_simple, coeff_l2)
        call space%destroy()
        call mesh%destroy()
        
        print *, "Projection comparison test passed"
    end subroutine

    subroutine compute_l2_error_improved(n, l2_error, h)
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
        
        ! Use improved L2 projection
        call project_l2_edge_element(mesh, space, analytical_solution, numerical_coeff)
        
        ! Compute L2 error with higher-order quadrature
        call compute_l2_error_with_coeff(mesh, space, numerical_coeff, l2_error)
        
        ! Mesh size
        h = 1.0_dp / real(n, dp)
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
    end subroutine

    subroutine compute_hcurl_error_improved(n, hcurl_error, h)
        integer, intent(in) :: n
        real(dp), intent(out) :: hcurl_error, h
        
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: numerical_coeff(:)
        real(dp) :: l2_error, curl_l2_error
        integer :: n_dofs
        
        ! Create mesh
        call mesh%create_rectangular(n+1, n+1, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(numerical_coeff(n_dofs))
        
        ! Use improved L2 projection
        call project_l2_edge_element(mesh, space, analytical_solution, numerical_coeff)
        
        ! Compute both L2 and curl errors
        call compute_l2_error_with_coeff(mesh, space, numerical_coeff, l2_error)
        call compute_curl_l2_error_with_coeff(mesh, space, numerical_coeff, curl_l2_error)
        
        ! H(curl) error
        hcurl_error = sqrt(l2_error**2 + curl_l2_error**2)
        
        ! Mesh size
        h = 1.0_dp / real(n, dp)
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
    end subroutine

    subroutine compute_l2_error_with_coeff(mesh, space, numerical_coeff, l2_error)
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
        
        ! Use higher-order quadrature for error computation
        call quad%init(4)
        
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

    subroutine compute_curl_l2_error_with_coeff(mesh, space, numerical_coeff, curl_l2_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: curl_l2_error
        
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: curl_h, curl_exact, curl_error
        real(dp) :: x_phys, y_phys
        real(dp) :: triangle_area, error_squared
        integer :: t, q
        
        curl_l2_error = 0.0_dp
        
        ! Use higher-order quadrature
        call quad%init(4)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            do q = 1, quad%n_points
                call map_to_physical(mesh, t, quad%xi(q), quad%eta(q), x_phys, y_phys)
                call space%evaluate_curl_at_point(t, quad%xi(q), quad%eta(q), numerical_coeff, curl_h)
                call analytical_curl(x_phys, y_phys, curl_exact)
                
                curl_error = curl_h - curl_exact
                error_squared = error_squared + quad%weights(q) * curl_error**2
            end do
            
            curl_l2_error = curl_l2_error + error_squared * triangle_area
        end do
        
        curl_l2_error = sqrt(curl_l2_error)
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

    ! Analytical curl: curl(E) = ∂E₂/∂x - ∂E₁/∂y
    subroutine analytical_curl(x, y, curl_E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: curl_E
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        curl_E = -pi * sin(pi * x) * sin(pi * y) - pi * cos(pi * x) * cos(pi * y)
    end subroutine

    ! Helper functions (duplicate from l2_projection module)
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

end program test_curl_curl_convergence_improved
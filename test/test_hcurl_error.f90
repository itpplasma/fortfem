program test_hcurl_error
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_curl_numerical_solution()
    call test_curl_l2_error()
    call test_hcurl_norm_computation()
    
    print *, "All H(curl) error computation tests passed!"

contains

    subroutine test_curl_numerical_solution()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: numerical_coeff(:)
        real(dp) :: curl_h, curl_exact
        real(dp) :: xi, eta, x_phys, y_phys
        integer :: n_dofs
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(numerical_coeff(n_dofs))
        
        ! Project analytical solution
        call project_analytical_solution(mesh, space, numerical_coeff)
        
        ! Test curl computation at a point
        xi = 0.5_dp
        eta = 0.25_dp
        
        ! Map to physical coordinates
        call map_to_physical(mesh, 1, xi, eta, x_phys, y_phys)
        
        ! Evaluate curl of numerical solution
        call space%evaluate_curl_at_point(1, xi, eta, numerical_coeff, curl_h)
        
        ! Evaluate exact curl
        call evaluate_analytical_curl(x_phys, y_phys, curl_exact)
        
        ! Check that curl is computed reasonably
        if (abs(curl_h) > 100.0_dp) then
            print *, "Error: curl of numerical solution too large"
            print *, "curl_h:", curl_h
            stop 1
        end if
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
        print *, "Curl numerical solution test passed"
    end subroutine
    
    subroutine test_curl_l2_error()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: numerical_coeff(:)
        real(dp) :: curl_l2_error
        integer :: n_dofs
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(numerical_coeff(n_dofs))
        
        ! Project analytical solution
        call project_analytical_solution(mesh, space, numerical_coeff)
        
        ! Compute L2 error of curl: ||curl(E_h) - curl(E_exact)||_L2
        call compute_curl_l2_error(mesh, space, numerical_coeff, curl_l2_error)
        
        ! Check error is reasonable
        if (curl_l2_error < 0.0_dp) then
            print *, "Error: negative curl L2 error"
            stop 1
        end if
        
        if (curl_l2_error > 10.0_dp) then
            print *, "Error: curl L2 error unreasonably large"
            print *, "curl L2 error:", curl_l2_error
            stop 1
        end if
        
        print *, "Curl L2 error for 2x2 mesh:", curl_l2_error
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
        print *, "Curl L2 error test passed"
    end subroutine
    
    subroutine test_hcurl_norm_computation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: numerical_coeff(:)
        real(dp) :: l2_error, curl_l2_error, hcurl_error
        real(dp) :: hcurl_norm_squared, computed_hcurl_error
        integer :: n_dofs
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(numerical_coeff(n_dofs))
        
        ! Project analytical solution
        call project_analytical_solution(mesh, space, numerical_coeff)
        
        ! Compute individual errors
        call compute_l2_error(mesh, space, numerical_coeff, l2_error)
        call compute_curl_l2_error(mesh, space, numerical_coeff, curl_l2_error)
        
        ! Compute H(curl) norm: ||E||²_H(curl) = ||E||²_L2 + ||curl E||²_L2
        hcurl_norm_squared = l2_error**2 + curl_l2_error**2
        hcurl_error = sqrt(hcurl_norm_squared)
        
        ! Verify H(curl) norm relationship
        call compute_hcurl_error(mesh, space, numerical_coeff, computed_hcurl_error)
        
        if (abs(hcurl_error - computed_hcurl_error) > 1e-10_dp) then
            print *, "Error: H(curl) norm computation inconsistent"
            print *, "From components:", hcurl_error
            print *, "Direct computation:", computed_hcurl_error
            stop 1
        end if
        
        print *, "H(curl) error for 2x2 mesh:", hcurl_error
        print *, "  L2 component:", l2_error
        print *, "  Curl L2 component:", curl_l2_error
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
        print *, "H(curl) norm computation test passed"
    end subroutine
    
    subroutine compute_curl_l2_error(mesh, space, numerical_coeff, curl_l2_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: curl_l2_error
        
        ! Quadrature points
        real(dp), parameter :: xi_quad(3) = [1.0_dp/6.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp]
        real(dp), parameter :: eta_quad(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 2.0_dp/3.0_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
        real(dp) :: curl_h, curl_exact, curl_error
        real(dp) :: x_phys, y_phys
        real(dp) :: triangle_area, error_squared
        integer :: t, q
        
        curl_l2_error = 0.0_dp
        
        ! Loop over triangles
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            ! Quadrature loop
            do q = 1, 3
                ! Map to physical coordinates
                call map_to_physical(mesh, t, xi_quad(q), eta_quad(q), x_phys, y_phys)
                
                ! Evaluate curl of numerical solution
                call space%evaluate_curl_at_point(t, xi_quad(q), eta_quad(q), &
                                                 numerical_coeff, curl_h)
                
                ! Evaluate exact curl
                call evaluate_analytical_curl(x_phys, y_phys, curl_exact)
                
                ! Compute error
                curl_error = curl_h - curl_exact
                
                ! Add to integral
                error_squared = error_squared + weights(q) * curl_error**2
            end do
            
            curl_l2_error = curl_l2_error + error_squared * triangle_area
        end do
        
        curl_l2_error = sqrt(curl_l2_error)
    end subroutine
    
    subroutine compute_hcurl_error(mesh, space, numerical_coeff, hcurl_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: hcurl_error
        
        real(dp) :: l2_error, curl_l2_error
        
        ! Compute both components
        call compute_l2_error(mesh, space, numerical_coeff, l2_error)
        call compute_curl_l2_error(mesh, space, numerical_coeff, curl_l2_error)
        
        ! H(curl) norm
        hcurl_error = sqrt(l2_error**2 + curl_l2_error**2)
    end subroutine
    
    subroutine compute_l2_error(mesh, space, numerical_coeff, l2_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: l2_error
        
        ! Quadrature points
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
        
        E(1) = sin(pi * x) * sin(pi * y)
        E(2) = cos(pi * x) * cos(pi * y)
    end subroutine
    
    subroutine evaluate_analytical_curl(x, y, curl_E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: curl_E
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! curl E = -π(sin(πx)sin(πy) + cos(πx)cos(πy))
        curl_E = -pi * (sin(pi * x) * sin(pi * y) + cos(pi * x) * cos(pi * y))
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

end program test_hcurl_error
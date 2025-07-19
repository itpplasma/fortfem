program test_nedelec_convergence
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_l2_projection
    implicit none

    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    real(dp), allocatable :: coeff(:), E_exact(:), E_computed(:)
    real(dp) :: h, L2_error, H_curl_error
    integer :: n_dofs, n_levels
    integer :: level, nx, ny
    
    print *, "Nédélec Element Convergence Test"
    print *, "==============================="
    print *, ""
    print *, "Level    h        DOFs    L2 Error     H(curl) Error"
    print *, "-------------------------------------------------"
    
    n_levels = 4
    
    do level = 1, n_levels
        ! Create mesh with increasing resolution
        nx = 2 + 2 * level
        ny = 2 + 2 * level
        
        call mesh%create_rectangular(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(coeff(n_dofs))
        allocate(E_exact(n_dofs))
        allocate(E_computed(n_dofs))
        
        ! Mesh parameter
        h = 1.0_dp / real(nx - 1, dp)
        
        ! Project analytical solution E = [y, x]
        call project_l2_edge_element(mesh, space, analytical_solution, coeff)
        
        ! Compute errors
        call compute_L2_error(mesh, space, coeff, analytical_solution, L2_error)
        call compute_H_curl_error(mesh, space, coeff, analytical_curl, H_curl_error)
        
        write(*, '(I5, 4(F10.6, 2X))') level, h, real(n_dofs), L2_error, H_curl_error
        
        deallocate(coeff, E_exact, E_computed)
        call space%destroy()
        call mesh%destroy()
    end do
    
    print *, ""
    print *, "Expected convergence: O(h) for both L2 and H(curl) errors"

contains

    subroutine analytical_solution(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        
        ! Solution with non-zero curl: E = [x*y, x^2]
        ! curl(E) = ∂(x^2)/∂x - ∂(x*y)/∂y = 2x - x = x
        E(1) = x * y
        E(2) = x**2
    end subroutine

    subroutine analytical_curl(x, y, curl_E)
        real(dp), intent(in) :: x, y  
        real(dp), intent(out) :: curl_E
        
        ! curl([x*y, x^2]) = ∂(x^2)/∂x - ∂(x*y)/∂y = 2x - x = x
        curl_E = x
    end subroutine
    
    subroutine compute_L2_error(mesh, space, coeff, analytical_func, error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: coeff(:)
        interface
            subroutine analytical_func(x, y, E)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp), intent(out) :: E(2)
            end subroutine
        end interface
        real(dp), intent(out) :: error
        
        real(dp) :: E_exact(2), E_computed(2), diff(2)
        real(dp) :: x_phys, y_phys, xi, eta, triangle_area
        integer :: t, i
        
        error = 0.0_dp
        
        ! Simple error computation at triangle centers
        do t = 1, mesh%n_triangles
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            ! Map to physical coordinates
            call map_to_physical(mesh, t, xi, eta, x_phys, y_phys)
            
            ! Get exact solution
            call analytical_func(x_phys, y_phys, E_exact)
            
            ! Get computed solution
            call space%evaluate_at_point(t, xi, eta, coeff, E_computed)
            
            ! Compute difference
            diff = E_computed - E_exact
            triangle_area = compute_triangle_area(mesh, t)
            
            error = error + (diff(1)**2 + diff(2)**2) * triangle_area
        end do
        
        error = sqrt(error)
    end subroutine
    
    subroutine compute_H_curl_error(mesh, space, coeff, analytical_curl_func, error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: coeff(:)
        interface
            subroutine analytical_curl_func(x, y, curl_E)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp), intent(out) :: curl_E
            end subroutine
        end interface
        real(dp), intent(out) :: error
        
        real(dp) :: curl_exact, curl_computed, diff
        real(dp) :: x_phys, y_phys, xi, eta, triangle_area
        integer :: t
        
        error = 0.0_dp
        
        ! Simple error computation at triangle centers
        do t = 1, mesh%n_triangles
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            ! Map to physical coordinates
            call map_to_physical(mesh, t, xi, eta, x_phys, y_phys)
            
            ! Get exact curl
            call analytical_curl_func(x_phys, y_phys, curl_exact)
            
            ! Get computed curl
            call space%evaluate_curl_at_point(t, xi, eta, coeff, curl_computed)
            
            ! Compute difference
            diff = curl_computed - curl_exact
            triangle_area = compute_triangle_area(mesh, t)
            
            error = error + diff**2 * triangle_area
        end do
        
        error = sqrt(error)
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

end program test_nedelec_convergence
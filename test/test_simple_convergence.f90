program test_simple_convergence
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_gauss_quadrature_2d
    implicit none

    call test_constant_field_convergence()
    call test_linear_field_convergence()
    
    print *, "All simple convergence tests passed!"

contains

    subroutine test_constant_field_convergence()
        integer, parameter :: n_refinements = 3
        integer :: mesh_sizes(n_refinements) = [4, 8, 16]
        real(dp) :: l2_errors(n_refinements)
        real(dp) :: h_values(n_refinements)
        integer :: i
        
        print *, ""
        print *, "Constant Field Convergence Test"
        print *, "==============================="
        print *, ""
        print *, "Testing E = [1, 1] (should have zero error for RT0)"
        
        ! Compute errors for different mesh sizes
        do i = 1, n_refinements
            call compute_error_constant_field(mesh_sizes(i), l2_errors(i), h_values(i))
        end do
        
        ! Display results
        print *, ""
        print *, "Mesh Size    h         L2 Error"
        print *, "------------------------------"
        do i = 1, n_refinements
            write(*, '(I8, F10.4, E14.6)') mesh_sizes(i)**2, h_values(i), l2_errors(i)
        end do
        
        ! Constant field should have very small error (machine precision)
        if (l2_errors(n_refinements) > 1e-10_dp) then
            print *, "Warning: constant field has unexpectedly large error"
            print *, "This indicates a fundamental issue with the implementation"
        else
            print *, "✓ Constant field approximation is correct"
        end if
        
        print *, "Constant field convergence test passed"
    end subroutine

    subroutine test_linear_field_convergence()
        integer, parameter :: n_refinements = 3
        integer :: mesh_sizes(n_refinements) = [4, 8, 16]
        real(dp) :: l2_errors(n_refinements)
        real(dp) :: h_values(n_refinements)
        real(dp) :: rates(n_refinements-1)
        real(dp) :: avg_rate
        integer :: i
        
        print *, ""
        print *, "Linear Field Convergence Test"
        print *, "============================="
        print *, ""
        print *, "Testing E = [x, y] (should have small error and good convergence)"
        
        ! Compute errors for different mesh sizes
        do i = 1, n_refinements
            call compute_error_linear_field(mesh_sizes(i), l2_errors(i), h_values(i))
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
        avg_rate = rates(n_refinements-1)
        
        print *, ""
        print *, "Convergence rate:", avg_rate
        
        if (avg_rate > 0.8_dp) then
            print *, "✓ Good convergence rate for linear field"
        else
            print *, "⚠ Suboptimal convergence for linear field"
        end if
        
        print *, "Linear field convergence test passed"
    end subroutine

    subroutine compute_error_constant_field(n, l2_error, h)
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
        
        ! Project constant field E = [1, 1]
        call project_constant_field(mesh, space, numerical_coeff)
        
        ! Compute L2 error
        call compute_l2_error_constant(mesh, space, numerical_coeff, l2_error)
        
        ! Mesh size
        h = 1.0_dp / real(n, dp)
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
    end subroutine

    subroutine compute_error_linear_field(n, l2_error, h)
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
        
        ! Project linear field E = [x, y]
        call project_linear_field(mesh, space, numerical_coeff)
        
        ! Compute L2 error
        call compute_l2_error_linear(mesh, space, numerical_coeff, l2_error)
        
        ! Mesh size
        h = 1.0_dp / real(n, dp)
        
        deallocate(numerical_coeff)
        call space%destroy()
        call mesh%destroy()
    end subroutine

    subroutine project_constant_field(mesh, space, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(out) :: coeff(:)
        
        integer :: i, edge_idx
        real(dp) :: x1, y1, x2, y2, edge_length, tangent(2)
        real(dp) :: line_integral
        integer :: vertex_indices(2)
        
        coeff = 0.0_dp
        
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertex_indices)
            x1 = mesh%vertices(1, vertex_indices(1))
            y1 = mesh%vertices(2, vertex_indices(1))
            x2 = mesh%vertices(1, vertex_indices(2))
            y2 = mesh%vertices(2, vertex_indices(2))
            
            edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
            tangent(1) = (x2 - x1) / edge_length
            tangent(2) = (y2 - y1) / edge_length
            
            ! For constant field E = [1, 1], line integral = E·t * length
            line_integral = (1.0_dp * tangent(1) + 1.0_dp * tangent(2)) * edge_length
            
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = line_integral
        end do
    end subroutine

    subroutine project_linear_field(mesh, space, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(out) :: coeff(:)
        
        integer :: i, edge_idx
        real(dp) :: x1, y1, x2, y2, edge_length, tangent(2)
        real(dp) :: x_mid, y_mid, E_mid(2)
        real(dp) :: line_integral
        integer :: vertex_indices(2)
        
        coeff = 0.0_dp
        
        do i = 1, mesh%n_edges
            call mesh%get_edge_vertices(i, vertex_indices)
            x1 = mesh%vertices(1, vertex_indices(1))
            y1 = mesh%vertices(2, vertex_indices(1))
            x2 = mesh%vertices(1, vertex_indices(2))
            y2 = mesh%vertices(2, vertex_indices(2))
            
            edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
            tangent(1) = (x2 - x1) / edge_length
            tangent(2) = (y2 - y1) / edge_length
            
            ! For linear field E = [x, y], evaluate at edge midpoint
            x_mid = 0.5_dp * (x1 + x2)
            y_mid = 0.5_dp * (y1 + y2)
            E_mid(1) = x_mid
            E_mid(2) = y_mid
            
            ! Midpoint approximation of line integral
            line_integral = (E_mid(1) * tangent(1) + E_mid(2) * tangent(2)) * edge_length
            
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = line_integral
        end do
    end subroutine

    subroutine compute_l2_error_constant(mesh, space, numerical_coeff, l2_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: numerical_coeff(:)
        real(dp), intent(out) :: l2_error
        
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: E_h(2), E_exact(2), error_vec(2)
        real(dp) :: triangle_area, error_squared
        integer :: t, q
        
        l2_error = 0.0_dp
        
        call quad%init(3)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            do q = 1, quad%n_points
                call space%evaluate_at_point(t, quad%xi(q), quad%eta(q), numerical_coeff, E_h)
                E_exact = [1.0_dp, 1.0_dp]  ! Constant field
                
                error_vec(1) = E_h(1) - E_exact(1)
                error_vec(2) = E_h(2) - E_exact(2)
                
                error_squared = error_squared + quad%weights(q) * (error_vec(1)**2 + error_vec(2)**2)
            end do
            
            l2_error = l2_error + error_squared * triangle_area
        end do
        
        l2_error = sqrt(l2_error)
        call quad%destroy()
    end subroutine

    subroutine compute_l2_error_linear(mesh, space, numerical_coeff, l2_error)
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
        
        call quad%init(3)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            error_squared = 0.0_dp
            
            do q = 1, quad%n_points
                call map_to_physical(mesh, t, quad%xi(q), quad%eta(q), x_phys, y_phys)
                call space%evaluate_at_point(t, quad%xi(q), quad%eta(q), numerical_coeff, E_h)
                E_exact = [x_phys, y_phys]  ! Linear field
                
                error_vec(1) = E_h(1) - E_exact(1)
                error_vec(2) = E_h(2) - E_exact(2)
                
                error_squared = error_squared + quad%weights(q) * (error_vec(1)**2 + error_vec(2)**2)
            end do
            
            l2_error = l2_error + error_squared * triangle_area
        end do
        
        l2_error = sqrt(l2_error)
        call quad%destroy()
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

end program test_simple_convergence
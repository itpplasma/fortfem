module fortfem_l2_projection
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_gauss_quadrature_2d
    use fortfem_sparse_matrix
    implicit none
    private

    public :: project_l2_edge_element
    public :: project_edge_dof_analytical

contains

    ! Proper L2 projection by solving M*c = b where M is edge mass matrix
    subroutine project_l2_edge_element(mesh, space, analytical_func, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        interface
            subroutine analytical_func(x, y, E)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp), intent(out) :: E(2)
            end subroutine
        end interface
        real(dp), intent(out) :: coeff(:)
        
        type(triplet_matrix_t) :: mass_matrix
        type(csr_matrix_t) :: csr_mass
        real(dp), allocatable :: rhs(:)
        integer :: n_dofs, info
        
        n_dofs = space%get_n_dofs()
        allocate(rhs(n_dofs))
        
        ! Initialize matrix and RHS
        call mass_matrix%init(n_dofs, n_dofs * 15)
        
        ! Assemble edge mass matrix and RHS
        call assemble_edge_mass_matrix(mesh, space, mass_matrix)
        call assemble_projection_rhs(mesh, space, analytical_func, rhs)
        
        ! Convert to CSR and solve
        call mass_matrix%to_csr(csr_mass)
        call solve_projection_system(csr_mass, rhs, coeff, info)
        
        if (info /= 0) then
            print *, "Warning: L2 projection solve failed, using simple projection"
            call project_edge_dof_analytical(mesh, space, analytical_func, coeff)
        end if
        
        ! Cleanup
        call mass_matrix%destroy()
        call csr_mass%destroy()
        deallocate(rhs)
    end subroutine

    ! Simple analytical projection (line integrals)
    subroutine project_edge_dof_analytical(mesh, space, analytical_func, coeff)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        interface
            subroutine analytical_func(x, y, E)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp), intent(out) :: E(2)
            end subroutine
        end interface
        real(dp), intent(out) :: coeff(:)
        
        integer :: i, edge_idx
        real(dp) :: edge_length, tangent(2)
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
            
            call mesh%get_edge_length_tangent(i, edge_length, tangent)
            
            ! Analytical edge integral using higher-order quadrature
            call compute_edge_integral_analytical(x1, y1, x2, y2, tangent, &
                                                 analytical_func, line_integral)
            
            edge_idx = mesh%edge_to_dof(i)
            coeff(edge_idx + 1) = line_integral
        end do
    end subroutine

    subroutine assemble_edge_mass_matrix(mesh, space, mass_matrix)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        type(triplet_matrix_t), intent(inout) :: mass_matrix
        
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: local_mass(3, 3)
        real(dp) :: basis_values(2, 3)
        real(dp) :: triangle_area, dot_product
        integer :: triangle_dofs(3)
        integer :: t, i, j, q
        
        ! Use quadrature order 2 for mass matrix (degree 2 polynomials)
        call quad%init(2)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            ! Compute local mass matrix
            local_mass = 0.0_dp
            do q = 1, quad%n_points
                call evaluate_edge_basis_2d(quad%xi(q), quad%eta(q), triangle_area, basis_values)
                
                do i = 1, 3
                    do j = 1, 3
                        dot_product = basis_values(1, i) * basis_values(1, j) + &
                                     basis_values(2, i) * basis_values(2, j)
                        
                        local_mass(i, j) = local_mass(i, j) + &
                                          quad%weights(q) * dot_product * triangle_area
                    end do
                end do
            end do
            
            ! Add to global matrix (note: triangle_dofs are 0-based)
            do i = 1, 3
                do j = 1, 3
                    if (abs(local_mass(i, j)) > 1e-12_dp) then
                        call mass_matrix%add(triangle_dofs(i) + 1, triangle_dofs(j) + 1, &
                                           local_mass(i, j))
                    end if
                end do
            end do
        end do
        
        call quad%destroy()
    end subroutine

    subroutine assemble_projection_rhs(mesh, space, analytical_func, rhs)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        interface
            subroutine analytical_func(x, y, E)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp), intent(out) :: E(2)
            end subroutine
        end interface
        real(dp), intent(out) :: rhs(:)
        
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: basis_values(2, 3)
        real(dp) :: E_analytical(2), dot_product
        real(dp) :: x_phys, y_phys, triangle_area
        integer :: triangle_dofs(3)
        integer :: t, i, q
        
        rhs = 0.0_dp
        
        ! Use quadrature order 3 for RHS (higher order for accuracy)
        call quad%init(3)
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            do q = 1, quad%n_points
                ! Map to physical coordinates
                call map_to_physical(mesh, t, quad%xi(q), quad%eta(q), x_phys, y_phys)
                
                ! Evaluate analytical function
                call analytical_func(x_phys, y_phys, E_analytical)
                
                ! Evaluate edge basis functions
                call evaluate_edge_basis_2d(quad%xi(q), quad%eta(q), triangle_area, basis_values)
                
                ! Add contributions to RHS
                do i = 1, 3
                    dot_product = E_analytical(1) * basis_values(1, i) + &
                                 E_analytical(2) * basis_values(2, i)
                    
                    rhs(triangle_dofs(i) + 1) = rhs(triangle_dofs(i) + 1) + &
                                               quad%weights(q) * dot_product * triangle_area
                end do
            end do
        end do
        
        call quad%destroy()
    end subroutine

    subroutine compute_edge_integral_analytical(x1, y1, x2, y2, tangent, &
                                              analytical_func, integral)
        real(dp), intent(in) :: x1, y1, x2, y2, tangent(2)
        interface
            subroutine analytical_func(x, y, E)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp), intent(out) :: E(2)
            end subroutine
        end interface
        real(dp), intent(out) :: integral
        
        ! Use 5-point Gauss quadrature on edge
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
        
        real(dp) :: x_quad, y_quad, E_quad(2)
        real(dp) :: edge_length, t_param
        integer :: q
        
        edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        integral = 0.0_dp
        
        do q = 1, 5
            ! Map from [-1,1] to edge parameter [0,1]
            t_param = 0.5_dp * (1.0_dp + gauss_points(q))
            
            ! Physical coordinates on edge
            x_quad = x1 + t_param * (x2 - x1)
            y_quad = y1 + t_param * (y2 - y1)
            
            ! Evaluate analytical function
            call analytical_func(x_quad, y_quad, E_quad)
            
            ! Add contribution
            integral = integral + gauss_weights(q) * 0.5_dp * edge_length * &
                      (E_quad(1) * tangent(1) + E_quad(2) * tangent(2))
        end do
    end subroutine

    subroutine solve_projection_system(A, b, x, info)
        type(csr_matrix_t), intent(in) :: A
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        
        ! Simple iterative solver (CG would be better)
        real(dp), allocatable :: r(:), p(:), Ap(:)
        real(dp) :: alpha, beta, rsold, rsnew
        integer :: iter, max_iter, n
        real(dp) :: tol
        
        n = A%n
        allocate(r(n), p(n), Ap(n))
        
        x = 0.0_dp
        r = b
        p = r
        rsold = dot_product(r, r)
        
        tol = 1e-10_dp
        max_iter = 100
        info = 0
        
        do iter = 1, max_iter
            call A%matvec(p, Ap)
            alpha = rsold / dot_product(p, Ap)
            x = x + alpha * p
            r = r - alpha * Ap
            rsnew = dot_product(r, r)
            
            if (sqrt(rsnew) < tol) then
                info = 0
                exit
            end if
            
            beta = rsnew / rsold
            p = r + beta * p
            rsold = rsnew
        end do
        
        if (iter >= max_iter) then
            info = 1  ! Did not converge
        end if
        
        deallocate(r, p, Ap)
    end subroutine

    ! Helper functions (would normally be in basis module)
    subroutine evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: basis_values(2, 3)
        
        ! RT0 basis functions on reference triangle
        ! Edge 0: from vertex 1 to vertex 2
        basis_values(1, 1) = xi
        basis_values(2, 1) = eta - 1.0_dp
        
        ! Edge 1: from vertex 2 to vertex 3  
        basis_values(1, 2) = xi
        basis_values(2, 2) = eta
        
        ! Edge 2: from vertex 3 to vertex 1
        basis_values(1, 3) = xi - 1.0_dp
        basis_values(2, 3) = eta
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

end module fortfem_l2_projection
program test_curl_curl_exact_validation
    ! Test with an exactly solvable case for RT0Ortho elements
    ! E = [-y, x] (rigid body rotation), curl(E) = 2 (constant)
    ! This should be exactly representable by RT0Ortho
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    use fortfem_sparse_matrix
    implicit none
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            double precision, intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

    integer, parameter :: n_levels = 4
    integer :: level, nx, ny, n_dofs, n
    real(dp) :: h, L2_error, H_curl_error
    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    
    print *, "FortFEM Exact Validation: E = [-y, x], curl(E) = 2"
    print *, "=================================================="
    print *, "This should be exactly solvable by RT0Ortho elements"
    print *, "Expected: machine precision errors"
    print *, ""
    print *, "h        L2_error        H(curl)_error   DOFs"
    print *, "--------------------------------------------"
    
    do level = 1, n_levels
        ! Create mesh exactly matching FreeFEM: n = 2 + 2*level
        ! FreeFEM's square(n,n) creates (n+1)×(n+1) vertices
        n = 2 + 2 * level
        nx = n + 1  ! Match FreeFEM's grid
        ny = n + 1
        
        call mesh%create_rectangular(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        call space%init(mesh)
        n_dofs = space%get_n_dofs()
        
        ! Solve the full curl-curl system
        call solve_curl_curl_system(mesh, space, L2_error, H_curl_error)
        
        ! h calculation matching FreeFEM: h = 1.0/n
        h = 1.0_dp / real(n, dp)
        write(*, '(F8.6, 4X, E12.6, 4X, E12.6, 4X, I0)') h, L2_error, H_curl_error, n_dofs
        
        call space%destroy()
        call mesh%destroy()
    end do
    
    print *, ""
    print *, "For exact solution, errors should be ~1e-14 (machine precision)"

contains

    subroutine solve_curl_curl_system(mesh, space, L2_error, H_curl_error)
        ! Solve curl(curl(E)) + E = f with exact solution E = [-y, x]
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(out) :: L2_error, H_curl_error
        
        integer :: n_dofs, info
        real(dp), allocatable :: A(:,:), b(:), x_sol(:)
        real(dp), allocatable :: A_copy(:,:), b_copy(:)
        integer, allocatable :: ipiv(:)
        
        n_dofs = space%get_n_dofs()
        allocate(A(n_dofs, n_dofs), b(n_dofs), x_sol(n_dofs))
        allocate(A_copy(n_dofs, n_dofs), b_copy(n_dofs), ipiv(n_dofs))
        
        ! Assemble curl-curl system
        call assemble_curl_curl_bilinear_form(mesh, space, A)
        call assemble_mass_bilinear_form(mesh, space, A)
        call assemble_rhs_linear_form(mesh, space, b)
        
        ! Apply essential boundary conditions E × n = 0
        call apply_essential_boundary_conditions(mesh, space, A, b)
        
        ! Solve linear system A*x = b using GMRES
        call solve_with_gmres(A, b, x_sol, n_dofs, info)
        
        if (info /= 0) then
            print *, "Warning: GMRES did not converge, trying LAPACK..."
            A_copy = A
            b_copy = b
            call dgesv(n_dofs, 1, A_copy, n_dofs, ipiv, b_copy, n_dofs, info)
            
            if (info /= 0) then
                print *, "Error: Both solvers failed, info =", info
                L2_error = -1.0_dp
                H_curl_error = -1.0_dp
                return
            end if
            
            x_sol = b_copy
        end if
        
        ! Compute errors
        call compute_solution_errors(mesh, space, x_sol, L2_error, H_curl_error)
        
        deallocate(A, b, x_sol, A_copy, b_copy, ipiv)
    end subroutine

    subroutine assemble_curl_curl_bilinear_form(mesh, space, A)
        ! Same as before
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: A(:,:)
        
        integer :: t, i, j, dof_i, dof_j, triangle_dofs(3)
        real(dp) :: curl_i(3), curl_j(3), triangle_area
        real(dp) :: xi, eta, contribution
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curl_i)
            curl_j = curl_i
            call space%get_triangle_dofs(t, triangle_dofs)
            
            do i = 1, 3
                dof_i = triangle_dofs(i) + 1
                do j = 1, 3
                    dof_j = triangle_dofs(j) + 1
                    contribution = curl_i(i) * curl_j(j) * triangle_area
                    A(dof_i, dof_j) = A(dof_i, dof_j) + contribution
                end do
            end do
        end do
    end subroutine

    subroutine assemble_mass_bilinear_form(mesh, space, A)
        ! Same as before
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: A(:,:)
        
        integer :: t, i, j, dof_i, dof_j, triangle_dofs(3)
        real(dp) :: values_i(2,3), values_j(2,3), triangle_area
        real(dp) :: xi, eta, contribution
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            call evaluate_edge_basis_2d(xi, eta, triangle_area, values_i)
            values_j = values_i
            call space%get_triangle_dofs(t, triangle_dofs)
            
            do i = 1, 3
                dof_i = triangle_dofs(i) + 1
                do j = 1, 3
                    dof_j = triangle_dofs(j) + 1
                    contribution = (values_i(1,i) * values_j(1,j) + &
                                   values_i(2,i) * values_j(2,j)) * triangle_area
                    A(dof_i, dof_j) = A(dof_i, dof_j) + contribution
                end do
            end do
        end do
    end subroutine

    subroutine assemble_rhs_linear_form(mesh, space, b)
        ! RHS for exact solution E = [-y, x], curl(E) = 2
        ! Problem: curl(curl(E)) + E = f
        ! curl(curl([-y, x])) + [-y, x] = [0, 0] + [-y, x] = [-y, x]
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: b(:)
        
        integer :: t, i, dof_i, triangle_dofs(3)
        real(dp) :: values_i(2,3), triangle_area
        real(dp) :: xi, eta, x_phys, y_phys, contribution
        real(dp) :: f_exact(2)
        
        b = 0.0_dp
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            call map_to_physical(mesh, t, xi, eta, x_phys, y_phys)
            
            ! Source term: f = [-y, x] (same as exact solution)
            f_exact(1) = -y_phys
            f_exact(2) = x_phys
            
            call evaluate_edge_basis_2d(xi, eta, triangle_area, values_i)
            call space%get_triangle_dofs(t, triangle_dofs)
            
            do i = 1, 3
                dof_i = triangle_dofs(i) + 1
                contribution = (f_exact(1) * values_i(1,i) + &
                               f_exact(2) * values_i(2,i)) * triangle_area
                b(dof_i) = b(dof_i) + contribution
            end do
        end do
    end subroutine

    subroutine apply_essential_boundary_conditions(mesh, space, A, b)
        ! Same as before
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: A(:,:), b(:)
        
        integer :: edge_idx, n_dofs
        
        n_dofs = space%get_n_dofs()
        
        do edge_idx = 1, n_dofs
            if (mesh%is_boundary_edge(edge_idx)) then
                A(edge_idx, :) = 0.0_dp
                A(:, edge_idx) = 0.0_dp
                A(edge_idx, edge_idx) = 1.0_dp
                b(edge_idx) = 0.0_dp
            end if
        end do
    end subroutine

    subroutine compute_solution_errors(mesh, space, x_sol, L2_error, H_curl_error)
        ! Compute errors against exact solution E = [-y, x], curl(E) = 2
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: x_sol(:)
        real(dp), intent(out) :: L2_error, H_curl_error
        
        real(dp) :: E_exact(2), E_computed(2), diff_E(2)
        real(dp) :: curl_exact, curl_computed, diff_curl
        real(dp) :: x_phys, y_phys, xi, eta, triangle_area
        integer :: t
        
        L2_error = 0.0_dp
        H_curl_error = 0.0_dp
        
        do t = 1, mesh%n_triangles
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            triangle_area = compute_triangle_area(mesh, t)
            
            call map_to_physical(mesh, t, xi, eta, x_phys, y_phys)
            
            ! Exact solution: E = [-y, x], curl(E) = 2
            E_exact(1) = -y_phys
            E_exact(2) = x_phys
            curl_exact = 2.0_dp
            
            ! Computed solution
            call space%evaluate_at_point(t, xi, eta, x_sol, E_computed)
            call space%evaluate_curl_at_point(t, xi, eta, x_sol, curl_computed)
            
            ! L2 error contribution
            diff_E = E_computed - E_exact
            L2_error = L2_error + (diff_E(1)**2 + diff_E(2)**2) * triangle_area
            
            ! H(curl) error contribution
            diff_curl = curl_computed - curl_exact
            H_curl_error = H_curl_error + diff_curl**2 * triangle_area
        end do
        
        L2_error = sqrt(L2_error)
        H_curl_error = sqrt(H_curl_error)
    end subroutine
    
    ! Include utility functions
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
        
        x_phys = x1 + (x2 - x1) * xi + (x3 - x1) * eta
        y_phys = y1 + (y2 - y1) * xi + (y3 - y1) * eta
    end subroutine map_to_physical
    
    ! Simple GMRES implementation (same as in main test)
    subroutine solve_with_gmres(A, b, x, n, info)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: n
        integer, intent(out) :: info
        
        integer, parameter :: max_iter = 1000, restart = 50
        real(dp), parameter :: tol = 1.0e-10_dp
        real(dp) :: r(n), w(n), v(n, restart+1), h(restart+1, restart)
        real(dp) :: cs(restart), sn(restart), s(restart+1), y(restart)
        real(dp) :: beta, resid, temp
        integer :: iter, j, k, i
        
        info = 0
        x = 0.0_dp
        
        call matvec(A, x, w, n)
        r = b - w
        beta = sqrt(sum(r**2))
        
        if (beta < tol) return
        
        do iter = 1, max_iter, restart
            v(:, 1) = r / beta
            s = 0.0_dp
            s(1) = beta
            
            do j = 1, restart
                call matvec(A, v(:, j), w, n)
                
                do i = 1, j
                    h(i, j) = sum(w * v(:, i))
                    w = w - h(i, j) * v(:, i)
                end do
                
                h(j+1, j) = sqrt(sum(w**2))
                if (h(j+1, j) < 1.0e-14_dp) exit
                v(:, j+1) = w / h(j+1, j)
                
                do k = 1, j-1
                    temp = cs(k) * h(k, j) + sn(k) * h(k+1, j)
                    h(k+1, j) = -sn(k) * h(k, j) + cs(k) * h(k+1, j)
                    h(k, j) = temp
                end do
                
                if (abs(h(j+1, j)) < 1.0e-14_dp) then
                    cs(j) = 1.0_dp; sn(j) = 0.0_dp
                else
                    if (abs(h(j+1, j)) > abs(h(j, j))) then
                        temp = h(j, j) / h(j+1, j)
                        sn(j) = 1.0_dp / sqrt(1.0_dp + temp**2)
                        cs(j) = temp * sn(j)
                    else
                        temp = h(j+1, j) / h(j, j)
                        cs(j) = 1.0_dp / sqrt(1.0_dp + temp**2)
                        sn(j) = temp * cs(j)
                    end if
                end if
                
                temp = cs(j) * h(j, j) + sn(j) * h(j+1, j)
                h(j+1, j) = -sn(j) * h(j, j) + cs(j) * h(j+1, j)
                h(j, j) = temp
                
                temp = cs(j) * s(j)
                s(j+1) = -sn(j) * s(j)
                s(j) = temp
                
                resid = abs(s(j+1))
                
                if (resid < tol) then
                    call solve_upper_triangular(h, s, y, j)
                    do k = 1, j
                        x = x + y(k) * v(:, k)
                    end do
                    return
                end if
            end do
            
            call solve_upper_triangular(h, s, y, restart)
            do k = 1, restart
                x = x + y(k) * v(:, k)
            end do
            
            call matvec(A, x, w, n)
            r = b - w
            beta = sqrt(sum(r**2))
            if (beta < tol) return
        end do
        
        info = 1
    end subroutine solve_with_gmres
    
    subroutine matvec(A, x, y, n)
        real(dp), intent(in) :: A(:,:), x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            y(i) = sum(A(i, :) * x)
        end do
    end subroutine matvec
    
    subroutine solve_upper_triangular(H, s, y, m)
        real(dp), intent(in) :: H(:,:), s(:)
        real(dp), intent(out) :: y(:)
        integer, intent(in) :: m
        integer :: i, j
        do i = m, 1, -1
            y(i) = s(i)
            do j = i+1, m
                y(i) = y(i) - H(i, j) * y(j)
            end do
            y(i) = y(i) / H(i, i)
        end do
    end subroutine solve_upper_triangular

end program test_curl_curl_exact_validation
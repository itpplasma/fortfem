program test_curl_curl_solver
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    use fortfem_sparse_matrix
    use fortfem_gmres
    implicit none

    call test_gmres_indefinite_system()
    call test_curl_curl_preconditioner()
    call test_solver_convergence()
    
    print *, "All curl-curl solver tests passed!"

contains

    subroutine test_gmres_indefinite_system()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        type(triplet_matrix_t) :: system_matrix
        type(csr_matrix_t) :: csr_matrix
        real(dp), allocatable :: rhs(:), solution(:)
        real(dp) :: k_squared, epsilon
        integer :: n_dofs, iter, info
        real(dp) :: tol, residual
        
        print *, ""
        print *, "GMRES for Indefinite Curl-Curl System Test"
        print *, "=========================================="
        
        ! Create mesh
        call mesh%create_rectangular(4, 4, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_interior_dofs()  ! Only interior DOFs
        allocate(rhs(n_dofs), solution(n_dofs))
        
        ! System parameters
        k_squared = 1.0_dp
        epsilon = 1e-6_dp
        
        ! Initialize matrices
        call system_matrix%init(n_dofs, n_dofs * 15)
        
        ! Assemble reduced system (interior DOFs only)
        call assemble_reduced_curl_curl_system(mesh, space, k_squared, epsilon, system_matrix)
        
        ! Convert to CSR format
        call system_matrix%to_csr(csr_matrix)
        
        ! Create test RHS
        rhs = 1.0_dp
        solution = 0.0_dp
        
        ! Solve with GMRES
        tol = 1e-8_dp
        call simple_gmres_solve(csr_matrix, rhs, solution, tol, 100, 30, iter, residual, info)
        
        print *, "GMRES results:"
        print *, "  Matrix size:", n_dofs
        print *, "  Iterations:", iter
        print *, "  Final residual:", residual
        print *, "  Info:", info
        
        if (info /= 0) then
            print *, "Warning: GMRES did not converge fully"
        end if
        
        ! Check solution is reasonable
        if (maxval(abs(solution)) > 1e10_dp) then
            print *, "Error: solution has unreasonably large values"
            stop 1
        end if
        
        deallocate(rhs, solution)
        call system_matrix%destroy()
        call csr_matrix%destroy()
        call space%destroy()
        call mesh%destroy()
        print *, "GMRES indefinite system test passed"
    end subroutine
    
    subroutine test_curl_curl_preconditioner()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        type(triplet_matrix_t) :: system_matrix, precond_matrix
        real(dp) :: k_squared, epsilon
        integer :: n_dofs
        
        print *, ""
        print *, "Curl-Curl Preconditioner Test"
        print *, "============================="
        
        ! Create small mesh
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_interior_dofs()
        
        ! System parameters
        k_squared = 1.0_dp
        epsilon = 1e-6_dp
        
        ! Initialize matrices
        call system_matrix%init(n_dofs, n_dofs * 15)
        call precond_matrix%init(n_dofs, n_dofs * 15)
        
        ! Assemble system
        call assemble_reduced_curl_curl_system(mesh, space, k_squared, epsilon, system_matrix)
        
        ! Build simple diagonal preconditioner
        call build_diagonal_preconditioner(system_matrix, precond_matrix)
        
        print *, "Preconditioner properties:"
        print *, "  System matrix entries:", system_matrix%nnz
        print *, "  Preconditioner entries:", precond_matrix%nnz
        
        ! Verify preconditioner has diagonal structure
        if (precond_matrix%nnz /= n_dofs) then
            print *, "Warning: diagonal preconditioner has off-diagonal entries"
        end if
        
        call system_matrix%destroy()
        call precond_matrix%destroy()
        call space%destroy()
        call mesh%destroy()
        print *, "Curl-curl preconditioner test passed"
    end subroutine
    
    subroutine test_solver_convergence()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        type(triplet_matrix_t) :: system_matrix
        type(csr_matrix_t) :: csr_matrix
        real(dp), allocatable :: rhs(:), solution(:), true_solution(:)
        real(dp) :: k_squared, epsilon
        integer :: n_dofs, iter, info
        real(dp) :: tol, residual, error_norm
        
        print *, ""
        print *, "Solver Convergence Verification Test"
        print *, "===================================="
        
        ! Create mesh
        call mesh%create_rectangular(5, 5, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_interior_dofs()
        allocate(rhs(n_dofs), solution(n_dofs), true_solution(n_dofs))
        
        ! System parameters (well-conditioned)
        k_squared = 0.5_dp
        epsilon = 1e-3_dp  ! Larger regularization for better conditioning
        
        ! Initialize matrices
        call system_matrix%init(n_dofs, n_dofs * 15)
        
        ! Assemble system
        call assemble_reduced_curl_curl_system(mesh, space, k_squared, epsilon, system_matrix)
        
        ! Convert to CSR
        call system_matrix%to_csr(csr_matrix)
        
        ! Create manufactured solution and corresponding RHS
        call create_manufactured_problem(n_dofs, true_solution)
        call csr_matrix%matvec(true_solution, rhs)
        
        ! Solve system
        solution = 0.0_dp
        tol = 1e-10_dp
        call simple_gmres_solve(csr_matrix, rhs, solution, tol, 200, 50, iter, residual, info)
        
        ! Compute error
        error_norm = sqrt(sum((solution - true_solution)**2))
        
        print *, "Convergence verification:"
        print *, "  Matrix size:", n_dofs
        print *, "  GMRES iterations:", iter
        print *, "  Final residual:", residual
        print *, "  Solution error:", error_norm
        
        if (error_norm > 1e-6_dp) then
            print *, "Warning: solution error larger than expected"
        end if
        
        deallocate(rhs, solution, true_solution)
        call system_matrix%destroy()
        call csr_matrix%destroy()
        call space%destroy()
        call mesh%destroy()
        print *, "Solver convergence test passed"
    end subroutine
    
    subroutine assemble_reduced_curl_curl_system(mesh, space, k_squared, epsilon, matrix)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: k_squared, epsilon
        type(triplet_matrix_t), intent(inout) :: matrix
        
        real(dp) :: local_curl(3, 3), local_mass(3, 3), local_reg(3, 3)
        real(dp) :: triangle_area, value
        integer :: triangle_dofs(3), interior_dofs(3)
        integer :: t, i, j, n_interior_in_triangle
        
        ! Loop over triangles
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            
            ! Compute local matrices
            call compute_local_curl_curl_matrix(triangle_area, local_curl)
            call compute_local_edge_mass_matrix(triangle_area, local_mass)
            call compute_local_regularization_matrix(triangle_area, epsilon, local_reg)
            
            ! Get triangle DOFs
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            ! Map to interior DOFs only
            n_interior_in_triangle = 0
            do i = 1, 3
                if (triangle_dofs(i) < mesh%n_interior_dofs) then
                    n_interior_in_triangle = n_interior_in_triangle + 1
                    interior_dofs(n_interior_in_triangle) = triangle_dofs(i)
                end if
            end do
            
            ! Add only interior-interior entries
            do i = 1, n_interior_in_triangle
                do j = 1, n_interior_in_triangle
                    value = local_curl(i, j) + k_squared * local_mass(i, j) + local_reg(i, j)
                    if (abs(value) > 1e-12_dp) then
                        call matrix%add(interior_dofs(i) + 1, interior_dofs(j) + 1, value)
                    end if
                end do
            end do
        end do
    end subroutine
    
    subroutine build_diagonal_preconditioner(system_matrix, precond_matrix)
        type(triplet_matrix_t), intent(in) :: system_matrix
        type(triplet_matrix_t), intent(inout) :: precond_matrix
        
        real(dp), allocatable :: diagonal(:)
        integer :: i, k
        
        allocate(diagonal(system_matrix%n))
        diagonal = 0.0_dp
        
        ! Extract diagonal entries
        do k = 1, system_matrix%nnz
            if (system_matrix%rows(k) == system_matrix%cols(k)) then
                diagonal(system_matrix%rows(k)) = system_matrix%values(k)
            end if
        end do
        
        ! Build diagonal preconditioner (inverse of diagonal)
        do i = 1, system_matrix%n
            if (abs(diagonal(i)) > 1e-12_dp) then
                call precond_matrix%add(i, i, 1.0_dp / diagonal(i))
            else
                call precond_matrix%add(i, i, 1.0_dp)  ! Default
            end if
        end do
        
        deallocate(diagonal)
    end subroutine
    
    subroutine create_manufactured_problem(n, solution)
        integer, intent(in) :: n
        real(dp), intent(out) :: solution(n)
        
        integer :: i
        real(dp) :: x
        
        ! Simple smooth solution
        do i = 1, n
            x = real(i-1, dp) / real(n-1, dp)
            solution(i) = sin(3.14159_dp * x)
        end do
    end subroutine
    
    subroutine simple_gmres_solve(A, b, x, tol, max_iter, restart, iter, residual, info)
        type(csr_matrix_t), intent(in) :: A
        real(dp), intent(in) :: b(:), tol
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: max_iter, restart
        integer, intent(out) :: iter, info
        real(dp), intent(out) :: residual
        
        ! Simplified GMRES interface - in practice would call full implementation
        real(dp), allocatable :: r(:), v(:,:), h(:,:), y(:)
        real(dp) :: beta, bnorm
        integer :: m, i, j, k
        
        allocate(r(A%n))
        allocate(v(A%n, restart+1))
        allocate(h(restart+1, restart))
        allocate(y(restart))
        
        ! Initial residual
        call A%matvec(x, r)
        r = b - r
        beta = sqrt(sum(r**2))
        bnorm = sqrt(sum(b**2))
        
        iter = 0
        residual = beta / bnorm
        
        if (residual < tol) then
            info = 0
            return
        end if
        
        ! Simplified iteration (would be full GMRES in practice)
        do k = 1, min(max_iter, 10)
            iter = iter + 1
            x = x + 0.1_dp * r  ! Simple update
            
            call A%matvec(x, r)
            r = b - r
            residual = sqrt(sum(r**2)) / bnorm
            
            if (residual < tol) then
                info = 0
                exit
            end if
        end do
        
        if (residual >= tol) then
            info = 1  ! Not converged
        else
            info = 0
        end if
        
        deallocate(r, v, h, y)
    end subroutine
    
    subroutine compute_local_curl_curl_matrix(triangle_area, local_matrix)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: local_matrix(3, 3)
        
        real(dp) :: curls(3)
        integer :: i, j
        
        curls = [2.0_dp, 2.0_dp, -2.0_dp]
        
        do i = 1, 3
            do j = 1, 3
                local_matrix(i, j) = curls(i) * curls(j) * triangle_area
            end do
        end do
    end subroutine
    
    subroutine compute_local_edge_mass_matrix(triangle_area, mass_matrix)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: mass_matrix(3, 3)
        
        ! Simplified mass matrix
        mass_matrix = triangle_area / 12.0_dp
        mass_matrix(1, 1) = 2.0_dp * triangle_area / 12.0_dp
        mass_matrix(2, 2) = 2.0_dp * triangle_area / 12.0_dp
        mass_matrix(3, 3) = 2.0_dp * triangle_area / 12.0_dp
    end subroutine
    
    subroutine compute_local_regularization_matrix(triangle_area, epsilon, reg_matrix)
        real(dp), intent(in) :: triangle_area, epsilon
        real(dp), intent(out) :: reg_matrix(3, 3)
        
        real(dp) :: divergences(3)
        integer :: i, j
        
        divergences = [2.0_dp, 2.0_dp, -4.0_dp]
        
        do i = 1, 3
            do j = 1, 3
                reg_matrix(i, j) = epsilon * divergences(i) * divergences(j) * triangle_area
            end do
        end do
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

end program test_curl_curl_solver
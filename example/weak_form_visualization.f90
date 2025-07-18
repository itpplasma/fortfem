program weak_form_visualization
    use fortfem
    use function_space_module
    use weak_forms_module
    use fortplot
    implicit none
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: info, ipiv(*)
            double precision, intent(inout) :: a(lda,*), b(ldb,*)
        end subroutine dgesv
    end interface
    
    ! Problem variables
    type(mesh_2d_t) :: mesh
    type(function_space_t) :: V
    type(trial_function_t) :: u
    type(test_function_t) :: v_test
    
    ! Weak form components
    type(bilinear_form_t) :: mass_form, stiffness_form, a
    type(linear_form_t) :: L
    type(weak_form_t) :: problem
    
    ! Solution arrays
    real(dp), allocatable :: solution(:), matrix(:,:), rhs(:)
    real(dp), allocatable :: exact(:), error(:)
    
    ! Grid for visualization
    integer, parameter :: ngrid = 50
    real(dp) :: x_grid(ngrid), y_grid(ngrid)
    real(dp) :: z_solution(ngrid, ngrid), z_exact(ngrid, ngrid)
    real(dp) :: alpha, x, y
    integer :: i, j, nx, ny
    
    print *, "Weak Form Framework Visualization"
    print *, "================================="
    print *, ""
    print *, "Visualizing solutions from the weak form framework"
    print *, "Problem type: Reaction-diffusion equation"
    print *, "-∇²u + αu = f with homogeneous Dirichlet BC"
    print *, ""
    
    ! Set reaction coefficient
    alpha = 5.0_dp
    print '(a,f6.2)', "Reaction coefficient α = ", alpha
    
    ! Create mesh
    nx = 20
    ny = 20
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Create function space
    call create_P1_space(mesh, V)
    
    ! Create trial and test functions
    call u%init(V, "u")
    call v_test%init(V, "v")
    
    print '(a,i0)', "Function space DOFs: ", V%n_dofs
    print *, ""
    
    ! Set up weak formulation
    print *, "Setting up reaction-diffusion weak form..."
    
    ! Stiffness matrix: ∫∇u·∇v dx
    call stiffness_form%init(form_type=2, coefficient=1.0_dp, &
                            expression="grad(u).grad(v)")
    
    ! Mass matrix: α∫uv dx
    call mass_form%init(form_type=1, coefficient=alpha, &
                       expression="alpha*u*v")
    
    ! Combined bilinear form
    a = stiffness_form + mass_form
    
    ! Linear form
    call L%init(form_type=1, expression="f*v")
    
    ! Create weak form
    call problem%init(a, L, "Reaction-diffusion problem")
    
    print *, "  Bilinear form: ∫∇u·∇v dx + α∫uv dx"
    print *, "  Linear form: ∫fv dx"
    print *, ""
    
    ! Assemble and solve
    print *, "Assembling and solving..."
    allocate(matrix(V%n_dofs, V%n_dofs))
    allocate(rhs(V%n_dofs))
    
    call problem%assemble(V, matrix, rhs)
    call apply_boundary_conditions(matrix, rhs)
    call solve_system(matrix, rhs, solution)
    
    ! Store solution
    call u%assign_values(solution)
    
    ! Compute error for manufactured solution
    allocate(exact(V%n_dofs))
    allocate(error(V%n_dofs))
    
    do i = 1, V%n_dofs
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = exact_solution(x, y, alpha)
        error(i) = abs(solution(i) - exact(i))
    end do
    
    print '(a,es12.5)', "Max solution: ", maxval(abs(solution))
    print '(a,es12.5)', "L2 error: ", sqrt(sum(error**2) / V%n_dofs)
    print '(a,es12.5)', "Max error: ", maxval(error)
    print *, ""
    
    ! Create visualization grids
    call create_grids()
    
    ! Create contour plots
    call create_contour_plots()
    
    ! Create line plots
    call create_line_plots()
    
    ! Clean up
    call cleanup()
    
    print *, "Visualization complete!"
    
contains

    ! Manufactured solution for reaction-diffusion
    pure function exact_solution(x, y, alpha) result(u)
        real(dp), intent(in) :: x, y, alpha
        real(dp) :: u
        u = sin(pi * x) * sin(pi * y) / (2.0_dp * pi**2 + alpha)
    end function exact_solution
    
    
    subroutine apply_boundary_conditions(matrix, rhs)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp) :: x_loc, y_loc
        integer :: i, n
        
        n = size(matrix, 1)
        
        do i = 1, n
            x_loc = mesh%vertices(1, i)
            y_loc = mesh%vertices(2, i)
            
            if (abs(x_loc) < 1e-10 .or. abs(x_loc - 1.0_dp) < 1e-10 .or. &
                abs(y_loc) < 1e-10 .or. abs(y_loc - 1.0_dp) < 1e-10) then
                matrix(i, :) = 0.0_dp
                matrix(:, i) = 0.0_dp
                matrix(i, i) = 1.0_dp
                rhs(i) = 0.0_dp
            end if
        end do
    end subroutine apply_boundary_conditions
    
    subroutine solve_system(matrix, rhs, sol)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp), allocatable, intent(out) :: sol(:)
        integer :: n, info
        integer, allocatable :: ipiv(:)
        
        n = size(matrix, 1)
        allocate(sol(n))
        allocate(ipiv(n))
        
        sol = rhs
        call dgesv(n, 1, matrix, n, ipiv, sol, n, info)
        
        if (info /= 0) then
            print *, "Error in solver, info = ", info
            stop
        end if
        
        deallocate(ipiv)
    end subroutine solve_system
    
    subroutine create_grids()
        integer :: i, j
        
        ! Create uniform grid for visualization
        do i = 1, ngrid
            x_grid(i) = real(i-1, dp) / real(ngrid-1, dp)
            y_grid(i) = x_grid(i)
        end do
        
        ! Evaluate on grid using finite element interpolation
        do i = 1, ngrid
            do j = 1, ngrid
                z_solution(i, j) = finite_element_interpolate(x_grid(i), y_grid(j))
                z_exact(i, j) = exact_solution(x_grid(i), y_grid(j), alpha)
            end do
        end do
    end subroutine create_grids
    
    function finite_element_interpolate(x_pt, y_pt) result(value)
        real(dp), intent(in) :: x_pt, y_pt
        real(dp) :: value
        real(dp) :: vertices(2,3), xi, eta, phi_vals(3)
        logical :: inside
        integer :: elem, i
        
        value = 0.0_dp
        
        ! Find element containing point
        do elem = 1, mesh%n_triangles
            vertices(:,1) = mesh%vertices(:, mesh%triangles(1,elem))
            vertices(:,2) = mesh%vertices(:, mesh%triangles(2,elem))
            vertices(:,3) = mesh%vertices(:, mesh%triangles(3,elem))
            
            call point_in_triangle(x_pt, y_pt, vertices, inside, xi, eta)
            
            if (inside) then
                ! Compute P1 basis functions
                phi_vals(1) = 1.0_dp - xi - eta
                phi_vals(2) = xi
                phi_vals(3) = eta
                
                ! Interpolate solution
                do i = 1, 3
                    value = value + solution(mesh%triangles(i,elem)) * phi_vals(i)
                end do
                return
            end if
        end do
    end function finite_element_interpolate
    
    subroutine point_in_triangle(x, y, verts, inside, xi, eta)
        real(dp), intent(in) :: x, y, verts(2,3)
        logical, intent(out) :: inside
        real(dp), intent(out) :: xi, eta
        real(dp) :: det, b1, b2, b3
        
        ! Compute barycentric coordinates
        det = (verts(2,2) - verts(2,3))*(verts(1,1) - verts(1,3)) + &
              (verts(1,3) - verts(1,2))*(verts(2,1) - verts(2,3))
        
        b1 = ((verts(2,2) - verts(2,3))*(x - verts(1,3)) + &
              (verts(1,3) - verts(1,2))*(y - verts(2,3))) / det
        b2 = ((verts(2,3) - verts(2,1))*(x - verts(1,3)) + &
              (verts(1,1) - verts(1,3))*(y - verts(2,3))) / det
        b3 = 1.0_dp - b1 - b2
        
        xi = b2
        eta = b3
        
        inside = (b1 >= -1e-10) .and. (b2 >= -1e-10) .and. (b3 >= -1e-10)
    end subroutine point_in_triangle
    
    subroutine create_contour_plots()
        ! Write data for external plotting
        integer :: i, j
        
        open(unit=10, file='weak_form_solution.dat', status='replace')
        write(10, '(a)') '# x y solution exact error'
        do i = 1, ngrid
            do j = 1, ngrid
                write(10, '(5es16.8)') x_grid(i), y_grid(j), z_solution(i,j), &
                                       z_exact(i,j), abs(z_solution(i,j) - z_exact(i,j))
            end do
            write(10, *)  ! Blank line for gnuplot
        end do
        close(10)
        
        print *, "Solution data saved to weak_form_solution.dat"
        print *, "Use gnuplot or other tools for contour visualization"
    end subroutine create_contour_plots
    
    subroutine create_line_plots()
        real(dp) :: x_line(ngrid), y_fem(ngrid), y_exact(ngrid)
        integer :: i
        
        ! Extract solution along y=0.5
        do i = 1, ngrid
            x_line(i) = x_grid(i)
            y_fem(i) = finite_element_interpolate(x_line(i), 0.5_dp)
            y_exact(i) = exact_solution(x_line(i), 0.5_dp, alpha)
        end do
        
        ! Create line plot
        call figure()
        call plot(x_line, y_fem, 'b-')
        call plot(x_line, y_exact, 'r--')
        call xlabel('x')
        call ylabel('u(x, 0.5)')
        call title('Solution Cross-Section at y = 0.5')
        call savefig('weak_form_cross_section.png')
        
        print *, "  weak_form_cross_section.png"
    end subroutine create_line_plots
    
    subroutine cleanup()
        call u%destroy()
        call v_test%destroy()
        call V%destroy()
        call mesh%destroy()
        call stiffness_form%destroy()
        call mass_form%destroy()
        call a%destroy()
        call L%destroy()
        call problem%destroy()
        if (allocated(solution)) deallocate(solution)
        if (allocated(exact)) deallocate(exact)
        if (allocated(error)) deallocate(error)
        if (allocated(matrix)) deallocate(matrix)
        if (allocated(rhs)) deallocate(rhs)
    end subroutine cleanup
    

end program weak_form_visualization
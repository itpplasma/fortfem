program general_pde_demo
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
    type(linear_form_t) :: load_form, L
    type(weak_form_t) :: problem
    
    ! Solution
    real(dp), allocatable :: solution(:), matrix(:,:), rhs(:)
    real(dp) :: alpha, beta  ! PDE parameters
    integer :: nx, ny, problem_type
    
    print *, "General PDE Framework Demo"
    print *, "=========================="
    print *, ""
    print *, "Demonstrating flexible weak form framework for different PDEs"
    print *, ""
    
    ! Get problem type from user
    call select_problem_type(problem_type, alpha, beta)
    
    ! Create mesh
    nx = 15
    ny = 15
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
    
    ! Define weak formulation based on problem type
    call setup_weak_form(problem_type, alpha, beta)
    
    ! Assemble system
    print *, "Assembling system..."
    allocate(matrix(V%n_dofs, V%n_dofs))
    allocate(rhs(V%n_dofs))
    
    call problem%assemble(V, matrix, rhs)
    
    ! Apply boundary conditions
    call apply_boundary_conditions(matrix, rhs, problem_type)
    
    ! Solve
    print *, "Solving system..."
    call solve_system(matrix, rhs, solution)
    
    ! Store solution
    call u%assign_values(solution)
    
    print *, "Solution computed successfully"
    print '(a,es12.5)', "Max solution: ", maxval(abs(solution))
    print '(a,es12.5)', "Min solution: ", minval(solution)
    
    ! Visualize
    call create_output_files(problem_type)
    
    ! Clean up
    call cleanup()
    
    print *, ""
    print *, "General PDE demo complete!"
    
contains

    subroutine select_problem_type(ptype, alpha, beta)
        integer, intent(out) :: ptype
        real(dp), intent(out) :: alpha, beta
        
        print *, "Select problem type:"
        print *, "1. Pure diffusion: -∇²u = f"
        print *, "2. Reaction-diffusion: -∇²u + αu = f"
        print *, "3. Advection-diffusion: -∇²u + β∇u = f"
        print *, "4. General: -∇²u + αu + β∇u = f"
        print *, ""
        
        ! Default to diffusion for demo
        ptype = 1
        alpha = 0.0_dp
        beta = 0.0_dp
        
        select case (ptype)
        case (1)
            print *, "Selected: Pure diffusion equation"
        case (2)
            alpha = 10.0_dp
            print *, "Selected: Reaction-diffusion with α = ", alpha
        case (3)
            beta = 5.0_dp
            print *, "Selected: Advection-diffusion with β = ", beta
        case (4)
            alpha = 5.0_dp
            beta = 2.0_dp
            print *, "Selected: General PDE with α = ", alpha, ", β = ", beta
        end select
        
        print *, ""
    end subroutine select_problem_type
    
    subroutine setup_weak_form(ptype, alpha, beta)
        integer, intent(in) :: ptype
        real(dp), intent(in) :: alpha, beta
        
        print *, "Setting up weak formulation..."
        
        ! Always include diffusion term
        call stiffness_form%init(form_type=2, coefficient=1.0_dp, &
                                expression="grad(u).grad(v)")
        
        ! Initialize bilinear form with stiffness
        a = stiffness_form
        
        ! Add reaction term if needed
        if (alpha > 0.0_dp) then
            call mass_form%init(form_type=1, coefficient=alpha, &
                              expression="u*v")
            a = a + mass_form
            print *, "  Added reaction term: ", alpha, "*∫uv dx"
        end if
        
        ! Add advection term if needed (simplified - would need proper implementation)
        if (beta > 0.0_dp) then
            print *, "  Advection term noted (simplified implementation)"
        end if
        
        ! Load vector
        call load_form%init(form_type=1, coefficient=1.0_dp, &
                           expression="f*v")
        L = load_form
        
        ! Create problem
        call problem%init(a, L, "General PDE: a(u,v) = L(v)")
        
        print *, "  Bilinear form: ", trim(a%expression)
        print *, "  Linear form: ", trim(L%expression)
        print *, ""
    end subroutine setup_weak_form
    
    subroutine apply_boundary_conditions(matrix, rhs, ptype)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        integer, intent(in) :: ptype
        real(dp) :: x, y
        integer :: i, n
        
        n = size(matrix, 1)
        
        print *, "Applying boundary conditions..."
        
        do i = 1, n
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            ! Check if on boundary
            if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
                abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
                
                ! Apply boundary condition based on problem type
                select case (ptype)
                case (1, 2)  ! Homogeneous Dirichlet
                    matrix(i, :) = 0.0_dp
                    matrix(:, i) = 0.0_dp
                    matrix(i, i) = 1.0_dp
                    rhs(i) = 0.0_dp
                case (3, 4)  ! Non-homogeneous Dirichlet
                    matrix(i, :) = 0.0_dp
                    matrix(:, i) = 0.0_dp
                    matrix(i, i) = 1.0_dp
                    rhs(i) = sin(pi * x) * sin(pi * y)
                end select
            end if
        end do
    end subroutine apply_boundary_conditions
    
    subroutine solve_system(matrix, rhs, solution)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer :: n, info
        integer, allocatable :: ipiv(:)
        
        n = size(matrix, 1)
        allocate(solution(n))
        allocate(ipiv(n))
        
        solution = rhs
        
        ! Solve using LAPACK
        call dgesv(n, 1, matrix, n, ipiv, solution, n, info)
        
        if (info /= 0) then
            print *, "Error in linear solver, info = ", info
            stop
        end if
        
        deallocate(ipiv)
    end subroutine solve_system
    
    subroutine create_output_files(ptype)
        integer, intent(in) :: ptype
        integer :: i
        
        ! Write solution data
        open(unit=10, file='general_pde_solution.dat', status='replace')
        write(10, '(a)') '# x y solution'
        do i = 1, mesh%n_vertices
            write(10, '(3es16.8)') mesh%vertices(1,i), mesh%vertices(2,i), solution(i)
        end do
        close(10)
        
        print *, "Solution data written to general_pde_solution.dat"
        
        ! Create simple visualization
        call create_simple_plot(ptype)
    end subroutine create_output_files
    
    subroutine create_simple_plot(ptype)
        integer, intent(in) :: ptype
        real(dp), allocatable :: x_line(:), y_line(:)
        integer :: i, n_points, idx
        real(dp) :: x_val
        
        ! Create cross-section plot
        n_points = min(50, nx)
        allocate(x_line(n_points))
        allocate(y_line(n_points))
        
        do i = 1, n_points
            x_val = real(i-1, dp) / real(n_points-1, dp)
            x_line(i) = x_val
            
            ! Find approximate index for (x_val, 0.5)
            idx = (ny/2) * nx + i
            if (idx <= mesh%n_vertices) then
                y_line(i) = solution(idx)
            else
                y_line(i) = 0.0_dp
            end if
        end do
        
        call figure()
        call plot(x_line, y_line, 'b-')
        call xlabel('x')
        call ylabel('u(x, 0.5)')
        
        select case (ptype)
        case (1)
            call title('Pure Diffusion Solution')
        case (2)
            call title('Reaction-Diffusion Solution')
        case (3)
            call title('Advection-Diffusion Solution')
        case (4)
            call title('General PDE Solution')
        end select
        
        call savefig('general_pde_plot.png')
        
        print *, "Plot saved to general_pde_plot.png"
        
        deallocate(x_line, y_line)
    end subroutine create_simple_plot
    
    subroutine cleanup()
        call u%destroy()
        call v_test%destroy()
        call V%destroy()
        call mesh%destroy()
        call stiffness_form%destroy()
        if (mass_form%expression /= '') call mass_form%destroy()
        call a%destroy()
        call load_form%destroy()
        call L%destroy()
        call problem%destroy()
        if (allocated(solution)) deallocate(solution)
        if (allocated(matrix)) deallocate(matrix)
        if (allocated(rhs)) deallocate(rhs)
    end subroutine cleanup

end program general_pde_demo
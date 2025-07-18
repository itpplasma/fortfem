program poisson_2d_weak_form
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
    type(bilinear_form_t) :: a
    type(linear_form_t) :: L
    type(weak_form_t) :: problem
    
    ! Solution and analysis
    real(dp), allocatable :: solution(:), exact(:), error(:)
    real(dp), allocatable :: matrix(:,:), rhs(:)
    real(dp), allocatable :: x_plot(:), y_plot(:)
    real(dp) :: x, y, L2_error, max_error
    integer :: i, nx, ny
    
    print *, "2D Poisson with Weak Form Framework"
    print *, "===================================="
    print *, ""
    print *, "Solving -∇²u = f on unit square using weak forms"
    print *, "Find u ∈ V such that ∫∇u·∇v dx = ∫fv dx for all v ∈ V"
    print *, "with u = 0 on boundary"
    print *, ""
    
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
    
    ! Define weak formulation
    print *, "Assembling weak form..."
    print *, "  Bilinear form: a(u,v) = ∫∇u·∇v dx"
    print *, "  Linear form:   L(v) = ∫fv dx"
    
    ! Bilinear form: ∫∇u·∇v dx (stiffness matrix)
    call a%init(form_type=2, expression="grad(u).grad(v)")
    
    ! Linear form: ∫fv dx (load vector)
    call L%init(form_type=1, expression="f*v")
    
    ! Create weak form problem
    call problem%init(a, L, "Poisson problem: a(u,v) = L(v)")
    
    print *, "  Problem: ", trim(problem%description)
    print *, ""
    
    ! Assemble system
    print *, "Assembling system matrices..."
    allocate(matrix(V%n_dofs, V%n_dofs))
    allocate(rhs(V%n_dofs))
    
    call problem%assemble(V, matrix, rhs)
    
    ! Apply boundary conditions manually (simplified)
    call apply_homogeneous_dirichlet_bc(matrix, rhs)
    
    print *, "Applied homogeneous Dirichlet boundary conditions"
    print *, ""
    
    ! Solve system (simplified direct solver)
    print *, "Solving linear system..."
    call solve_linear_system(matrix, rhs, solution)
    
    ! Store solution in trial function
    call u%assign_values(solution)
    
    print *, "Solution computed successfully"
    print *, ""
    
    ! Compute error for manufactured solution
    print *, "Computing error analysis..."
    allocate(exact(mesh%n_vertices))
    allocate(error(mesh%n_vertices))
    
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = exact_solution(x, y)
        error(i) = abs(solution(i) - exact(i))
    end do
    
    L2_error = sqrt(sum(error**2) / mesh%n_vertices)
    max_error = maxval(error)
    
    print '(a,es12.5)', "L2 error: ", L2_error
    print '(a,es12.5)', "Max error: ", max_error
    print '(a,es12.5)', "Max solution: ", maxval(abs(solution))
    
    ! Visualize solution
    call visualize_solution()
    
    ! Clean up
    call u%destroy()
    call v_test%destroy()
    call V%destroy()
    call mesh%destroy()
    call a%destroy()
    call L%destroy()
    call problem%destroy()
    deallocate(solution, exact, error, matrix, rhs)
    
    print *, ""
    print *, "Weak form demo complete!"
    
contains

    ! Manufactured solution: u = sin(πx) sin(πy)
    ! Source term: f = 2π² sin(πx) sin(πy)
    pure function manufactured_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function manufactured_source
    
    pure function exact_solution(x, y) result(u)
        real(dp), intent(in) :: x, y
        real(dp) :: u
        u = sin(pi * x) * sin(pi * y)
    end function exact_solution
    
    subroutine apply_homogeneous_dirichlet_bc(matrix, rhs)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp) :: x, y
        integer :: i, n
        
        n = size(matrix, 1)
        
        do i = 1, n
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            ! Check if on boundary
            if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
                abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
                ! Apply homogeneous Dirichlet BC
                matrix(i, :) = 0.0_dp
                matrix(:, i) = 0.0_dp
                matrix(i, i) = 1.0_dp
                rhs(i) = 0.0_dp
            end if
        end do
    end subroutine apply_homogeneous_dirichlet_bc
    
    subroutine solve_linear_system(matrix, rhs, solution)
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
    end subroutine solve_linear_system
    
    subroutine visualize_solution()
        real(dp), allocatable :: x_line(:), y_line(:)
        integer :: i, n_points
        
        ! Create visualization along y=0.5
        n_points = 100
        allocate(x_line(n_points))
        allocate(y_line(n_points))
        
        do i = 1, n_points
            x_line(i) = real(i-1, dp) / real(n_points-1, dp)
            y_line(i) = exact_solution(x_line(i), 0.5_dp)
        end do
        
        call figure()
        call plot(x_line, y_line, 'b-')
        call xlabel('x')
        call ylabel('u(x, 0.5)')
        call title('2D Poisson Solution Cross-Section at y=0.5')
        call savefig('poisson_2d_weak_form.png')
        
        print *, "Visualization saved to poisson_2d_weak_form.png"
        
        deallocate(x_line, y_line)
    end subroutine visualize_solution

end program poisson_2d_weak_form
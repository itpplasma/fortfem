program elasticity_demo
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
    type(trial_function_t) :: u1, u2  ! displacement components
    type(test_function_t) :: v1, v2   ! test functions
    
    ! Weak form components
    type(bilinear_form_t) :: a11, a12, a21, a22, a_total
    type(linear_form_t) :: L1, L2
    type(weak_form_t) :: problem
    
    ! Solution and material properties
    real(dp), allocatable :: solution(:), matrix(:,:), rhs(:)
    real(dp), allocatable :: disp_x(:), disp_y(:)
    real(dp) :: E, nu, lambda, mu  ! Material parameters
    integer :: nx, ny, n_dofs
    
    print *, "2D Linear Elasticity Demo"
    print *, "========================="
    print *, ""
    print *, "Demonstrating weak form framework for vector-valued problems"
    print *, "Problem: Linear elasticity with body forces"
    print *, ""
    
    ! Material properties (steel-like)
    E = 200.0e9_dp   ! Young's modulus (Pa)
    nu = 0.3_dp      ! Poisson's ratio
    lambda = E * nu / ((1.0_dp + nu) * (1.0_dp - 2.0_dp * nu))
    mu = E / (2.0_dp * (1.0_dp + nu))
    
    print '(a,es12.5)', "Young's modulus E = ", E
    print '(a,f6.3)', "Poisson's ratio ν = ", nu
    print '(a,es12.5)', "Lamé parameter λ = ", lambda
    print '(a,es12.5)', "Shear modulus μ = ", mu
    print *, ""
    
    ! Create mesh
    nx = 10
    ny = 10
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Create function space (for each component)
    call create_P1_space(mesh, V)
    
    ! Create trial and test functions for each displacement component
    call u1%init(V, "u1")  ! x-displacement
    call u2%init(V, "u2")  ! y-displacement
    call v1%init(V, "v1")  ! x-test function
    call v2%init(V, "v2")  ! y-test function
    
    n_dofs = 2 * V%n_dofs  ! Total DOFs for vector problem
    print '(a,i0)', "Total DOFs (2 × vertices): ", n_dofs
    print *, ""
    
    ! Define weak formulation for linear elasticity
    call setup_elasticity_weak_form()
    
    ! Assemble system
    print *, "Assembling elasticity system..."
    allocate(matrix(n_dofs, n_dofs))
    allocate(rhs(n_dofs))
    
    call assemble_elasticity_system()
    
    ! Apply boundary conditions
    call apply_elasticity_bc()
    
    ! Solve
    print *, "Solving elasticity system..."
    call solve_elasticity_system()
    
    ! Extract displacement components
    call extract_displacements()
    
    ! Analyze results
    call analyze_results()
    
    ! Visualize
    call visualize_elasticity()
    
    ! Clean up
    call cleanup()
    
    print *, ""
    print *, "Elasticity demo complete!"
    
contains

    subroutine setup_elasticity_weak_form()
        print *, "Setting up elasticity weak formulation..."
        print *, "  Bilinear form: a(u,v) = ∫[μ(∇u:∇v) + λ(∇·u)(∇·v)]dx"
        print *, "  Linear form:   L(v) = ∫f·v dx"
        print *, ""
        
        ! For simplicity, we'll create placeholder forms
        ! In a full implementation, these would properly compute the elasticity tensor
        
        ! (1,1) block: μ ∂u1/∂x1 ∂v1/∂x1 + μ ∂u1/∂x2 ∂v1/∂x2 + λ ∂u1/∂x1 ∂v1/∂x1
        call a11%init(form_type=2, coefficient=mu, &
                      expression="(mu + lambda)*grad(u1).grad(v1)")
        
        ! (1,2) block: λ ∂u2/∂x2 ∂v1/∂x1 + μ ∂u2/∂x1 ∂v1/∂x2
        call a12%init(form_type=3, coefficient=lambda, &
                      expression="lambda*grad(u2).grad(v1)")
        
        ! (2,1) block: μ ∂u1/∂x2 ∂v2/∂x1 + λ ∂u1/∂x1 ∂v2/∂x2
        call a21%init(form_type=3, coefficient=mu, &
                      expression="mu*grad(u1).grad(v2)")
        
        ! (2,2) block: μ ∂u2/∂x1 ∂v2/∂x1 + μ ∂u2/∂x2 ∂v2/∂x2 + λ ∂u2/∂x2 ∂v2/∂x2
        call a22%init(form_type=2, coefficient=mu, &
                      expression="(mu + lambda)*grad(u2).grad(v2)")
        
        ! Load vectors (body forces)
        call L1%init(form_type=1, coefficient=1.0_dp, &
                     expression="f1*v1")  ! x-component
        call L2%init(form_type=1, coefficient=1.0_dp, &
                     expression="f2*v2")  ! y-component
        
        print *, "Weak form components initialized"
    end subroutine setup_elasticity_weak_form
    
    subroutine assemble_elasticity_system()
        ! Simplified assembly - in practice this would use proper elasticity assembly
        real(dp), allocatable :: K11(:,:), K12(:,:), K21(:,:), K22(:,:)
        real(dp), allocatable :: f1(:), f2(:)
        integer :: n, i, j
        
        n = V%n_dofs
        allocate(K11(n,n), K12(n,n), K21(n,n), K22(n,n))
        allocate(f1(n), f2(n))
        
        ! Initialize matrices (placeholder - would use proper assembly)
        K11 = 0.0_dp
        K12 = 0.0_dp
        K21 = 0.0_dp
        K22 = 0.0_dp
        f1 = 0.0_dp
        f2 = 0.0_dp
        
        ! Simple diagonal entries for demo
        do i = 1, n
            K11(i,i) = mu + lambda
            K22(i,i) = mu + lambda
            f1(i) = 0.0_dp     ! No body force in x
            f2(i) = -1.0_dp    ! Gravity in y
        end do
        
        ! Assemble block matrix
        matrix = 0.0_dp
        rhs = 0.0_dp
        
        ! [K11  K12] [u1]   [f1]
        ! [K21  K22] [u2] = [f2]
        
        matrix(1:n, 1:n) = K11
        matrix(1:n, n+1:2*n) = K12
        matrix(n+1:2*n, 1:n) = K21
        matrix(n+1:2*n, n+1:2*n) = K22
        
        rhs(1:n) = f1
        rhs(n+1:2*n) = f2
        
        deallocate(K11, K12, K21, K22, f1, f2)
    end subroutine assemble_elasticity_system
    
    subroutine apply_elasticity_bc()
        real(dp) :: x, y
        integer :: i, n
        
        print *, "Applying boundary conditions..."
        print *, "  Fixed support at x = 0"
        print *, "  Free surface elsewhere"
        
        n = V%n_dofs
        
        do i = 1, n
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            ! Fixed support at x = 0 (both u1 and u2 = 0)
            if (abs(x) < 1e-10) then
                ! u1 = 0
                matrix(i, :) = 0.0_dp
                matrix(:, i) = 0.0_dp
                matrix(i, i) = 1.0_dp
                rhs(i) = 0.0_dp
                
                ! u2 = 0
                matrix(n+i, :) = 0.0_dp
                matrix(:, n+i) = 0.0_dp
                matrix(n+i, n+i) = 1.0_dp
                rhs(n+i) = 0.0_dp
            end if
        end do
    end subroutine apply_elasticity_bc
    
    subroutine solve_elasticity_system()
        integer :: info
        integer, allocatable :: ipiv(:)
        
        allocate(solution(n_dofs))
        allocate(ipiv(n_dofs))
        
        solution = rhs
        
        ! Solve using LAPACK
        call dgesv(n_dofs, 1, matrix, n_dofs, ipiv, solution, n_dofs, info)
        
        if (info /= 0) then
            print *, "Error in elasticity solver, info = ", info
            stop
        end if
        
        deallocate(ipiv)
    end subroutine solve_elasticity_system
    
    subroutine extract_displacements()
        integer :: n
        
        n = V%n_dofs
        allocate(disp_x(n), disp_y(n))
        
        disp_x = solution(1:n)
        disp_y = solution(n+1:2*n)
        
        ! Store in trial functions
        call u1%assign_values(disp_x)
        call u2%assign_values(disp_y)
    end subroutine extract_displacements
    
    subroutine analyze_results()
        real(dp) :: max_disp_x, max_disp_y, max_disp_total
        
        max_disp_x = maxval(abs(disp_x))
        max_disp_y = maxval(abs(disp_y))
        max_disp_total = maxval(sqrt(disp_x**2 + disp_y**2))
        
        print *, ""
        print *, "Results analysis:"
        print '(a,es12.5)', "  Max x-displacement: ", max_disp_x
        print '(a,es12.5)', "  Max y-displacement: ", max_disp_y
        print '(a,es12.5)', "  Max total displacement: ", max_disp_total
        print '(a,es12.5)', "  Min y-displacement: ", minval(disp_y)
    end subroutine analyze_results
    
    subroutine visualize_elasticity()
        real(dp), allocatable :: x_line(:), y_line(:)
        integer :: i, n_points, idx
        real(dp) :: x_val
        
        ! Create displacement plot along top edge
        n_points = min(20, nx)
        allocate(x_line(n_points))
        allocate(y_line(n_points))
        
        do i = 1, n_points
            x_val = real(i-1, dp) / real(n_points-1, dp)
            x_line(i) = x_val
            
            ! Find approximate index for (x_val, 1.0)
            idx = (ny-1) * nx + i
            if (idx <= V%n_dofs) then
                y_line(i) = disp_y(idx)
            else
                y_line(i) = 0.0_dp
            end if
        end do
        
        call figure()
        call plot(x_line, y_line, 'r-')
        call xlabel('x')
        call ylabel('y-displacement')
        call title('Elasticity: y-displacement along top edge')
        call savefig('elasticity_displacement.png')
        
        print *, "Displacement plot saved to elasticity_displacement.png"
        
        ! Write data file
        call write_displacement_data()
        
        deallocate(x_line, y_line)
    end subroutine visualize_elasticity
    
    subroutine write_displacement_data()
        integer :: i
        
        open(unit=10, file='elasticity_solution.dat', status='replace')
        write(10, '(a)') '# x y disp_x disp_y magnitude'
        do i = 1, V%n_dofs
            write(10, '(5es16.8)') mesh%vertices(1,i), mesh%vertices(2,i), &
                                   disp_x(i), disp_y(i), &
                                   sqrt(disp_x(i)**2 + disp_y(i)**2)
        end do
        close(10)
        
        print *, "Displacement data written to elasticity_solution.dat"
    end subroutine write_displacement_data
    
    subroutine cleanup()
        call u1%destroy()
        call u2%destroy()
        call v1%destroy()
        call v2%destroy()
        call V%destroy()
        call mesh%destroy()
        call a11%destroy()
        call a12%destroy()
        call a21%destroy()
        call a22%destroy()
        call L1%destroy()
        call L2%destroy()
        if (allocated(solution)) deallocate(solution)
        if (allocated(matrix)) deallocate(matrix)
        if (allocated(rhs)) deallocate(rhs)
        if (allocated(disp_x)) deallocate(disp_x)
        if (allocated(disp_y)) deallocate(disp_y)
    end subroutine cleanup

end program elasticity_demo
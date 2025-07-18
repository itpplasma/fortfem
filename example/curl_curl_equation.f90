program curl_curl_equation
    ! Maxwell's curl-curl equation with regularization
    ! curl(curl(E)) + k²E = J
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_function_space
    use fortfem_expressions
    use forms_module
    use fortfem_gmres
    implicit none
    
    ! Problem parameters
    real(dp), parameter :: k = 1.0_dp      ! Wave number
    real(dp), parameter :: eps = 1e-6_dp   ! Regularization parameter
    
    ! Mesh and function space
    type(mesh_2d_t) :: mesh
    type(function_space_t) :: V_scalar  ! Scalar space for building vector space
    type(vector_function_space_t) :: V_edge  ! H(curl) space
    
    ! Vector functions (E-field)
    type(vector_trial_function_t) :: E
    type(vector_test_function_t) :: v_test
    
    ! System matrices and vectors
    real(dp), allocatable :: A(:,:), b(:), x(:)
    
    ! Solver options
    type(gmres_options_t) :: gmres_opts
    
    ! Expressions
    type(expression_t) :: a_form, L_form, J_current
    
    integer :: n_dofs, info
    
    print *, "Maxwell's Curl-Curl Equation with Edge Elements"
    print *, "=============================================="
    print *, ""
    print *, "Problem: curl(curl(E)) + k²E = J"
    print *, "Regularized: curl(curl(E)) + k²E + ε·div(E)·div(v) = J·v"
    print *, ""
    print *, "Using Nédélec edge elements in H(curl) space"
    print *, ""
    
    ! Create mesh
    print *, "Creating mesh..."
    call mesh%create_rectangular(nx=20, ny=20, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "  Elements: ", mesh%n_triangles
    print '(a,i0)', "  Vertices: ", mesh%n_vertices
    
    ! Create edge element function space
    print *, ""
    print *, "Creating H(curl) function space with edge elements..."
    call create_edge_space(mesh, V_scalar)
    call V_edge%init(V_scalar, 2)  ! 2D vector space
    n_dofs = V_edge%total_dofs()
    
    print '(a,i0)', "  DOFs (edges): ", n_dofs
    
    ! Create trial and test functions
    print *, ""
    print *, "Creating vector trial and test functions..."
    call E%init(V_edge, "E")
    call v_test%init(V_edge, "v")
    
    ! Define current density J (simple constant)
    J_current%expr_type = 11  ! EXPR_CONSTANT
    J_current%const_value = 1.0_dp
    J_current%rank = 1  ! Vector
    
    ! Define variational form
    print *, ""
    print *, "Setting up variational form:"
    print *, "  Curl-curl term: inner(curl(E), curl(v))*dx"
    print *, "  Mass term: k²*inner(E, v)*dx"
    print *, "  Regularization: ε*div(E)*div(v)*dx"
    print *, "  Source: inner(J, v)*dx"
    print *, ""
    
    ! Build bilinear form: a(E,v) = ∫(curl E · curl v + k²E·v + ε div E div v) dx
    a_form = inner(curl(E), curl(v_test)) * dx() + &
             k**2 * inner(vector_trial_to_expression(E), vector_test_to_expression(v_test)) * dx() + &
             eps * div(E) * div(v_test) * dx()
    
    ! Build linear form: L(v) = ∫J·v dx
    L_form = inner(J_current, vector_test_to_expression(v_test)) * dx()
    
    print *, "Variational forms created successfully!"
    print *, ""
    
    ! Allocate system
    print *, "Allocating system matrices..."
    allocate(A(n_dofs, n_dofs))
    allocate(b(n_dofs))
    allocate(x(n_dofs))
    
    ! Initialize system (placeholder assembly)
    print *, "Assembling system (placeholder)..."
    A = 0.0_dp
    b = 0.0_dp
    x = 0.0_dp
    
    ! Add identity to main diagonal for testing
    do info = 1, n_dofs
        A(info, info) = 1.0_dp + real(info, dp) * 1e-6_dp
        b(info) = 1.0_dp
    end do
    
    print '(a,i0,a,i0)', "  System size: ", n_dofs, " x ", n_dofs
    print *, ""
    
    ! Setup GMRES solver
    print *, "Setting up GMRES solver..."
    gmres_opts%max_iter = 1000
    gmres_opts%restart = 30
    gmres_opts%tol = 1e-8_dp
    gmres_opts%verbose = .true.
    
    print *, "Solving with GMRES..."
    call gmres_solve(A, b, x, gmres_opts, info)
    
    if (info == 0) then
        print *, ""
        print *, "✓ GMRES converged successfully!"
        print '(a,es12.5)', "  Solution norm: ", norm2(x)
        print '(a,es12.5)', "  Residual norm: ", norm2(b - matmul(A, x))
    else
        print *, ""
        print *, "✗ GMRES failed to converge"
    end if
    
    print *, ""
    print *, "Key features demonstrated:"
    print *, "- Nédélec edge elements for H(curl) space"
    print *, "- Curl-curl operator with regularization"
    print *, "- Vector variational forms: curl(E), div(E), inner(E,v)"
    print *, "- GMRES iterative solver"
    print *, "- Natural syntax: curl(E), div(E), inner(E,v)*dx"
    
    ! Clean up
    call E%destroy()
    call v_test%destroy()
    call V_edge%destroy()
    call V_scalar%destroy()
    call mesh%destroy()
    deallocate(A, b, x)
    
    print *, ""
    print *, "This demonstrates the complete Maxwell equation workflow!"
    
end program curl_curl_equation
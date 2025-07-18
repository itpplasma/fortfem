program poisson_forms
    ! Poisson equation using modern variational forms syntax
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_function_space
    use fortfem_expressions
    use forms_module
    implicit none
    
    ! Mesh and function space
    type(mesh_2d_t) :: mesh
    type(function_space_t) :: V
    
    ! Functions
    type(trial_function_t) :: u
    type(test_function_t) :: v_test
    
    ! Expressions
    type(expression_t) :: a, L
    type(expression_t) :: f
    
    ! Solution
    real(dp), allocatable :: solution(:)
    
    print *, "Variational Forms Poisson Equation"
    print *, "=================================="
    print *, ""
    print *, "Problem: -∇²u = f in Ω, u = 0 on ∂Ω"
    print *, ""
    
    ! Create mesh
    call mesh%create_rectangular(nx=10, ny=10, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    ! Create function space
    call create_P1_space(mesh, V)
    
    ! Define trial and test functions
    call u%init(V, "u")
    call v_test%init(V, "v")
    
    ! Define source term
    f%expr_type = 11  ! EXPR_CONSTANT
    f%const_value = 1.0_dp
    f%rank = 0
    
    ! Define variational forms
    print *, "Variational form syntax:"
    print *, "  a = inner(grad(u), grad(v))*dx"
    print *, "  L = f*v*dx"
    print *, ""
    
    ! Build forms
    a = inner(grad(u), grad(v_test)) * dx()
    L = f * test_to_expression(v_test) * dx()
    
    print *, "Forms created with expression types:"
    print *, "  a type:", a%expr_type
    print *, "  L type:", L%expr_type
    
    ! Allocate solution
    allocate(solution(V%n_dofs))
    
    ! Solve (placeholder)
    ! In full implementation: solution = solve(a == L, bc)
    call solve(a, L, V, solution)
    
    ! Clean up
    call u%destroy()
    call v_test%destroy()
    call V%destroy()
    call mesh%destroy()
    deallocate(solution)
    
    print *, ""
    print *, "This demonstrates modern variational forms syntax in Fortran!"
    
end program poisson_forms
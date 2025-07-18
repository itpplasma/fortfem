program forms_examples
    ! Comprehensive variational forms examples
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
    
    ! Physical parameters
    real(dp), parameter :: nu = 0.01_dp
    real(dp), parameter :: k = 2.0_dp
    real(dp), parameter :: dt = 0.1_dp
    
    print *, "Variational Forms Examples"
    print *, "=========================="
    print *, ""
    
    ! Create mesh and function space
    call mesh%create_rectangular(nx=10, ny=10, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    call create_P1_space(mesh, V)
    call u%init(V, "u")
    call v_test%init(V, "v")
    
    ! Example 1: Poisson equation
    print *, "1. Poisson equation: -∇²u = f"
    print *, "   Modern syntax:"
    print *, "   a = inner(grad(u), grad(v))*dx"
    print *, "   L = f*v*dx"
    print *, ""
    
    ! Example 2: Helmholtz equation
    print *, "2. Helmholtz equation: -∇²u - k²u = f"
    print *, "   Modern syntax:"
    print *, "   a = (inner(grad(u), grad(v)) - k**2*u*v)*dx"
    print *, "   L = f*v*dx"
    print *, ""
    
    ! Example 3: Heat equation
    print *, "3. Heat equation: ∂u/∂t - ∇²u = f"
    print *, "   FEniCS syntax (backward Euler):"
    print *, "   a = (u*v + dt*inner(grad(u), grad(v)))*dx"
    print *, "   L = (u_n + dt*f)*v*dx"
    print *, ""
    
    ! Example 4: Convection-diffusion
    print *, "4. Convection-diffusion: ∂u/∂t + b·∇u - ν∇²u = f"
    print *, "   Modern syntax:"
    print *, "   a = (u*v + dt*inner(b, grad(u))*v + dt*nu*inner(grad(u), grad(v)))*dx"
    print *, "   L = u_n*v*dx + dt*f*v*dx"
    print *, ""
    
    ! Example 5: Neumann boundary conditions
    print *, "5. Poisson with Neumann BC: -∇²u = f, ∂u/∂n = g"
    print *, "   Modern syntax:"
    print *, "   a = inner(grad(u), grad(v))*dx"
    print *, "   L = f*v*dx + g*v*ds"
    print *, ""
    
    ! Example 6: Robin boundary conditions
    print *, "6. Robin BC: -∇²u = f, ∂u/∂n + αu = g"
    print *, "   Modern syntax:"
    print *, "   a = inner(grad(u), grad(v))*dx + alpha*u*v*ds"
    print *, "   L = f*v*dx + g*v*ds"
    print *, ""
    
    ! Example 7: Mixed boundary conditions
    print *, "7. Mixed BC on different boundaries:"
    print *, "   Modern syntax:"
    print *, "   a = inner(grad(u), grad(v))*dx + alpha*u*v*ds(1)"
    print *, "   L = f*v*dx + g*v*ds(1) + h*v*ds(2)"
    print *, ""
    
    ! Example 8: Subdomain integration
    print *, "8. Different materials (subdomains):"
    print *, "   Modern syntax:"
    print *, "   a = kappa1*inner(grad(u), grad(v))*dx(1) + kappa2*inner(grad(u), grad(v))*dx(2)"
    print *, "   L = f1*v*dx(1) + f2*v*dx(2)"
    print *, ""
    
    ! Example 9: Nonlinear problem (conceptual)
    print *, "9. Nonlinear diffusion (conceptual):"
    print *, "   FEniCS-style would need nonlinear solve:"
    print *, "   F = inner(pow(inner(grad(u), grad(u)), (p-2)/2)*grad(u), grad(v))*dx - f*v*dx"
    print *, "   solve(F == 0, u)"
    print *, ""
    
    ! Example 10: Vector problems (elasticity)
    print *, "10. Linear elasticity (vector-valued, conceptual):"
    print *, "    Modern syntax:"
    print *, "    a = inner(sigma(u), epsilon(v))*dx"
    print *, "    L = inner(f, v)*dx"
    print *, "    where sigma(u) = 2*mu*epsilon(u) + lmbda*tr(epsilon(u))*I"
    print *, ""
    
    ! Clean up
    call u%destroy()
    call v_test%destroy()
    call V%destroy()
    call mesh%destroy()
    
    print *, "These examples show the power and flexibility of modern variational forms syntax!"
    print *, "Key features:"
    print *, "- Natural mathematical notation"
    print *, "- Flexible integration measures dx, ds"
    print *, "- Support for subdomains and boundaries"
    print *, "- Type-safe operator overloading"
    print *, "- Extensible to new PDE types"
    
end program forms_examples
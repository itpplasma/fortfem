program weak_form_demo
    use fortfem
    use function_space_module
    use weak_forms_module
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(function_space_t) :: V
    type(trial_function_t) :: u
    type(test_function_t) :: v_test
    type(bilinear_form_t) :: a, mass_form, stiffness_form
    type(linear_form_t) :: L
    type(weak_form_t) :: problem
    
    print *, "Weak Form Framework Demo"
    print *, "========================"
    print *, ""
    
    ! Create mesh and function space
    call mesh%create_rectangular(nx=10, ny=10, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    ! Create P1 function space
    V = P1(mesh)
    
    print '(a,i0)', "Function space DOFs: ", V%n_dofs
    
    ! Create trial and test functions
    call u%init(V, "u")
    call v_test%init(V, "v")
    
    print *, "Created trial function: ", trim(u%name)
    print *, "Created test function: ", trim(v_test%name)
    
    ! Demonstrate operator overloading
    print *, ""
    print *, "Creating bilinear forms..."
    
    ! Mass matrix: ∫ u*v dx
    call mass_form%init(1, expression="u*v")
    
    ! Stiffness matrix: ∫ ∇u·∇v dx
    call stiffness_form%init(2, expression="grad(u).grad(v)")
    
    ! Combined form: mass + stiffness
    a = mass_form + stiffness_form
    
    print *, "Mass form: ", trim(mass_form%expression)
    print *, "Stiffness form: ", trim(stiffness_form%expression)
    print *, "Combined form: ", trim(a%expression)
    
    ! Create linear form
    print *, ""
    print *, "Creating linear form..."
    
    call L%init(1, expression="f*v")
    
    print *, "Linear form: ", trim(L%expression)
    
    ! Create weak form
    problem = a - L
    
    print *, ""
    print *, "Weak form: ", trim(problem%description)
    
    ! Demonstrate assembly (placeholder)
    print *, ""
    print *, "Assembling weak form..."
    
    call test_assembly()
    
    ! Clean up
    call u%destroy()
    call v_test%destroy()
    call V%destroy()
    call mesh%destroy()
    call mass_form%destroy()
    call stiffness_form%destroy()
    call a%destroy()
    call L%destroy()
    call problem%destroy()
    
    print *, ""
    print *, "Demo complete!"
    
contains

    subroutine test_assembly()
        real(dp), allocatable :: matrix(:,:), vector(:)
        
        allocate(matrix(V%n_dofs, V%n_dofs))
        allocate(vector(V%n_dofs))
        
        call problem%assemble(V, matrix, vector)
        
        print '(a,i0,a,i0)', "Assembled matrix: ", size(matrix,1), " x ", size(matrix,2)
        print '(a,i0)', "Assembled vector: ", size(vector)
        
        deallocate(matrix, vector)
    end subroutine test_assembly

end program weak_form_demo
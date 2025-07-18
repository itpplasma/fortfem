module test_weak_forms
    use fortfem
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private
    
    public :: collect_weak_forms
    
contains

    subroutine collect_weak_forms(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        
        testsuite = [ &
            new_unittest("function_space_creation", test_function_space_creation), &
            new_unittest("trial_function_init", test_trial_function_init), &
            new_unittest("test_function_init", test_test_function_init), &
            new_unittest("bilinear_form_creation", test_bilinear_form_creation), &
            new_unittest("linear_form_creation", test_linear_form_creation), &
            new_unittest("form_addition", test_form_addition), &
            new_unittest("weak_form_creation", test_weak_form_creation), &
            new_unittest("deep_copy_assignment", test_deep_copy_assignment) &
        ]
        
    end subroutine collect_weak_forms
    
    subroutine test_function_space_creation(error)
        type(error_type), allocatable, intent(out) :: error
        type(mesh_2d_t) :: mesh
        type(function_space_t) :: space
        
        ! Create mesh
        call mesh%create_rectangular(nx=3, ny=3, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        ! Create P1 function space
        call create_P1_space(mesh, space)
        
        ! Test properties
        call check(error, space%n_dofs == mesh%n_vertices, &
                   "P1 space should have n_dofs == n_vertices")
        call check(error, space%element_type == 1, &
                   "P1 space should have element_type == 1")
        call check(error, space%n_components == 1, &
                   "P1 space should be scalar (n_components == 1)")
        
        ! Clean up
        call space%destroy()
        call mesh%destroy()
        
    end subroutine test_function_space_creation
    
    subroutine test_trial_function_init(error)
        type(error_type), allocatable, intent(out) :: error
        type(mesh_2d_t) :: mesh
        type(function_space_t) :: space
        type(trial_function_t) :: u
        
        ! Create mesh and space
        call mesh%create_rectangular(nx=3, ny=3, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        call create_P1_space(mesh, space)
        
        ! Create trial function
        call u%init(space, "u")
        
        ! Test properties
        call check(error, u%name == "u", "Trial function name should be set")
        call check(error, allocated(u%coefficients), &
                   "Trial function coefficients should be allocated")
        call check(error, size(u%coefficients) == space%n_dofs, &
                   "Trial function coefficients size should match space DOFs")
        
        ! Clean up
        call u%destroy()
        call space%destroy()
        call mesh%destroy()
        
    end subroutine test_trial_function_init
    
    subroutine test_test_function_init(error)
        type(error_type), allocatable, intent(out) :: error
        type(mesh_2d_t) :: mesh
        type(function_space_t) :: space
        type(test_function_t) :: v
        
        ! Create mesh and space
        call mesh%create_rectangular(nx=3, ny=3, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        call create_P1_space(mesh, space)
        
        ! Create test function
        call v%init(space, "v")
        
        ! Test properties
        call check(error, v%name == "v", "Test function name should be set")
        
        ! Clean up
        call v%destroy()
        call space%destroy()
        call mesh%destroy()
        
    end subroutine test_test_function_init
    
    subroutine test_bilinear_form_creation(error)
        type(error_type), allocatable, intent(out) :: error
        type(bilinear_form_t) :: form
        
        ! Create bilinear form
        call form%init(form_type=1, coefficient=2.0_dp, expression="u*v")
        
        ! Test properties
        call check(error, form%form_type == 1, "Form type should be set")
        call check(error, abs(form%coefficient - 2.0_dp) < 1e-14, &
                   "Coefficient should be set correctly")
        call check(error, form%expression == "u*v", &
                   "Expression should be set correctly")
        
        ! Clean up
        call form%destroy()
        
    end subroutine test_bilinear_form_creation
    
    subroutine test_linear_form_creation(error)
        type(error_type), allocatable, intent(out) :: error
        type(linear_form_t) :: form
        
        ! Create linear form
        call form%init(form_type=1, coefficient=3.0_dp, expression="f*v")
        
        ! Test properties
        call check(error, form%form_type == 1, "Form type should be set")
        call check(error, abs(form%coefficient - 3.0_dp) < 1e-14, &
                   "Coefficient should be set correctly")
        call check(error, form%expression == "f*v", &
                   "Expression should be set correctly")
        
        ! Clean up
        call form%destroy()
        
    end subroutine test_linear_form_creation
    
    subroutine test_form_addition(error)
        type(error_type), allocatable, intent(out) :: error
        type(bilinear_form_t) :: form1, form2, result
        
        ! Create two forms
        call form1%init(form_type=1, expression="u*v")
        call form2%init(form_type=2, expression="grad(u).grad(v)")
        
        ! Add them
        call add_bilinear_forms(form1, form2, result)
        
        ! Test result
        call check(error, result%operation == 1, &
                   "Result should have addition operation")
        call check(error, allocated(result%left_form), &
                   "Left form should be allocated")
        call check(error, allocated(result%right_form), &
                   "Right form should be allocated")
        
        ! Clean up
        call form1%destroy()
        call form2%destroy()
        call result%destroy()
        
    end subroutine test_form_addition
    
    subroutine test_weak_form_creation(error)
        type(error_type), allocatable, intent(out) :: error
        type(bilinear_form_t) :: a
        type(linear_form_t) :: L
        type(weak_form_t) :: problem
        
        ! Create forms
        call a%init(form_type=2, expression="grad(u).grad(v)")
        call L%init(form_type=1, expression="f*v")
        
        ! Create weak form
        call problem%init(a, L, "Poisson problem")
        
        ! Test properties
        call check(error, problem%description == "Poisson problem", &
                   "Description should be set")
        
        ! Clean up
        call problem%destroy()
        call a%destroy()
        call L%destroy()
        
    end subroutine test_weak_form_creation
    
    subroutine test_deep_copy_assignment(error)
        type(error_type), allocatable, intent(out) :: error
        type(bilinear_form_t) :: original, copy
        
        ! Create original with nested structure
        call original%init(form_type=1, expression="original")
        allocate(original%left_form)
        call original%left_form%init(form_type=2, expression="left")
        
        ! Test deep copy assignment
        copy = original
        
        ! Modify original
        original%expression = "modified"
        original%left_form%expression = "modified_left"
        
        ! Test that copy is independent
        call check(error, copy%expression == "original", &
                   "Copy should retain original expression")
        call check(error, copy%left_form%expression == "left", &
                   "Copy should retain original left form expression")
        
        ! Clean up
        call original%destroy()
        call copy%destroy()
        
    end subroutine test_deep_copy_assignment

end module test_weak_forms
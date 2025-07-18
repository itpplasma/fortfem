module natural_forms_module
    use fortfem_kinds
    use function_space_module
    use expressions_module
    implicit none
    private
    
    public :: natural_form_t, form_t
    public :: grad, div, curl, inner, dot  ! Re-export from expressions
    public :: operator(*), operator(+), operator(-), operator(/)  ! Re-export operators
    public :: operator(==)  ! For boundary conditions
    public :: dx, ds  ! Re-export integration measures
    
    ! Form type - represents a linear/bilinear form (FEniCS style)
    type :: form_t
        type(expression_t) :: expr  ! The complete expression including measure
    contains
        procedure :: destroy => form_destroy
        generic :: assignment(=) => form_assign
        procedure :: form_assign
    end type form_t
    
    ! Problem type: represents a FEM problem  
    type :: problem_t
        type(form_t), allocatable :: forms(:)
        integer :: n_forms = 0
    contains
        procedure :: add_form
        procedure :: destroy => problem_destroy
        procedure :: assemble => problem_assemble
        generic :: assignment(=) => problem_assign
        procedure :: problem_assign
    end type problem_t
    
    ! Boundary condition type
    type :: bc_t
        type(expression_t) :: lhs, rhs
        integer :: bc_type = 1  ! 1=Dirichlet, 2=Neumann, 3=Robin
        integer :: boundary_id = 0
    end type bc_t
    
    ! Interface for creating integrals
    interface integral
        module procedure create_volume_integral
        module procedure create_surface_integral
    end interface
    
    ! Interface for boundary conditions
    interface operator(==)
        module procedure expression_equals_expression
        module procedure expression_equals_scalar
        module procedure trial_equals_scalar
    end interface
    
    ! Re-export grad interface with proper types
    interface grad
        module procedure grad_trial
        module procedure grad_test
    end interface
    
contains

    ! Create volume integral
    function create_volume_integral(expr) result(integ)
        type(expression_t), intent(in) :: expr
        type(integral_t) :: integ
        
        integ%integrand = expr
        integ%domain_type = 1  ! Volume
        integ%domain_id = 0    ! Whole domain
    end function create_volume_integral
    
    ! Create surface integral
    function create_surface_integral(expr, boundary_id) result(integ)
        type(expression_t), intent(in) :: expr
        integer, intent(in) :: boundary_id
        type(integral_t) :: integ
        
        integ%integrand = expr
        integ%domain_type = 2  ! Surface
        integ%domain_id = boundary_id
    end function create_surface_integral
    
    ! Gradient of trial function
    function grad_trial(u) result(grad_expr)
        type(trial_function_t), intent(in) :: u
        type(expression_t) :: grad_expr
        type(expression_t) :: u_expr
        
        u_expr = trial_to_expression(u)
        grad_expr = gradient_of_expression(u_expr)
    end function grad_trial
    
    ! Gradient of test function
    function grad_test(v) result(grad_expr)
        type(test_function_t), intent(in) :: v
        type(expression_t) :: grad_expr
        type(expression_t) :: v_expr
        
        v_expr = test_to_expression(v)
        grad_expr = gradient_of_expression(v_expr)
    end function grad_test
    
    ! Boundary condition operators
    function expression_equals_expression(lhs, rhs) result(bc)
        type(expression_t), intent(in) :: lhs, rhs
        type(bc_t) :: bc
        
        bc%lhs = lhs
        bc%rhs = rhs
        bc%bc_type = 1  ! Dirichlet by default
    end function expression_equals_expression
    
    function expression_equals_scalar(expr, value) result(bc)
        type(expression_t), intent(in) :: expr
        real(dp), intent(in) :: value
        type(bc_t) :: bc
        type(expression_t) :: const_expr
        
        bc%lhs = expr
        
        ! Create constant expression
        const_expr%expr_type = 11  ! EXPR_CONSTANT
        const_expr%rank = 0
        const_expr%const_value = value
        bc%rhs = const_expr
        
        bc%bc_type = 1  ! Dirichlet
    end function expression_equals_scalar
    
    function trial_equals_scalar(u, value) result(bc)
        type(trial_function_t), intent(in) :: u
        real(dp), intent(in) :: value
        type(bc_t) :: bc
        type(expression_t) :: u_expr, const_expr
        
        u_expr = trial_to_expression(u)
        bc%lhs = u_expr
        
        const_expr%expr_type = 11  ! EXPR_CONSTANT
        const_expr%rank = 0
        const_expr%const_value = value
        bc%rhs = const_expr
        
        bc%bc_type = 1  ! Dirichlet
    end function trial_equals_scalar
    
    ! Add integral to weak form
    subroutine add_integral(this, integ)
        class(natural_form_t), intent(inout) :: this
        type(integral_t), intent(in) :: integ
        type(integral_t), allocatable :: temp(:)
        
        if (allocated(this%integrals)) then
            allocate(temp(this%n_integrals + 1))
            temp(1:this%n_integrals) = this%integrals
            temp(this%n_integrals + 1) = integ
            
            deallocate(this%integrals)
            this%integrals = temp
        else
            allocate(this%integrals(1))
            this%integrals(1) = integ
        end if
        
        this%n_integrals = this%n_integrals + 1
    end subroutine add_integral
    
    ! Assemble weak form into matrix and vector
    subroutine natural_form_assemble(this, space, matrix, rhs)
        class(natural_form_t), intent(in) :: this
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:)
        real(dp), intent(out) :: rhs(:)
        integer :: i
        
        ! Initialize
        matrix = 0.0_dp
        rhs = 0.0_dp
        
        ! Assemble each integral
        do i = 1, this%n_integrals
            call assemble_integral(this%integrals(i), space, matrix, rhs)
        end do
    end subroutine natural_form_assemble
    
    ! Assemble a single integral (placeholder)
    subroutine assemble_integral(integ, space, matrix, rhs)
        type(integral_t), intent(in) :: integ
        type(function_space_t), intent(in) :: space
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        
        ! This would perform the actual finite element assembly
        ! based on the expression tree in integ%integrand
        ! For now, it's a placeholder
    end subroutine assemble_integral
    
    ! Destructor for integral
    subroutine integral_destroy(this)
        class(integral_t), intent(inout) :: this
        
        call this%integrand%destroy()
        this%domain_type = 0
        this%domain_id = 0
    end subroutine integral_destroy
    
    ! Destructor for weak form
    subroutine natural_form_destroy(this)
        class(natural_form_t), intent(inout) :: this
        integer :: i
        
        if (allocated(this%integrals)) then
            do i = 1, this%n_integrals
                call this%integrals(i)%destroy()
            end do
            deallocate(this%integrals)
        end if
        
        this%n_integrals = 0
    end subroutine natural_form_destroy
    
    ! Deep copy for integral
    subroutine integral_assign(this, other)
        class(integral_t), intent(inout) :: this
        type(integral_t), intent(in) :: other
        
        call this%destroy()
        
        this%integrand = other%integrand
        this%domain_type = other%domain_type
        this%domain_id = other%domain_id
    end subroutine integral_assign
    
    ! Deep copy for weak form
    subroutine natural_form_assign(this, other)
        class(natural_form_t), intent(inout) :: this
        type(natural_form_t), intent(in) :: other
        integer :: i
        
        call this%destroy()
        
        if (allocated(other%integrals)) then
            allocate(this%integrals(other%n_integrals))
            do i = 1, other%n_integrals
                this%integrals(i) = other%integrals(i)
            end do
        end if
        
        this%n_integrals = other%n_integrals
    end subroutine natural_form_assign

end module natural_forms_module
module expressions_module
    use fortfem_kinds
    use function_space_module
    implicit none
    private
    
    public :: expression_t, grad, div, curl, inner, dot
    public :: operator(*), operator(+), operator(-), operator(/)
    public :: dx, ds  ! Integration measures (dx=domain, ds=boundary)
    public :: trial_to_expression, test_to_expression, gradient_of_expression
    public :: vector_trial_to_expression, vector_test_to_expression
    
    ! Expression types
    integer, parameter :: EXPR_FUNCTION = 1
    integer, parameter :: EXPR_GRADIENT = 2
    integer, parameter :: EXPR_DIVERGENCE = 3
    integer, parameter :: EXPR_CURL = 4
    integer, parameter :: EXPR_PRODUCT = 5
    integer, parameter :: EXPR_SUM = 6
    integer, parameter :: EXPR_DIFFERENCE = 7
    integer, parameter :: EXPR_SCALAR_MULT = 8
    integer, parameter :: EXPR_DOT_PRODUCT = 9
    integer, parameter :: EXPR_CROSS_PRODUCT = 10
    integer, parameter :: EXPR_CONSTANT = 11
    integer, parameter :: EXPR_COORDINATE = 12
    integer, parameter :: EXPR_INNER_PRODUCT = 13
    integer, parameter :: EXPR_MEASURE = 14
    
    ! Base expression type - represents any mathematical expression
    type :: expression_t
        integer :: expr_type = 0
        integer :: rank = 0  ! 0=scalar, 1=vector, 2=tensor
        
        ! For function expressions
        type(trial_function_t), allocatable :: trial_func
        type(test_function_t), allocatable :: test_func
        type(vector_trial_function_t), allocatable :: vector_trial_func
        type(vector_test_function_t), allocatable :: vector_test_func
        logical :: is_trial = .false.
        logical :: is_test = .false.
        logical :: is_vector = .false.
        
        ! For constant values
        real(dp) :: const_value = 0.0_dp
        
        ! For binary operations
        type(expression_t), allocatable :: left, right
        
        ! For coordinate expressions
        integer :: coord_index = 0  ! 1=x, 2=y, 3=z
        
        ! For measure expressions (dx, ds)
        integer :: measure_type = 0  ! 1=dx (domain), 2=ds (boundary)
        integer :: subdomain_id = -1  ! -1=all, >=0 specific subdomain
        
    contains
        procedure :: destroy => expression_destroy
        procedure :: evaluate => expression_evaluate
        procedure :: is_constant => expression_is_constant
        generic :: assignment(=) => expression_assign
        procedure :: expression_assign
    end type expression_t
    
    ! Operator interfaces
    interface operator(*)
        module procedure multiply_expressions
        module procedure scalar_times_expression
        module procedure expression_times_scalar
    end interface
    
    interface operator(+)
        module procedure add_expressions
    end interface
    
    interface operator(-)
        module procedure subtract_expressions
        module procedure negate_expression
    end interface
    
    interface operator(/)
        module procedure divide_expression_by_scalar
    end interface
    
    interface grad
        module procedure gradient_of_expression
    end interface
    
    interface div
        module procedure divergence_of_expression
        module procedure div_vector_trial
        module procedure div_vector_test
    end interface
    
    interface curl
        module procedure curl_of_expression
        module procedure curl_vector_trial
        module procedure curl_vector_test
    end interface
    
    interface inner
        module procedure inner_product
    end interface
    
    interface dot
        module procedure inner_product  ! dot is alias for inner
    end interface
    
contains

    ! Create expression from trial function
    function trial_to_expression(u) result(expr)
        type(trial_function_t), intent(in) :: u
        type(expression_t) :: expr
        
        expr%expr_type = EXPR_FUNCTION
        expr%rank = 0  ! Scalar function
        allocate(expr%trial_func)
        expr%trial_func = u
        expr%is_trial = .true.
    end function trial_to_expression
    
    ! Create expression from test function
    function test_to_expression(v) result(expr)
        type(test_function_t), intent(in) :: v
        type(expression_t) :: expr
        
        expr%expr_type = EXPR_FUNCTION
        expr%rank = 0  ! Scalar function
        allocate(expr%test_func)
        expr%test_func = v
        expr%is_test = .true.
    end function test_to_expression
    
    ! Create expression from vector trial function
    function vector_trial_to_expression(u) result(expr)
        type(vector_trial_function_t), intent(in) :: u
        type(expression_t) :: expr
        
        expr%expr_type = EXPR_FUNCTION
        expr%rank = 1  ! Vector function
        allocate(expr%vector_trial_func)
        expr%vector_trial_func = u
        expr%is_trial = .true.
        expr%is_vector = .true.
    end function vector_trial_to_expression
    
    ! Create expression from vector test function
    function vector_test_to_expression(v) result(expr)
        type(vector_test_function_t), intent(in) :: v
        type(expression_t) :: expr
        
        expr%expr_type = EXPR_FUNCTION
        expr%rank = 1  ! Vector function
        allocate(expr%vector_test_func)
        expr%vector_test_func = v
        expr%is_test = .true.
        expr%is_vector = .true.
    end function vector_test_to_expression
    
    ! Gradient operator
    function gradient_of_expression(expr) result(grad_expr)
        type(expression_t), intent(in) :: expr
        type(expression_t) :: grad_expr
        
        grad_expr%expr_type = EXPR_GRADIENT
        grad_expr%rank = expr%rank + 1  ! Gradient increases rank by 1
        allocate(grad_expr%left)
        grad_expr%left = expr
    end function gradient_of_expression
    
    ! Multiplication of expressions
    function multiply_expressions(expr1, expr2) result(prod)
        type(expression_t), intent(in) :: expr1, expr2
        type(expression_t) :: prod
        
        ! Check for dot product of gradients
        if (expr1%expr_type == EXPR_GRADIENT .and. expr2%expr_type == EXPR_GRADIENT) then
            prod%expr_type = EXPR_DOT_PRODUCT
            prod%rank = 0  ! Dot product of vectors is scalar
        else
            prod%expr_type = EXPR_PRODUCT
            prod%rank = expr1%rank + expr2%rank
        end if
        
        allocate(prod%left, prod%right)
        prod%left = expr1
        prod%right = expr2
    end function multiply_expressions
    
    ! Scalar multiplication
    function scalar_times_expression(scalar, expr) result(prod)
        real(dp), intent(in) :: scalar
        type(expression_t), intent(in) :: expr
        type(expression_t) :: prod
        
        prod%expr_type = EXPR_SCALAR_MULT
        prod%rank = expr%rank
        prod%const_value = scalar
        allocate(prod%left)
        prod%left = expr
    end function scalar_times_expression
    
    function expression_times_scalar(expr, scalar) result(prod)
        type(expression_t), intent(in) :: expr
        real(dp), intent(in) :: scalar
        type(expression_t) :: prod
        
        prod = scalar_times_expression(scalar, expr)
    end function expression_times_scalar
    
    ! Addition of expressions
    function add_expressions(expr1, expr2) result(sum_expr)
        type(expression_t), intent(in) :: expr1, expr2
        type(expression_t) :: sum_expr
        
        if (expr1%rank /= expr2%rank) then
            error stop "Cannot add expressions of different ranks"
        end if
        
        sum_expr%expr_type = EXPR_SUM
        sum_expr%rank = expr1%rank
        allocate(sum_expr%left, sum_expr%right)
        sum_expr%left = expr1
        sum_expr%right = expr2
    end function add_expressions
    
    ! Subtraction of expressions
    function subtract_expressions(expr1, expr2) result(diff)
        type(expression_t), intent(in) :: expr1, expr2
        type(expression_t) :: diff
        
        if (expr1%rank /= expr2%rank) then
            error stop "Cannot subtract expressions of different ranks"
        end if
        
        diff%expr_type = EXPR_DIFFERENCE
        diff%rank = expr1%rank
        allocate(diff%left, diff%right)
        diff%left = expr1
        diff%right = expr2
    end function subtract_expressions
    
    ! Negation
    function negate_expression(expr) result(neg)
        type(expression_t), intent(in) :: expr
        type(expression_t) :: neg
        
        neg = scalar_times_expression(-1.0_dp, expr)
    end function negate_expression
    
    ! Division by scalar
    function divide_expression_by_scalar(expr, scalar) result(quot)
        type(expression_t), intent(in) :: expr
        real(dp), intent(in) :: scalar
        type(expression_t) :: quot
        
        if (abs(scalar) < 1e-14) then
            error stop "Division by zero in expression"
        end if
        
        quot = scalar_times_expression(1.0_dp/scalar, expr)
    end function divide_expression_by_scalar
    
    ! Divergence operator
    function divergence_of_expression(expr) result(div_expr)
        type(expression_t), intent(in) :: expr
        type(expression_t) :: div_expr
        
        if (expr%rank < 1) then
            error stop "Cannot take divergence of scalar"
        end if
        
        div_expr%expr_type = EXPR_DIVERGENCE
        div_expr%rank = expr%rank - 1  ! Divergence decreases rank by 1
        allocate(div_expr%left)
        div_expr%left = expr
    end function divergence_of_expression
    
    ! Curl operator
    function curl_of_expression(expr) result(curl_expr)
        type(expression_t), intent(in) :: expr
        type(expression_t) :: curl_expr
        
        if (expr%rank /= 1) then
            error stop "Curl only defined for vector fields"
        end if
        
        curl_expr%expr_type = EXPR_CURL
        curl_expr%rank = 1  ! Curl of vector is vector (in 3D)
        allocate(curl_expr%left)
        curl_expr%left = expr
    end function curl_of_expression
    
    ! Curl of vector trial function
    function curl_vector_trial(u) result(curl_expr)
        type(vector_trial_function_t), intent(in) :: u
        type(expression_t) :: curl_expr
        
        curl_expr%expr_type = EXPR_CURL
        curl_expr%rank = 0  ! Curl of 2D vector is scalar
        allocate(curl_expr%left)
        curl_expr%left = vector_trial_to_expression(u)
    end function curl_vector_trial
    
    ! Curl of vector test function
    function curl_vector_test(v) result(curl_expr)
        type(vector_test_function_t), intent(in) :: v
        type(expression_t) :: curl_expr
        
        curl_expr%expr_type = EXPR_CURL
        curl_expr%rank = 0  ! Curl of 2D vector is scalar
        allocate(curl_expr%left)
        curl_expr%left = vector_test_to_expression(v)
    end function curl_vector_test
    
    ! Div of vector trial function
    function div_vector_trial(u) result(div_expr)
        type(vector_trial_function_t), intent(in) :: u
        type(expression_t) :: div_expr
        
        div_expr%expr_type = EXPR_DIVERGENCE
        div_expr%rank = 0  ! Div of vector is scalar
        allocate(div_expr%left)
        div_expr%left = vector_trial_to_expression(u)
    end function div_vector_trial
    
    ! Div of vector test function
    function div_vector_test(v) result(div_expr)
        type(vector_test_function_t), intent(in) :: v
        type(expression_t) :: div_expr
        
        div_expr%expr_type = EXPR_DIVERGENCE
        div_expr%rank = 0  ! Div of vector is scalar
        allocate(div_expr%left)
        div_expr%left = vector_test_to_expression(v)
    end function div_vector_test
    
    ! Inner product (FEniCS style)
    function inner_product(expr1, expr2) result(prod)
        type(expression_t), intent(in) :: expr1, expr2
        type(expression_t) :: prod
        
        ! For scalars, inner is just multiplication
        if (expr1%rank == 0 .and. expr2%rank == 0) then
            prod = multiply_expressions(expr1, expr2)
            return
        end if
        
        ! For vectors, inner is dot product
        if (expr1%rank == 1 .and. expr2%rank == 1) then
            prod%expr_type = EXPR_INNER_PRODUCT
            prod%rank = 0  ! Inner product of vectors is scalar
            allocate(prod%left, prod%right)
            prod%left = expr1
            prod%right = expr2
            return
        end if
        
        ! For tensors, inner is contraction
        if (expr1%rank == 2 .and. expr2%rank == 2) then
            prod%expr_type = EXPR_INNER_PRODUCT
            prod%rank = 0  ! Full contraction gives scalar
            allocate(prod%left, prod%right)
            prod%left = expr1
            prod%right = expr2
            return
        end if
        
        error stop "Inner product not defined for these ranks"
    end function inner_product
    
    ! Integration measures (FEniCS style)
    function dx(subdomain_id) result(expr)
        integer, optional, intent(in) :: subdomain_id
        type(expression_t) :: expr
        
        expr%expr_type = EXPR_MEASURE
        expr%rank = 0
        expr%measure_type = 1  ! Domain integral
        
        if (present(subdomain_id)) then
            expr%subdomain_id = subdomain_id
        else
            expr%subdomain_id = -1  ! All domains
        end if
    end function dx
    
    function ds(boundary_id) result(expr)
        integer, optional, intent(in) :: boundary_id
        type(expression_t) :: expr
        
        expr%expr_type = EXPR_MEASURE
        expr%rank = 0
        expr%measure_type = 2  ! Boundary integral
        
        if (present(boundary_id)) then
            expr%subdomain_id = boundary_id
        else
            expr%subdomain_id = -1  ! All boundaries
        end if
    end function ds
    
    ! Evaluate expression at a point (placeholder for actual implementation)
    function expression_evaluate(this, x, y, z) result(value)
        class(expression_t), intent(in) :: this
        real(dp), intent(in), optional :: x, y, z
        real(dp) :: value
        
        ! This would evaluate the expression tree
        value = 0.0_dp
    end function expression_evaluate
    
    ! Check if expression is constant
    function expression_is_constant(this) result(is_const)
        class(expression_t), intent(in) :: this
        logical :: is_const
        
        is_const = (this%expr_type == EXPR_CONSTANT)
    end function expression_is_constant
    
    ! Destroy expression
    subroutine expression_destroy(this)
        class(expression_t), intent(inout) :: this
        
        if (allocated(this%trial_func)) deallocate(this%trial_func)
        if (allocated(this%test_func)) deallocate(this%test_func)
        if (allocated(this%vector_trial_func)) deallocate(this%vector_trial_func)
        if (allocated(this%vector_test_func)) deallocate(this%vector_test_func)
        if (allocated(this%left)) then
            call this%left%destroy()
            deallocate(this%left)
        end if
        if (allocated(this%right)) then
            call this%right%destroy()
            deallocate(this%right)
        end if
        
        this%expr_type = 0
        this%rank = 0
    end subroutine expression_destroy
    
    ! Deep copy assignment
    subroutine expression_assign(this, other)
        class(expression_t), intent(inout) :: this
        type(expression_t), intent(in) :: other
        
        call this%destroy()
        
        this%expr_type = other%expr_type
        this%rank = other%rank
        this%is_trial = other%is_trial
        this%is_test = other%is_test
        this%is_vector = other%is_vector
        this%const_value = other%const_value
        this%coord_index = other%coord_index
        this%measure_type = other%measure_type
        this%subdomain_id = other%subdomain_id
        
        if (allocated(other%trial_func)) then
            allocate(this%trial_func)
            this%trial_func = other%trial_func
        end if
        
        if (allocated(other%test_func)) then
            allocate(this%test_func)
            this%test_func = other%test_func
        end if
        
        if (allocated(other%vector_trial_func)) then
            allocate(this%vector_trial_func)
            this%vector_trial_func = other%vector_trial_func
        end if
        
        if (allocated(other%vector_test_func)) then
            allocate(this%vector_test_func)
            this%vector_test_func = other%vector_test_func
        end if
        
        if (allocated(other%left)) then
            allocate(this%left)
            this%left = other%left
        end if
        
        if (allocated(other%right)) then
            allocate(this%right)
            this%right = other%right
        end if
    end subroutine expression_assign

end module expressions_module
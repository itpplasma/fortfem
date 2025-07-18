module weak_forms_module
    use fortfem_kinds
    use function_space_module
    implicit none
    private
    
    public :: bilinear_form_t, linear_form_t, weak_form_t
    public :: operator(+), operator(-), operator(*)
    public :: grad
    public :: operator(.dot.)
    public :: constant, zero
    
    ! Forward declarations for operator overloading
    interface operator(+)
        module procedure add_bilinear_forms
        module procedure add_linear_forms
    end interface
    
    interface operator(-)
        module procedure subtract_bilinear_forms
        module procedure subtract_linear_forms
        module procedure bilinear_minus_linear  ! a(u,v) - L(v) = 0
    end interface
    
    interface operator(*)
        module procedure scalar_times_trial_function
        module procedure trial_function_times_scalar
        module procedure scalar_times_test_function
        module procedure test_function_times_scalar
    end interface
    
    interface operator(.dot.)
        module procedure dot_product_trial_test
    end interface
    
    interface grad
        module procedure grad_trial_function
        module procedure grad_test_function
    end interface
    
    ! Bilinear form: a(u,v) - represents integrals involving both trial and test functions
    type :: bilinear_form_t
        character(len=256) :: expression = ""
        integer :: form_type = 0  ! 1=mass, 2=stiffness, 3=custom
        real(dp) :: coefficient = 1.0_dp
        logical :: is_symmetric = .false.
        
        ! For composite forms
        type(bilinear_form_t), allocatable :: left_form, right_form
        integer :: operation = 0  ! 1=add, 2=subtract
    contains
        procedure :: init => bilinear_form_init
        procedure :: destroy => bilinear_form_destroy
        procedure :: assemble => assemble_bilinear_form
        procedure :: is_composite
        generic :: assignment(=) => bilinear_form_assign
        procedure :: bilinear_form_assign
    end type bilinear_form_t
    
    ! Linear form: L(v) - represents integrals involving only test functions
    type :: linear_form_t
        character(len=256) :: expression = ""
        integer :: form_type = 0  ! 1=load, 2=boundary, 3=custom
        real(dp) :: coefficient = 1.0_dp
        
        ! For composite forms
        type(linear_form_t), allocatable :: left_form, right_form
        integer :: operation = 0  ! 1=add, 2=subtract
        
        ! Source function (if applicable)
        procedure(source_function_interface), pointer, nopass :: source_func => null()
    contains
        procedure :: init => linear_form_init
        procedure :: destroy => linear_form_destroy
        procedure :: assemble => assemble_linear_form
        procedure :: is_composite => linear_form_is_composite
        generic :: assignment(=) => linear_form_assign
        procedure :: linear_form_assign
    end type linear_form_t
    
    ! Complete weak form: Find u such that a(u,v) = L(v) for all v
    type :: weak_form_t
        type(bilinear_form_t) :: bilinear_form
        type(linear_form_t) :: linear_form
        character(len=256) :: description = ""
    contains
        procedure :: init => weak_form_init
        procedure :: destroy => weak_form_destroy
        procedure :: assemble => assemble_weak_form
    end type weak_form_t
    
    ! Differential operator results
    type :: gradient_t
        type(trial_function_t), allocatable :: trial_func
        type(test_function_t), allocatable :: test_func
        logical :: is_trial = .false.
        logical :: is_test = .false.
    contains
        procedure :: destroy => gradient_destroy
        generic :: assignment(=) => gradient_assign
        procedure :: gradient_assign
    end type gradient_t
    
    ! Abstract interface for source functions
    abstract interface
        pure function source_function_interface(x, y) result(val)
            import :: dp
            real(dp), intent(in) :: x, y
            real(dp) :: val
        end function source_function_interface
    end interface
    
    ! Form type constants
    integer, parameter :: FORM_MASS = 1
    integer, parameter :: FORM_STIFFNESS = 2
    integer, parameter :: FORM_CUSTOM = 3
    integer, parameter :: FORM_LOAD = 1
    integer, parameter :: FORM_BOUNDARY = 2
    
    ! Operation constants
    integer, parameter :: OP_ADD = 1
    integer, parameter :: OP_SUBTRACT = 2
    
contains

    ! Bilinear form methods
    subroutine bilinear_form_init(this, form_type, coefficient, expression)
        class(bilinear_form_t), intent(inout) :: this
        integer, intent(in) :: form_type
        real(dp), intent(in), optional :: coefficient
        character(len=*), intent(in), optional :: expression
        
        this%form_type = form_type
        if (present(coefficient)) this%coefficient = coefficient
        if (present(expression)) this%expression = expression
        
        ! Set symmetry based on form type
        select case (form_type)
        case (FORM_MASS, FORM_STIFFNESS)
            this%is_symmetric = .true.
        case default
            this%is_symmetric = .false.
        end select
    end subroutine bilinear_form_init
    
    subroutine bilinear_form_destroy(this)
        class(bilinear_form_t), intent(inout) :: this
        
        if (allocated(this%left_form)) then
            call this%left_form%destroy()
            deallocate(this%left_form)
        end if
        
        if (allocated(this%right_form)) then
            call this%right_form%destroy()
            deallocate(this%right_form)
        end if
        
        this%form_type = 0
        this%coefficient = 0.0_dp
        this%operation = 0
    end subroutine bilinear_form_destroy
    
    function is_composite(this) result(composite)
        class(bilinear_form_t), intent(in) :: this
        logical :: composite
        
        composite = (this%operation > 0)
    end function is_composite
    
    subroutine assemble_bilinear_form(this, space, matrix)
        class(bilinear_form_t), intent(in) :: this
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:)
        
        ! This would assemble the bilinear form into a matrix
        ! For now, just a placeholder
        matrix = 0.0_dp
        
        if (this%is_composite()) then
            ! Handle composite forms recursively
            print *, "Assembling composite bilinear form: ", trim(this%expression)
        else
            ! Handle simple forms
            select case (this%form_type)
            case (FORM_MASS)
                print *, "Assembling mass matrix"
            case (FORM_STIFFNESS)
                print *, "Assembling stiffness matrix"
            case (FORM_CUSTOM)
                print *, "Assembling custom bilinear form: ", trim(this%expression)
            end select
        end if
    end subroutine assemble_bilinear_form
    
    ! Linear form methods
    subroutine linear_form_init(this, form_type, coefficient, expression)
        class(linear_form_t), intent(inout) :: this
        integer, intent(in) :: form_type
        real(dp), intent(in), optional :: coefficient
        character(len=*), intent(in), optional :: expression
        
        this%form_type = form_type
        if (present(coefficient)) this%coefficient = coefficient
        if (present(expression)) this%expression = expression
    end subroutine linear_form_init
    
    subroutine linear_form_destroy(this)
        class(linear_form_t), intent(inout) :: this
        
        if (allocated(this%left_form)) then
            call this%left_form%destroy()
            deallocate(this%left_form)
        end if
        
        if (allocated(this%right_form)) then
            call this%right_form%destroy()
            deallocate(this%right_form)
        end if
        
        this%form_type = 0
        this%coefficient = 0.0_dp
        this%operation = 0
        this%source_func => null()
    end subroutine linear_form_destroy
    
    function linear_form_is_composite(this) result(composite)
        class(linear_form_t), intent(in) :: this
        logical :: composite
        
        composite = (this%operation > 0)
    end function linear_form_is_composite
    
    subroutine assemble_linear_form(this, space, vector)
        class(linear_form_t), intent(in) :: this
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: vector(:)
        
        ! This would assemble the linear form into a vector
        vector = 0.0_dp
        
        if (this%is_composite()) then
            print *, "Assembling composite linear form: ", trim(this%expression)
        else
            select case (this%form_type)
            case (FORM_LOAD)
                print *, "Assembling load vector"
            case (FORM_BOUNDARY)
                print *, "Assembling boundary linear form"
            case (FORM_CUSTOM)
                print *, "Assembling custom linear form: ", trim(this%expression)
            end select
        end if
    end subroutine assemble_linear_form
    
    ! Weak form methods
    subroutine weak_form_init(this, bilinear_form, linear_form, description)
        class(weak_form_t), intent(inout) :: this
        type(bilinear_form_t), intent(in) :: bilinear_form
        type(linear_form_t), intent(in) :: linear_form
        character(len=*), intent(in), optional :: description
        
        this%bilinear_form = bilinear_form
        this%linear_form = linear_form
        if (present(description)) this%description = description
    end subroutine weak_form_init
    
    subroutine weak_form_destroy(this)
        class(weak_form_t), intent(inout) :: this
        
        call this%bilinear_form%destroy()
        call this%linear_form%destroy()
        this%description = ""
    end subroutine weak_form_destroy
    
    subroutine assemble_weak_form(this, space, matrix, vector)
        class(weak_form_t), intent(in) :: this
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:)
        real(dp), intent(out) :: vector(:)
        
        call this%bilinear_form%assemble(space, matrix)
        call this%linear_form%assemble(space, vector)
    end subroutine assemble_weak_form
    
    ! Operator overloading implementations
    function add_bilinear_forms(a, b) result(c)
        type(bilinear_form_t), intent(in) :: a, b
        type(bilinear_form_t) :: c
        
        call c%init(FORM_CUSTOM, expression="(" // trim(a%expression) // ") + (" // trim(b%expression) // ")")
        
        allocate(c%left_form, c%right_form)
        c%left_form = a
        c%right_form = b
        c%operation = OP_ADD
    end function add_bilinear_forms
    
    function subtract_bilinear_forms(a, b) result(c)
        type(bilinear_form_t), intent(in) :: a, b
        type(bilinear_form_t) :: c
        
        call c%init(FORM_CUSTOM, expression="(" // trim(a%expression) // ") - (" // trim(b%expression) // ")")
        
        allocate(c%left_form, c%right_form)
        c%left_form = a
        c%right_form = b
        c%operation = OP_SUBTRACT
    end function subtract_bilinear_forms
    
    function add_linear_forms(a, b) result(c)
        type(linear_form_t), intent(in) :: a, b
        type(linear_form_t) :: c
        
        call c%init(FORM_CUSTOM, expression="(" // trim(a%expression) // ") + (" // trim(b%expression) // ")")
        
        allocate(c%left_form, c%right_form)
        c%left_form = a
        c%right_form = b
        c%operation = OP_ADD
    end function add_linear_forms
    
    function subtract_linear_forms(a, b) result(c)
        type(linear_form_t), intent(in) :: a, b
        type(linear_form_t) :: c
        
        call c%init(FORM_CUSTOM, expression="(" // trim(a%expression) // ") - (" // trim(b%expression) // ")")
        
        allocate(c%left_form, c%right_form)
        c%left_form = a
        c%right_form = b
        c%operation = OP_SUBTRACT
    end function subtract_linear_forms
    
    function bilinear_minus_linear(a, L) result(weak_form)
        type(bilinear_form_t), intent(in) :: a
        type(linear_form_t), intent(in) :: L
        type(weak_form_t) :: weak_form
        
        call weak_form%init(a, L, "a(u,v) - L(v) = 0")
    end function bilinear_minus_linear
    
    ! Scalar multiplication
    function scalar_times_trial_function(alpha, u) result(form)
        real(dp), intent(in) :: alpha
        type(trial_function_t), intent(in) :: u
        type(bilinear_form_t) :: form
        
        call form%init(FORM_CUSTOM, coefficient=alpha, expression=trim(u%name))
    end function scalar_times_trial_function
    
    function trial_function_times_scalar(u, alpha) result(form)
        type(trial_function_t), intent(in) :: u
        real(dp), intent(in) :: alpha
        type(bilinear_form_t) :: form
        
        call form%init(FORM_CUSTOM, coefficient=alpha, expression=trim(u%name))
    end function trial_function_times_scalar
    
    function scalar_times_test_function(alpha, v) result(form)
        real(dp), intent(in) :: alpha
        type(test_function_t), intent(in) :: v
        type(linear_form_t) :: form
        
        call form%init(FORM_CUSTOM, coefficient=alpha, expression=trim(v%name))
    end function scalar_times_test_function
    
    function test_function_times_scalar(v, alpha) result(form)
        type(test_function_t), intent(in) :: v
        real(dp), intent(in) :: alpha
        type(linear_form_t) :: form
        
        call form%init(FORM_CUSTOM, coefficient=alpha, expression=trim(v%name))
    end function test_function_times_scalar
    
    ! Differential operators
    function grad_trial_function(u) result(grad_u)
        type(trial_function_t), intent(in) :: u
        type(gradient_t) :: grad_u
        
        allocate(grad_u%trial_func)
        grad_u%trial_func = u
        grad_u%is_trial = .true.
    end function grad_trial_function
    
    function grad_test_function(v) result(grad_v)
        type(test_function_t), intent(in) :: v
        type(gradient_t) :: grad_v
        
        allocate(grad_v%test_func)
        grad_v%test_func = v
        grad_v%is_test = .true.
    end function grad_test_function
    
    subroutine gradient_destroy(this)
        class(gradient_t), intent(inout) :: this
        
        if (allocated(this%trial_func)) deallocate(this%trial_func)
        if (allocated(this%test_func)) deallocate(this%test_func)
        this%is_trial = .false.
        this%is_test = .false.
    end subroutine gradient_destroy
    
    ! Dot product of gradients
    function dot_product_trial_test(grad_u, grad_v) result(form)
        type(gradient_t), intent(in) :: grad_u, grad_v
        type(bilinear_form_t) :: form
        
        if (grad_u%is_trial .and. grad_v%is_test) then
            call form%init(FORM_STIFFNESS, expression="grad(" // &
                           trim(grad_u%trial_func%name) // ") . grad(" // &
                           trim(grad_v%test_func%name) // ")")
        else
            error stop "Invalid gradient dot product: need trial and test function"
        end if
    end function dot_product_trial_test
    
    
    ! Utility functions
    function constant(value) result(f)
        real(dp), intent(in) :: value
        real(dp) :: f
        f = value
    end function constant
    
    function zero() result(f)
        real(dp) :: f
        f = 0.0_dp
    end function zero

    ! Deep-copy assignment operators
    subroutine bilinear_form_assign(this, other)
        class(bilinear_form_t), intent(inout) :: this
        type(bilinear_form_t), intent(in) :: other
        
        call this%destroy()
        
        this%expression = other%expression
        this%form_type = other%form_type
        this%coefficient = other%coefficient
        this%is_symmetric = other%is_symmetric
        this%operation = other%operation
        
        if (allocated(other%left_form)) then
            allocate(this%left_form)
            this%left_form = other%left_form
        end if
        
        if (allocated(other%right_form)) then
            allocate(this%right_form)
            this%right_form = other%right_form
        end if
    end subroutine bilinear_form_assign

    subroutine linear_form_assign(this, other)
        class(linear_form_t), intent(inout) :: this
        type(linear_form_t), intent(in) :: other
        
        call this%destroy()
        
        this%expression = other%expression
        this%form_type = other%form_type
        this%coefficient = other%coefficient
        this%operation = other%operation
        this%source_func => other%source_func
        
        if (allocated(other%left_form)) then
            allocate(this%left_form)
            this%left_form = other%left_form
        end if
        
        if (allocated(other%right_form)) then
            allocate(this%right_form)
            this%right_form = other%right_form
        end if
    end subroutine linear_form_assign

    subroutine gradient_assign(this, other)
        class(gradient_t), intent(inout) :: this
        type(gradient_t), intent(in) :: other
        
        call this%destroy()
        
        if (allocated(other%trial_func)) then
            allocate(this%trial_func)
            this%trial_func = other%trial_func
        end if
        
        if (allocated(other%test_func)) then
            allocate(this%test_func)
            this%test_func = other%test_func
        end if
        
        this%is_trial = other%is_trial
        this%is_test = other%is_test
    end subroutine gradient_assign

end module weak_forms_module
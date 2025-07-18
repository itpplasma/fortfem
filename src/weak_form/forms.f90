module forms_module
    use fortfem_kinds
    use function_space_module
    use expressions_module
    implicit none
    private
    
    public :: grad, div, curl, inner, dot
    public :: operator(*), operator(+), operator(-), operator(/)
    public :: dx, ds
    public :: solve
    public :: trial_to_expression, test_to_expression
    public :: vector_trial_to_expression, vector_test_to_expression
    
    ! Gradient interfaces
    interface grad
        module procedure grad_trial
        module procedure grad_test
    end interface
    
    ! Note: curl and div interfaces are already provided by expressions_module
    
contains

    ! Gradient of trial function
    function grad_trial(u) result(grad_expr)
        type(trial_function_t), intent(in) :: u
        type(expression_t) :: grad_expr
        
        grad_expr = grad(trial_to_expression(u))
    end function grad_trial
    
    ! Gradient of test function
    function grad_test(v) result(grad_expr)
        type(test_function_t), intent(in) :: v
        type(expression_t) :: grad_expr
        
        grad_expr = grad(test_to_expression(v))
    end function grad_test
    
    ! Solve a variational problem (placeholder)
    subroutine solve(a_expr, L_expr, space, solution, bc)
        type(expression_t), intent(in) :: a_expr  ! Bilinear form a(u,v)
        type(expression_t), intent(in) :: L_expr  ! Linear form L(v)
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: solution(:)
        real(dp), optional, intent(in) :: bc(:)
        
        ! This would assemble and solve the system
        print *, "Solving variational problem..."
        print *, "  a = ", a_expr%expr_type
        print *, "  L = ", L_expr%expr_type
        
        ! Placeholder
        solution = 0.0_dp
    end subroutine solve

end module forms_module
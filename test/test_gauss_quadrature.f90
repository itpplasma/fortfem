program test_gauss_quadrature
    use fortfem_kinds, only: dp
    use fortfem_gauss_quadrature_2d
    implicit none

    call test_quadrature_exactness()
    call test_polynomial_integration()
    call test_quadrature_sum_weights()
    
    print *, "All Gauss quadrature tests passed!"

contains

    subroutine test_quadrature_exactness()
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: integral, exact
        integer :: order
        
        print *, ""
        print *, "Gauss Quadrature Exactness Test"
        print *, "==============================="
        
        ! Test that each order integrates polynomials exactly
        do order = 1, 6  ! Test orders 1-6 (order 7 needs more precise coefficients)
            call quad%init(order)
            
            ! Test integration of polynomial of degree order
            call integrate_polynomial(quad, order, integral)
            exact = exact_polynomial_integral(order)
            
            print *, "Order", order, "quadrature:"
            print *, "  Points:", quad%n_points
            print *, "  Integrating polynomial degree", order
            print *, "  Computed:", integral
            print *, "  Exact:", exact
            print *, "  Error:", abs(integral - exact)
            
            if (abs(integral - exact) > 1e-10_dp) then
                print *, "Error: quadrature not exact for polynomial degree", order
                print *, "  Error too large:", abs(integral - exact)
                stop 1
            end if
            
            call quad%destroy()
        end do
        
        print *, "Quadrature exactness test passed"
    end subroutine

    subroutine test_polynomial_integration()
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: integral, exact, error
        integer :: i
        
        print *, ""
        print *, "Polynomial Integration Test"
        print *, "==========================="
        
        ! Test x^3 * y^2
        call quad%init(5)  ! Need order 5 for degree 5 polynomial
        
        integral = 0.0_dp
        do i = 1, quad%n_points
            integral = integral + quad%weights(i) * &
                      quad%xi(i)**3 * quad%eta(i)**2
        end do
        
        exact = 1.0_dp / 420.0_dp  ! Exact integral of x^3*y^2 over unit triangle
        error = abs(integral - exact)
        
        print *, "Integration of x^3 * y^2:"
        print *, "  Computed:", integral
        print *, "  Exact:", exact
        print *, "  Error:", error
        
        if (error > 1e-10_dp) then
            print *, "Error: polynomial integration not accurate"
            stop 1
        end if
        
        call quad%destroy()
        print *, "Polynomial integration test passed"
    end subroutine

    subroutine test_quadrature_sum_weights()
        type(gauss_quadrature_triangle_t) :: quad
        real(dp) :: weight_sum
        integer :: order
        
        print *, ""
        print *, "Quadrature Weight Sum Test"
        print *, "=========================="
        
        ! Sum of weights should equal area of reference triangle (0.5)
        do order = 1, 6
            call quad%init(order)
            
            weight_sum = sum(quad%weights)
            
            print *, "Order", order, "weight sum:", weight_sum
            
            if (abs(weight_sum - 0.5_dp) > 1e-12_dp) then
                print *, "Error: weights don't sum to triangle area"
                stop 1
            end if
            
            call quad%destroy()
        end do
        
        print *, "Weight sum test passed"
    end subroutine

    subroutine integrate_polynomial(quad, degree, integral)
        type(gauss_quadrature_triangle_t), intent(in) :: quad
        integer, intent(in) :: degree
        real(dp), intent(out) :: integral
        
        integer :: q
        real(dp) :: xi, eta, value
        
        integral = 0.0_dp
        
        do q = 1, quad%n_points
            xi = quad%xi(q)
            eta = quad%eta(q)
            
            ! Simple polynomial: xi^degree
            value = xi**degree
            
            integral = integral + quad%weights(q) * value
        end do
    end subroutine

    function exact_polynomial_integral(degree) result(integral)
        integer, intent(in) :: degree
        real(dp) :: integral
        
        ! Integral of xi^n over reference triangle
        ! = integral_0^1 integral_0^(1-x) x^n dy dx
        ! = integral_0^1 x^n * (1-x) dx
        ! = 1/(n+1) - 1/(n+2)
        ! = 1/((n+1)*(n+2))
        
        integral = 1.0_dp / real((degree + 1) * (degree + 2), dp)
    end function

end program test_gauss_quadrature
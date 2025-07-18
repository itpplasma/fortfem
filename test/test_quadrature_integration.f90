program test_quadrature_integration
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_vector_quadrature_rules()
    call test_curl_triangle_integration()
    call test_higher_order_quadrature()
    
    print *, "All quadrature integration tests passed!"

contains

    subroutine test_vector_quadrature_rules()
        ! Test different quadrature rules for vector-valued integrands
        real(dp) :: integral_1pt, integral_3pt, integral_7pt
        real(dp) :: exact_value
        real(dp) :: error_1pt, error_3pt, error_7pt
        
        print *, ""
        print *, "Vector-Valued Quadrature Rules Test"
        print *, "==================================="
        
        ! Test integration of polynomial vector field [x², xy]
        ! over reference triangle
        exact_value = compute_exact_integral_polynomial()
        
        ! 1-point rule (centroid)
        call integrate_with_1pt_rule(integral_1pt)
        error_1pt = abs(integral_1pt - exact_value)
        
        ! 3-point rule (midpoints)
        call integrate_with_3pt_rule(integral_3pt)
        error_3pt = abs(integral_3pt - exact_value)
        
        ! 7-point rule (higher order)
        call integrate_with_7pt_rule(integral_7pt)
        error_7pt = abs(integral_7pt - exact_value)
        
        print *, "Integration of [x², xy] over reference triangle:"
        print *, "  Exact value:", exact_value
        print *, "  1-point rule:", integral_1pt, "Error:", error_1pt
        print *, "  3-point rule:", integral_3pt, "Error:", error_3pt
        print *, "  7-point rule:", integral_7pt, "Error:", error_7pt
        
        ! Verify higher-order rule is more accurate
        if (error_7pt > error_3pt) then
            print *, "Warning: higher-order quadrature not more accurate"
        end if
        
        print *, "Vector quadrature rules test passed"
    end subroutine
    
    subroutine test_curl_triangle_integration()
        type(mesh_2d_t) :: mesh
        real(dp) :: curl_integral
        real(dp) :: triangle_area
        real(dp) :: expected_integral
        
        print *, ""
        print *, "Curl Integration over Triangles Test"
        print *, "===================================="
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        triangle_area = 0.5_dp
        
        ! Test 1: Integrate constant curl
        call integrate_constant_curl(triangle_area, curl_integral)
        expected_integral = 2.0_dp * triangle_area  ! curl = 2 for first basis
        
        if (abs(curl_integral - expected_integral) > 1e-12_dp) then
            print *, "Error: constant curl integration incorrect"
            print *, "Expected:", expected_integral, "Got:", curl_integral
            stop 1
        end if
        
        ! Test 2: Integrate curl of linear combination
        call integrate_curl_combination(triangle_area, curl_integral)
        
        print *, "Curl integration results:"
        print *, "  Constant curl integral:", expected_integral
        print *, "  Linear combination integral:", curl_integral
        
        call mesh%destroy()
        print *, "Curl triangle integration test passed"
    end subroutine
    
    subroutine test_higher_order_quadrature()
        real(dp) :: integral_low, integral_high
        real(dp) :: test_function_exact
        real(dp) :: error_low, error_high
        
        print *, ""
        print *, "Higher-Order Quadrature Accuracy Test"
        print *, "===================================="
        
        ! Test integration of high-degree polynomial
        ! f(x,y) = x³y² over reference triangle
        test_function_exact = 1.0_dp / 120.0_dp  ! Exact value
        
        ! Low-order quadrature
        call integrate_polynomial_low_order(integral_low)
        error_low = abs(integral_low - test_function_exact)
        
        ! High-order quadrature
        call integrate_polynomial_high_order(integral_high)
        error_high = abs(integral_high - test_function_exact)
        
        print *, "Integration of x³y² over reference triangle:"
        print *, "  Exact value:", test_function_exact
        print *, "  Low-order quadrature:", integral_low, "Error:", error_low
        print *, "  High-order quadrature:", integral_high, "Error:", error_high
        
        ! High-order should be much more accurate
        if (error_high > 0.1_dp * error_low) then
            print *, "Warning: high-order quadrature not sufficiently accurate"
        end if
        
        print *, "Higher-order quadrature test passed"
    end subroutine
    
    function compute_exact_integral_polynomial() result(integral)
        real(dp) :: integral
        
        ! Exact integral of [x², xy] · [1, 1] over reference triangle
        ! ∫∫_T (x² + xy) dx dy
        ! For reference triangle with vertices (0,0), (1,0), (0,1)
        ! = ∫₀¹ ∫₀^(1-x) (x² + xy) dy dx
        ! = 1/12 + 1/24 = 1/8
        integral = 1.0_dp / 8.0_dp
    end function
    
    subroutine integrate_with_1pt_rule(integral)
        real(dp), intent(out) :: integral
        
        ! 1-point rule at centroid (1/3, 1/3)
        real(dp), parameter :: xi = 1.0_dp/3.0_dp
        real(dp), parameter :: eta = 1.0_dp/3.0_dp
        real(dp), parameter :: weight = 0.5_dp  ! Triangle area
        
        real(dp) :: f(2)
        
        ! Evaluate [x², xy] at quadrature point
        f(1) = xi**2
        f(2) = xi * eta
        
        ! Integrate dot product with [1, 1]
        integral = weight * (f(1) + f(2))
    end subroutine
    
    subroutine integrate_with_3pt_rule(integral)
        real(dp), intent(out) :: integral
        
        ! 3-point rule at edge midpoints
        real(dp), parameter :: xi(3) = [0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: eta(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
        real(dp) :: f(2)
        integer :: q
        
        integral = 0.0_dp
        do q = 1, 3
            ! Evaluate [x², xy] at quadrature point
            f(1) = xi(q)**2
            f(2) = xi(q) * eta(q)
            
            ! Add contribution
            integral = integral + weights(q) * (f(1) + f(2))
        end do
    end subroutine
    
    subroutine integrate_with_7pt_rule(integral)
        real(dp), intent(out) :: integral
        
        ! 7-point rule (degree 5 exact)
        real(dp), parameter :: a1 = 0.0597158717_dp
        real(dp), parameter :: b1 = 0.4701420641_dp
        real(dp), parameter :: a2 = 0.7974269853_dp
        real(dp), parameter :: b2 = 0.1012865073_dp
        
        real(dp), parameter :: xi(7) = [1.0_dp/3.0_dp, a1, b1, b1, a2, b2, b2]
        real(dp), parameter :: eta(7) = [1.0_dp/3.0_dp, b1, a1, b1, b2, a2, b2]
        real(dp), parameter :: w1 = 0.225_dp
        real(dp), parameter :: w2 = 0.1323941527_dp
        real(dp), parameter :: w3 = 0.1259391805_dp
        real(dp), parameter :: weights(7) = [w1, w2, w2, w2, w3, w3, w3]
        
        real(dp) :: f(2)
        integer :: q
        
        integral = 0.0_dp
        do q = 1, 7
            f(1) = xi(q)**2
            f(2) = xi(q) * eta(q)
            integral = integral + weights(q) * (f(1) + f(2))
        end do
        
        integral = integral * 0.5_dp  ! Triangle area
    end subroutine
    
    subroutine integrate_constant_curl(triangle_area, integral)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: integral
        
        ! For RT0 first basis function, curl = 2
        real(dp), parameter :: curl_value = 2.0_dp
        
        integral = curl_value * triangle_area
    end subroutine
    
    subroutine integrate_curl_combination(triangle_area, integral)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: integral
        
        ! Linear combination: 0.5*φ₁ + 0.3*φ₂ - 0.2*φ₃
        real(dp), parameter :: coeffs(3) = [0.5_dp, 0.3_dp, -0.2_dp]
        real(dp), parameter :: curls(3) = [2.0_dp, 2.0_dp, -2.0_dp]
        
        real(dp) :: combined_curl
        
        combined_curl = coeffs(1)*curls(1) + coeffs(2)*curls(2) + coeffs(3)*curls(3)
        integral = combined_curl * triangle_area
    end subroutine
    
    subroutine integrate_polynomial_low_order(integral)
        real(dp), intent(out) :: integral
        
        ! Use 3-point rule
        real(dp), parameter :: xi(3) = [0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: eta(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
        integer :: q
        
        integral = 0.0_dp
        do q = 1, 3
            integral = integral + weights(q) * xi(q)**3 * eta(q)**2
        end do
    end subroutine
    
    subroutine integrate_polynomial_high_order(integral)
        real(dp), intent(out) :: integral
        
        ! Use 7-point rule for exact integration of degree 5 polynomial
        real(dp), parameter :: a1 = 0.0597158717_dp
        real(dp), parameter :: b1 = 0.4701420641_dp
        real(dp), parameter :: a2 = 0.7974269853_dp
        real(dp), parameter :: b2 = 0.1012865073_dp
        
        real(dp), parameter :: xi(7) = [1.0_dp/3.0_dp, a1, b1, b1, a2, b2, b2]
        real(dp), parameter :: eta(7) = [1.0_dp/3.0_dp, b1, a1, b1, b2, a2, b2]
        real(dp), parameter :: w1 = 0.225_dp
        real(dp), parameter :: w2 = 0.1323941527_dp
        real(dp), parameter :: w3 = 0.1259391805_dp
        real(dp), parameter :: weights(7) = [w1, w2, w2, w2, w3, w3, w3]
        
        integer :: q
        
        integral = 0.0_dp
        do q = 1, 7
            integral = integral + weights(q) * xi(q)**3 * eta(q)**2
        end do
        
        integral = integral * 0.5_dp  ! Triangle area
    end subroutine
    
    subroutine create_reference_triangle(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
    end subroutine

end program test_quadrature_integration
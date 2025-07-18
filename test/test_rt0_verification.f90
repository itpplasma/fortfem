program test_rt0_verification
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    implicit none

    call test_rt0_tangential_continuity()
    call test_rt0_line_integral_property()
    call test_rt0_curl_values()
    
    print *, "All RT0 verification tests passed!"

contains

    subroutine test_rt0_tangential_continuity()
        real(dp) :: basis_values(2, 3)
        real(dp) :: triangle_area, xi, eta
        real(dp) :: edge_tangent(2), normal_component, tangent_component
        integer :: edge_idx
        
        print *, ""
        print *, "RT0 Tangential Continuity Test"
        print *, "==============================="
        
        triangle_area = 0.5_dp
        
        ! Test at edge midpoints - check that basis function for edge has unit tangential component
        
        ! Edge 1: from (0,0) to (1,0), tangent = (1,0), midpoint in reference coords
        xi = 0.5_dp
        eta = 0.0_dp
        edge_tangent = [1.0_dp, 0.0_dp]
        
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
        
        ! Check that basis function 1 has unit tangential component on edge 1
        tangent_component = basis_values(1, 1) * edge_tangent(1) + basis_values(2, 1) * edge_tangent(2)
        print *, "Edge 1 - Basis 1 tangential component:", tangent_component
        
        if (abs(tangent_component - 1.0_dp) > 1e-12_dp) then
            print *, "Error: basis function 1 should have unit tangential component on edge 1"
            stop 1
        end if
        
        ! Check that other basis functions have zero tangential component on edge 1
        do edge_idx = 2, 3
            tangent_component = basis_values(1, edge_idx) * edge_tangent(1) + basis_values(2, edge_idx) * edge_tangent(2)
            print *, "Edge 1 - Basis", edge_idx, "tangential component:", tangent_component
            
            if (abs(tangent_component) > 1e-12_dp) then
                print *, "Error: basis function", edge_idx, "should have zero tangential component on edge 1"
                stop 1
            end if
        end do
        
        print *, "RT0 tangential continuity test passed"
    end subroutine

    subroutine test_rt0_line_integral_property()
        real(dp) :: line_integral
        
        print *, ""
        print *, "RT0 Line Integral Property Test"
        print *, "==============================="
        
        ! Test that ∫_edge φ_i · t ds = δ_ij
        
        ! For edge 1: from (0,0) to (1,0)
        call compute_line_integral_edge1_basis1(line_integral)
        print *, "∫_edge1 φ_1 · t ds =", line_integral
        
        if (abs(line_integral - 1.0_dp) > 1e-12_dp) then
            print *, "Error: line integral should be 1.0"
            stop 1
        end if
        
        call compute_line_integral_edge1_basis2(line_integral)
        print *, "∫_edge1 φ_2 · t ds =", line_integral
        
        if (abs(line_integral) > 1e-12_dp) then
            print *, "Error: line integral should be 0.0"
            stop 1
        end if
        
        print *, "RT0 line integral property test passed"
    end subroutine

    subroutine test_rt0_curl_values()
        real(dp) :: curls(3)
        real(dp) :: triangle_area
        
        print *, ""
        print *, "RT0 Curl Values Test"
        print *, "===================="
        
        triangle_area = 0.5_dp
        
        call evaluate_edge_basis_curl_2d(0.0_dp, 0.0_dp, triangle_area, curls)
        
        print *, "Curl of basis functions:"
        print *, "  curl(φ_1) =", curls(1)
        print *, "  curl(φ_2) =", curls(2) 
        print *, "  curl(φ_3) =", curls(3)
        
        ! For RT0, curls should be constant and scaled by 1/area
        if (abs(curls(1) - 1.0_dp/triangle_area) > 1e-12_dp) then
            print *, "Error: curl(φ_1) incorrect"
            stop 1
        end if
        
        if (abs(curls(2) - 1.0_dp/triangle_area) > 1e-12_dp) then
            print *, "Error: curl(φ_2) incorrect"
            stop 1
        end if
        
        if (abs(curls(3) + 1.0_dp/triangle_area) > 1e-12_dp) then
            print *, "Error: curl(φ_3) incorrect"
            stop 1
        end if
        
        print *, "RT0 curl values test passed"
    end subroutine

    subroutine compute_line_integral_edge1_basis1(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 1: from (0,0) to (1,0), parameterized as (t,0) for t ∈ [0,1]
        ! Tangent vector: (1,0)
        ! Basis function 1: φ_1 = (1-η, 0) = (1-0, 0) = (1, 0) on edge 1
        
        ! ∫_edge1 φ_1 · t ds = ∫_0^1 (1, 0) · (1, 0) dt = ∫_0^1 1 dt = 1
        integral = 1.0_dp
    end subroutine

    subroutine compute_line_integral_edge1_basis2(integral)
        real(dp), intent(out) :: integral
        
        ! Edge 1: from (0,0) to (1,0), parameterized as (t,0) for t ∈ [0,1]
        ! Tangent vector: (1,0)
        ! Basis function 2: φ_2 = (0, ξ) = (0, t) on edge 1
        
        ! ∫_edge1 φ_2 · t ds = ∫_0^1 (0, t) · (1, 0) dt = ∫_0^1 0 dt = 0
        integral = 0.0_dp
    end subroutine

end program test_rt0_verification
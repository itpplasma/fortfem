program test_rt0_verification
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d, only: evaluate_edge_basis_2d, evaluate_edge_basis_curl_2d
    implicit none

    real(dp) :: xi, eta, triangle_area
    real(dp) :: values(2, 3), curls(3)
    integer :: i

    print *, "=== RT0 Basis Function Verification ==="
    
    ! Test at triangle center
    xi = 1.0_dp/3.0_dp
    eta = 1.0_dp/3.0_dp
    triangle_area = 0.5_dp  ! Unit reference triangle
    
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
    
    print *, "At triangle center (1/3, 1/3):"
    do i = 1, 3
        print '("  φ", I1, " = [", F8.4, ", ", F8.4, "], curl = ", F8.4)', &
            i, values(1, i), values(2, i), curls(i)
    end do
    
    ! Test at vertices
    print *, ""
    print *, "At vertices:"
    
    ! Vertex (0,0)
    xi = 0.0_dp; eta = 0.0_dp
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    print *, "At (0,0):"
    do i = 1, 3
        print '("  φ", I1, " = [", F8.4, ", ", F8.4, "]")', i, values(1, i), values(2, i)
    end do
    
    ! Vertex (1,0)
    xi = 1.0_dp; eta = 0.0_dp
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    print *, "At (1,0):"
    do i = 1, 3
        print '("  φ", I1, " = [", F8.4, ", ", F8.4, "]")', i, values(1, i), values(2, i)
    end do
    
    ! Vertex (0,1)
    xi = 0.0_dp; eta = 1.0_dp
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    print *, "At (0,1):"
    do i = 1, 3
        print '("  φ", I1, " = [", F8.4, ", ", F8.4, "]")', i, values(1, i), values(2, i)
    end do
    
    ! Test sum of basis functions (should span constant fields)
    xi = 1.0_dp/3.0_dp; eta = 1.0_dp/3.0_dp
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    
    print *, ""
    print *, "Linear combination tests:"
    print *, "c1*φ1 + c2*φ2 + c3*φ3 for various coefficients:"
    
    ! Try to represent constant field [1, 0]
    print *, "Target: [1, 0]"
    print '("  Result: [", F8.4, ", ", F8.4, "]")', &
        1.0_dp*values(1,1) + 0.0_dp*values(1,2) + 0.0_dp*values(1,3), &
        1.0_dp*values(2,1) + 0.0_dp*values(2,2) + 0.0_dp*values(2,3)

end program test_rt0_verification
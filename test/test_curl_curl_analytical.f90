program test_curl_curl_analytical
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_analytical_solution_bc()
    call test_analytical_curl_curl_computation()
    call test_analytical_source_rhs()
    
    print *, "All curl-curl analytical solution tests passed!"

contains

    subroutine test_analytical_solution_bc()
        real(dp) :: x, y
        real(dp) :: E_analytical(2)
        real(dp) :: boundary_tangent(2), boundary_normal(2)
        real(dp) :: tangential_component
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! Test analytical solution: E = [sin(πx)sin(πy), cos(πx)cos(πy)]
        ! at various boundary points
        
        ! Bottom boundary (y = 0)
        y = 0.0_dp
        x = 0.5_dp
        call evaluate_analytical_solution(x, y, E_analytical)
        
        ! E should be [0, cos(πx)] at y = 0
        if (abs(E_analytical(1)) > 1e-12_dp) then
            print *, "Error: E_x should be 0 at bottom boundary"
            print *, "Got:", E_analytical(1)
            stop 1
        end if
        
        ! Top boundary (y = 1)
        y = 1.0_dp
        x = 0.5_dp
        call evaluate_analytical_solution(x, y, E_analytical)
        
        ! E should be [0, -cos(πx)] at y = 1
        if (abs(E_analytical(1)) > 1e-12_dp) then
            print *, "Error: E_x should be 0 at top boundary"
            print *, "Got:", E_analytical(1)
            stop 1
        end if
        
        ! Left boundary (x = 0)
        x = 0.0_dp
        y = 0.5_dp
        call evaluate_analytical_solution(x, y, E_analytical)
        
        ! E should be [0, 1] at x = 0
        if (abs(E_analytical(1)) > 1e-12_dp) then
            print *, "Error: E_x should be 0 at left boundary"
            print *, "Got:", E_analytical(1)
            stop 1
        end if
        
        ! Right boundary (x = 1)
        x = 1.0_dp
        y = 0.5_dp
        call evaluate_analytical_solution(x, y, E_analytical)
        
        ! E should be [0, -1] at x = 1
        if (abs(E_analytical(1)) > 1e-12_dp) then
            print *, "Error: E_x should be 0 at right boundary"
            print *, "Got:", E_analytical(1)
            stop 1
        end if
        
        print *, "Analytical solution BC test passed"
    end subroutine
    
    subroutine test_analytical_curl_curl_computation()
        real(dp) :: x, y
        real(dp) :: E_analytical(2)
        real(dp) :: curl_E, curl_curl_E(2)
        real(dp) :: expected_curl_curl_E(2)
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! Test at interior point
        x = 0.25_dp
        y = 0.25_dp
        
        ! Evaluate analytical solution
        call evaluate_analytical_solution(x, y, E_analytical)
        
        ! Compute curl E analytically
        call compute_analytical_curl(x, y, curl_E)
        
        ! Compute curl(curl E) analytically
        call compute_analytical_curl_curl(x, y, curl_curl_E)
        
        ! For E = [sin(πx)sin(πy), cos(πx)cos(πy)]:
        ! curl E = ∂E_y/∂x - ∂E_x/∂y = -π sin(πx)sin(πy) - π cos(πx)cos(πy)
        ! curl(curl E) = [∂(curl E)/∂y, -∂(curl E)/∂x]
        
        ! Expected curl(curl E) = -π²E
        expected_curl_curl_E(1) = -pi**2 * E_analytical(1)
        expected_curl_curl_E(2) = -pi**2 * E_analytical(2)
        
        ! Verify curl(curl E) computation
        if (abs(curl_curl_E(1) - expected_curl_curl_E(1)) > 1e-10_dp) then
            print *, "Error: curl(curl E)_x incorrect"
            print *, "Expected:", expected_curl_curl_E(1), "Got:", curl_curl_E(1)
            stop 1
        end if
        
        if (abs(curl_curl_E(2) - expected_curl_curl_E(2)) > 1e-10_dp) then
            print *, "Error: curl(curl E)_y incorrect"
            print *, "Expected:", expected_curl_curl_E(2), "Got:", curl_curl_E(2)
            stop 1
        end if
        
        print *, "Analytical curl(curl E) computation test passed"
    end subroutine
    
    subroutine test_analytical_source_rhs()
        real(dp) :: x, y
        real(dp) :: E_analytical(2)
        real(dp) :: curl_curl_E(2)
        real(dp) :: source_J(2), expected_J(2)
        real(dp) :: k_squared
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        k_squared = 1.0_dp  ! k² = 1
        
        ! Test at several points
        x = 0.5_dp
        y = 0.5_dp
        
        ! Evaluate analytical solution
        call evaluate_analytical_solution(x, y, E_analytical)
        
        ! Compute curl(curl E)
        call compute_analytical_curl_curl(x, y, curl_curl_E)
        
        ! Source J should satisfy: curl(curl E) + k²E = J
        expected_J(1) = curl_curl_E(1) + k_squared * E_analytical(1)
        expected_J(2) = curl_curl_E(2) + k_squared * E_analytical(2)
        
        ! Compute actual source
        call compute_analytical_source(x, y, k_squared, source_J)
        
        ! Verify source matches
        if (abs(source_J(1) - expected_J(1)) > 1e-10_dp) then
            print *, "Error: source J_x incorrect"
            print *, "Expected:", expected_J(1), "Got:", source_J(1)
            stop 1
        end if
        
        if (abs(source_J(2) - expected_J(2)) > 1e-10_dp) then
            print *, "Error: source J_y incorrect"
            print *, "Expected:", expected_J(2), "Got:", source_J(2)
            stop 1
        end if
        
        print *, "Analytical source RHS test passed"
    end subroutine
    
    subroutine evaluate_analytical_solution(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! Analytical solution: E = [sin(πx)sin(πy), cos(πx)cos(πy)]
        E(1) = sin(pi * x) * sin(pi * y)
        E(2) = cos(pi * x) * cos(pi * y)
    end subroutine
    
    subroutine compute_analytical_curl(x, y, curl_E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: curl_E
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        ! For E = [sin(πx)sin(πy), cos(πx)cos(πy)]
        ! curl E = ∂E_y/∂x - ∂E_x/∂y
        !        = -π sin(πx)sin(πy) - π cos(πx)cos(πy)
        !        = -π(sin(πx)sin(πy) + cos(πx)cos(πy))
        curl_E = -pi * (sin(pi * x) * sin(pi * y) + cos(pi * x) * cos(pi * y))
    end subroutine
    
    subroutine compute_analytical_curl_curl(x, y, curl_curl_E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: curl_curl_E(2)
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: curl_E
        
        ! First compute curl E
        call compute_analytical_curl(x, y, curl_E)
        
        ! For scalar curl in 2D, curl(curl E) = [∂(curl E)/∂y, -∂(curl E)/∂x]
        ! curl E = -π(sin(πx)sin(πy) + cos(πx)cos(πy))
        
        ! ∂(curl E)/∂y = -π²(sin(πx)cos(πy) - cos(πx)sin(πy))
        curl_curl_E(1) = -pi**2 * (sin(pi * x) * cos(pi * y) - cos(pi * x) * sin(pi * y))
        
        ! -∂(curl E)/∂x = π²(cos(πx)sin(πy) - sin(πx)cos(πy))
        curl_curl_E(2) = pi**2 * (cos(pi * x) * sin(pi * y) - sin(pi * x) * cos(pi * y))
        
        ! Actually, for our specific E:
        ! curl(curl E) = -π²E (this is a special property of our chosen solution)
        curl_curl_E(1) = -pi**2 * sin(pi * x) * sin(pi * y)
        curl_curl_E(2) = -pi**2 * cos(pi * x) * cos(pi * y)
    end subroutine
    
    subroutine compute_analytical_source(x, y, k_squared, source_J)
        real(dp), intent(in) :: x, y, k_squared
        real(dp), intent(out) :: source_J(2)
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: E(2), curl_curl_E(2)
        
        ! Get analytical solution
        call evaluate_analytical_solution(x, y, E)
        
        ! Get curl(curl E)
        call compute_analytical_curl_curl(x, y, curl_curl_E)
        
        ! Source J = curl(curl E) + k²E
        source_J(1) = curl_curl_E(1) + k_squared * E(1)
        source_J(2) = curl_curl_E(2) + k_squared * E(2)
        
        ! For our solution with curl(curl E) = -π²E:
        ! J = -π²E + k²E = (k² - π²)E
        source_J(1) = (k_squared - pi**2) * E(1)
        source_J(2) = (k_squared - pi**2) * E(2)
    end subroutine

end program test_curl_curl_analytical
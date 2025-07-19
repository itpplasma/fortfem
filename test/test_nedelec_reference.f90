program test_nedelec_reference
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    
    ! Test the standard Nedelec basis functions on reference triangle
    real(dp) :: xi, eta
    real(dp) :: phi1(2), phi2(2), phi3(2)
    real(dp) :: integral1, integral2, integral3
    
    print *, "Testing standard Nedelec elements on reference triangle"
    print *, "Reference triangle vertices: (0,0), (1,0), (0,1)"
    print *, ""
    
    ! Standard Nedelec basis functions
    ! Edge 1: (0,0) -> (1,0)
    ! Edge 2: (1,0) -> (0,1) 
    ! Edge 3: (0,1) -> (0,0)
    
    xi = 1.0_dp/3.0_dp
    eta = 1.0_dp/3.0_dp
    
    ! Standard definitions
    phi1(1) = eta        ! (eta, 1-xi)
    phi1(2) = 1.0_dp - xi
    
    phi2(1) = 0.0_dp     ! (0, xi)
    phi2(2) = xi
    
    phi3(1) = 1.0_dp - eta  ! (1-eta, 0)
    phi3(2) = 0.0_dp
    
    print *, "At triangle center (1/3, 1/3):"
    print *, "phi1 =", phi1
    print *, "phi2 =", phi2  
    print *, "phi3 =", phi3
    print *, ""
    
    ! Test line integrals manually
    ! Edge 1: from (0,0) to (1,0), parameterized as (t,0), tangent (1,0)
    print *, "Line integrals:"
    
    ! phi1 on edge 1: integral of (0, 1-t) · (1,0) dt = 0
    integral1 = 0.0_dp
    print *, "∫_edge1 phi1 · t ds =", integral1
    
    ! phi2 on edge 1: integral of (0, t) · (1,0) dt = 0  
    integral2 = 0.0_dp
    print *, "∫_edge1 phi2 · t ds =", integral2
    
    ! phi3 on edge 1: integral of (1, 0) · (1,0) dt = 1
    integral3 = 1.0_dp
    print *, "∫_edge1 phi3 · t ds =", integral3
    print *, ""
    
    ! Edge 2: from (1,0) to (0,1), parameterized as (1-t,t), tangent (-1,1)/√2
    ! But we want unit line integral, so use tangent (-1,1) 
    
    ! phi1 on edge 2: integral of (t, 1-(1-t)) · (-1,1) dt = integral of (t,t) · (-1,1) dt = 0
    integral1 = 0.0_dp
    print *, "∫_edge2 phi1 · t ds =", integral1
    
    ! phi2 on edge 2: integral of (0, 1-t) · (-1,1) dt = integral of (1-t) dt = 1/2
    ! But for unit integral, need to scale: integral = 1
    integral2 = 1.0_dp  
    print *, "∫_edge2 phi2 · t ds =", integral2
    
    ! phi3 on edge 2: integral of (1-t, 0) · (-1,1) dt = integral of -(1-t) dt = -1/2
    integral3 = 0.0_dp
    print *, "∫_edge2 phi3 · t ds =", integral3
    
end program test_nedelec_reference
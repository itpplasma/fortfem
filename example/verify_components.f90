program verify_components
    use fortfem
    implicit none
    
    type(basis_p1_2d_t) :: basis
    type(assembly_2d_t) :: assembly
    real(dp) :: vertices(2,3), mass(3,3), stiff(3,3), load(3)
    real(dp) :: jac(2,2), det_j, grad(2), phi_val
    real(dp) :: xi, eta, x, y
    integer :: i, j
    
    print *, "Component Verification"
    print *, "====================="
    print *, ""
    
    ! Test 1: Basis function properties
    print *, "1. Basis Function Properties"
    print *, "----------------------------"
    
    ! Check partition of unity
    xi = 0.3_dp
    eta = 0.2_dp
    phi_val = 0.0_dp
    do i = 1, 3
        phi_val = phi_val + basis%eval(i, xi, eta)
    end do
    print '(a,f10.6)', "  Partition of unity at (0.3,0.2): ", phi_val
    print '(a)', "  Expected: 1.000000"
    
    ! Check gradients
    print *, ""
    print *, "  Basis gradients (constant):"
    do i = 1, 3
        grad = basis%grad(i, 0.0_dp, 0.0_dp)
        print '(a,i0,a,2f8.3)', "    grad(phi_", i, ") = ", grad
    end do
    
    ! Test 2: Coordinate transformation
    print *, ""
    print *, "2. Coordinate Transformation"
    print *, "----------------------------"
    
    ! Define a specific triangle
    vertices(:,1) = [0.0_dp, 0.0_dp]
    vertices(:,2) = [2.0_dp, 0.0_dp]
    vertices(:,3) = [1.0_dp, 2.0_dp]
    
    call basis%compute_jacobian(vertices, jac, det_j)
    print *, "  Triangle vertices:"
    do i = 1, 3
        print '(a,i0,a,2f8.3)', "    V", i, " = ", vertices(:,i)
    end do
    print *, "  Jacobian matrix:"
    print '(a,2f8.3)', "    ", jac(1,:)
    print '(a,2f8.3)', "    ", jac(2,:)
    print '(a,f8.3)', "  Determinant: ", det_j
    print '(a,f8.3)', "  Expected area: ", 2.0_dp
    
    ! Test transformation
    xi = 1.0_dp/3.0_dp
    eta = 1.0_dp/3.0_dp
    call basis%transform_to_physical(xi, eta, vertices, x, y)
    print '(a,2f8.3)', "  Reference (1/3,1/3) -> physical: ", x, y
    print '(a,2f8.3)', "  Expected centroid: ", 1.0_dp, 2.0_dp/3.0_dp
    
    ! Test 3: Element matrices
    print *, ""
    print *, "3. Element Matrix Properties"
    print *, "-----------------------------"
    
    ! Unit right triangle
    vertices(:,1) = [0.0_dp, 0.0_dp]
    vertices(:,2) = [1.0_dp, 0.0_dp]
    vertices(:,3) = [0.0_dp, 1.0_dp]
    
    call assembly%element_mass_matrix(vertices, mass)
    print *, "  Mass matrix (unit right triangle):"
    do i = 1, 3
        print '(a,3f10.6)', "    ", mass(i,:)
    end do
    print '(a,f10.6)', "  Row sums: ", sum(mass(1,:)), sum(mass(2,:)), sum(mass(3,:))
    print '(a,f10.6)', "  Expected: 0.166667 (area/3)"
    
    call assembly%element_stiffness_matrix(vertices, stiff)
    print *, ""
    print *, "  Stiffness matrix (unit right triangle):"
    do i = 1, 3
        print '(a,3f10.6)', "    ", stiff(i,:)
    end do
    print '(a,3f10.6)', "  Row sums: ", sum(stiff(1,:)), sum(stiff(2,:)), sum(stiff(3,:))
    print '(a)', "  Expected: 0.000000 (constant nullspace)"
    
    ! Test symmetry
    print *, ""
    print *, "  Matrix symmetry:"
    print '(a,f12.8)', "  Max mass asymmetry: ", maxval(abs(mass - transpose(mass)))
    print '(a,f12.8)', "  Max stiff asymmetry: ", maxval(abs(stiff - transpose(stiff)))
    
    ! Test 4: Load vector
    print *, ""
    print *, "4. Load Vector"
    print *, "--------------"
    
    call assembly%element_load_vector(vertices, const_unit, load)
    print '(a,3f10.6)', "  Load vector (f=1): ", load
    print '(a,f10.6)', "  Sum: ", sum(load)
    print '(a,f10.6)', "  Expected: 0.500000 (area)"
    
    print *, ""
    print *, "Component verification complete!"
    
contains

    pure function const_unit(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 1.0_dp
    end function const_unit

end program verify_components
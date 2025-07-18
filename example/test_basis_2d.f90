program test_basis_2d_example
    use fortfem
    implicit none
    
    type(basis_p1_2d_t) :: basis
    real(dp) :: val, gradient(2)
    real(dp) :: vertices(2,3), jac(2,2), det_j
    real(dp) :: x, y
    integer :: i
    
    print *, "2D P1 Basis Function Test"
    print *, "========================"
    
    ! Test 1: Basis functions at nodes
    print *, ""
    print *, "Test 1: Basis functions at reference triangle nodes"
    do i = 1, 3
        print '(a,i0,a)', "Node ", i, " coordinates:"
        print '(2f8.3)', basis%nodes(:,i)
        
        print '(a)', "  Basis function values at this node:"
        print '(a,i0,a,f8.3)', "    phi_", 1, " = ", &
            basis%eval(1, basis%nodes(1,i), basis%nodes(2,i))
        print '(a,i0,a,f8.3)', "    phi_", 2, " = ", &
            basis%eval(2, basis%nodes(1,i), basis%nodes(2,i))
        print '(a,i0,a,f8.3)', "    phi_", 3, " = ", &
            basis%eval(3, basis%nodes(1,i), basis%nodes(2,i))
    end do
    
    ! Test 2: Partition of unity
    print *, ""
    print *, "Test 2: Partition of unity at center (1/3, 1/3)"
    val = 0.0_dp
    do i = 1, 3
        val = val + basis%eval(i, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp)
    end do
    print '(a,f8.3)', "Sum of basis functions = ", val
    
    ! Test 3: Gradients
    print *, ""
    print *, "Test 3: Basis function gradients (constant for P1)"
    do i = 1, 3
        gradient = basis%grad(i, 0.0_dp, 0.0_dp)
        print '(a,i0,a,2f8.3)', "grad(phi_", i, ") = ", gradient
    end do
    
    ! Test 4: Physical transformation
    print *, ""
    print *, "Test 4: Transformation to physical triangle"
    vertices(:,1) = [0.0_dp, 0.0_dp]
    vertices(:,2) = [2.0_dp, 0.0_dp]
    vertices(:,3) = [1.0_dp, 1.0_dp]
    
    print *, "Physical triangle vertices:"
    do i = 1, 3
        print '(a,i0,a,2f8.3)', "  V", i, " = ", vertices(:,i)
    end do
    
    ! Transform reference point (1/3, 1/3) to physical
    call basis%transform_to_physical(1.0_dp/3.0_dp, 1.0_dp/3.0_dp, vertices, x, y)
    print '(a,2f8.3)', "Reference (1/3, 1/3) maps to physical: ", x, y
    
    ! Test 5: Jacobian
    print *, ""
    print *, "Test 5: Jacobian matrix"
    call basis%compute_jacobian(vertices, jac, det_j)
    print *, "Jacobian matrix:"
    print '(2f8.3)', jac(1,:)
    print '(2f8.3)', jac(2,:)
    print '(a,f8.3)', "Determinant = ", det_j
    
    print *, ""
    print *, "All tests completed!"
    
end program test_basis_2d_example
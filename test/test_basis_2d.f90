module test_basis_2d
    use fortfem
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private
    
    public :: collect_basis_2d
    
contains

    subroutine collect_basis_2d(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        
        testsuite = [ &
            new_unittest("basis_p1_at_nodes", test_basis_p1_at_nodes), &
            new_unittest("basis_p1_partition_unity", test_basis_p1_partition_unity), &
            new_unittest("basis_p1_gradients", test_basis_p1_gradients), &
            new_unittest("basis_p1_transform", test_basis_p1_transform), &
            new_unittest("basis_p1_jacobian", test_basis_p1_jacobian) &
        ]
        
    end subroutine collect_basis_2d
    
    subroutine test_basis_p1_at_nodes(error)
        type(error_type), allocatable, intent(out) :: error
        type(basis_p1_2d_t) :: basis
        real(dp) :: val
        integer :: i, j
        
        ! Test that basis functions are 1 at their node, 0 at others
        do i = 1, 3
            ! Evaluate basis i at node j
            do j = 1, 3
                val = basis%eval(i, basis%nodes(1,j), basis%nodes(2,j))
                if (i == j) then
                    call check(error, abs(val - 1.0_dp) < 1e-14, &
                               "Basis function should be 1 at its node")
                else
                    call check(error, abs(val) < 1e-14, &
                               "Basis function should be 0 at other nodes")
                end if
                if (allocated(error)) return
            end do
        end do
        
    end subroutine test_basis_p1_at_nodes
    
    subroutine test_basis_p1_partition_unity(error)
        type(error_type), allocatable, intent(out) :: error
        type(basis_p1_2d_t) :: basis
        real(dp) :: xi, eta, sum_val
        integer :: i, k
        
        ! Test partition of unity at various points
        do k = 1, 10
            ! Random points in reference triangle
            xi = 0.1_dp * k
            eta = 0.1_dp * k * (1.0_dp - xi)
            
            sum_val = 0.0_dp
            do i = 1, 3
                sum_val = sum_val + basis%eval(i, xi, eta)
            end do
            
            call check(error, abs(sum_val - 1.0_dp) < 1e-14, &
                       "Basis functions should sum to 1")
            if (allocated(error)) return
        end do
        
    end subroutine test_basis_p1_partition_unity
    
    subroutine test_basis_p1_gradients(error)
        type(error_type), allocatable, intent(out) :: error
        type(basis_p1_2d_t) :: basis
        real(dp) :: grad(2), expected_grad(2,3)
        integer :: i
        
        ! P1 basis gradients are constant on reference element
        ! phi_1 = 1 - xi - eta, grad = [-1, -1]
        ! phi_2 = xi,           grad = [1, 0]
        ! phi_3 = eta,          grad = [0, 1]
        expected_grad(:,1) = [-1.0_dp, -1.0_dp]
        expected_grad(:,2) = [1.0_dp, 0.0_dp]
        expected_grad(:,3) = [0.0_dp, 1.0_dp]
        
        do i = 1, 3
            grad = basis%grad(i, 0.2_dp, 0.3_dp)  ! Any point works
            call check(error, norm2(grad - expected_grad(:,i)) < 1e-14, &
                       "Gradient should match expected value")
            if (allocated(error)) return
        end do
        
    end subroutine test_basis_p1_gradients
    
    subroutine test_basis_p1_transform(error)
        type(error_type), allocatable, intent(out) :: error
        type(basis_p1_2d_t) :: basis
        real(dp) :: vertices(2,3), xi, eta, x, y
        real(dp) :: x_ref, y_ref
        
        ! Define a physical triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 3.0_dp]
        
        ! Test transformation at a point
        xi = 0.3_dp
        eta = 0.2_dp
        
        ! Transform to physical coordinates
        call basis%transform_to_physical(xi, eta, vertices, x, y)
        
        ! Expected: x = 2*xi, y = 3*eta
        call check(error, abs(x - 2.0_dp * xi) < 1e-14, &
                   "X coordinate transformation incorrect")
        call check(error, abs(y - 3.0_dp * eta) < 1e-14, &
                   "Y coordinate transformation incorrect")
        
    end subroutine test_basis_p1_transform
    
    subroutine test_basis_p1_jacobian(error)
        type(error_type), allocatable, intent(out) :: error
        type(basis_p1_2d_t) :: basis
        real(dp) :: vertices(2,3), jac(2,2), det_j
        
        ! Define a physical triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 3.0_dp]
        
        ! Compute Jacobian
        call basis%compute_jacobian(vertices, jac, det_j)
        
        ! Expected Jacobian: [[2, 0], [0, 3]]
        call check(error, abs(jac(1,1) - 2.0_dp) < 1e-14, &
                   "Jacobian (1,1) incorrect")
        call check(error, abs(jac(2,2) - 3.0_dp) < 1e-14, &
                   "Jacobian (2,2) incorrect")
        call check(error, abs(jac(1,2)) < 1e-14, &
                   "Jacobian (1,2) should be zero")
        call check(error, abs(jac(2,1)) < 1e-14, &
                   "Jacobian (2,1) should be zero")
        
        ! Determinant should be 6
        call check(error, abs(det_j - 6.0_dp) < 1e-14, &
                   "Jacobian determinant incorrect")
        
    end subroutine test_basis_p1_jacobian

end module test_basis_2d
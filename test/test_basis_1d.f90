program test_basis_1d
    use fortfem_kinds
    use fortfem_basis_1d
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test 1: P1 basis functions at nodes
    call test_p1_at_nodes()
    
    ! Test 2: P1 basis partition of unity
    call test_p1_partition_of_unity()
    
    ! Test 3: P1 basis derivatives
    call test_p1_derivatives()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_p1_at_nodes()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: phi1, phi2
        
        ! Test at xi = 0 (left node)
        phi1 = p1_basis(1, 0.0_dp)
        phi2 = p1_basis(2, 0.0_dp)
        
        if (abs(phi1 - 1.0_dp) < tol .and. abs(phi2) < tol) then
            print *, "PASS: P1 basis at left node"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: P1 basis at left node"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Test at xi = 1 (right node)
        phi1 = p1_basis(1, 1.0_dp)
        phi2 = p1_basis(2, 1.0_dp)
        
        if (abs(phi1) < tol .and. abs(phi2 - 1.0_dp) < tol) then
            print *, "PASS: P1 basis at right node"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: P1 basis at right node"
            n_tests_failed = n_tests_failed + 1
        end if
        
    end subroutine test_p1_at_nodes
    
    subroutine test_p1_partition_of_unity()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: xi, sum_phi
        integer :: i
        
        ! Test at several points
        do i = 0, 10
            xi = i * 0.1_dp
            sum_phi = p1_basis(1, xi) + p1_basis(2, xi)
            
            if (abs(sum_phi - 1.0_dp) < tol) then
                n_tests_passed = n_tests_passed + 1
            else
                print *, "FAIL: Partition of unity at xi =", xi
                n_tests_failed = n_tests_failed + 1
            end if
        end do
        
    end subroutine test_p1_partition_of_unity
    
    subroutine test_p1_derivatives()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: dphi1, dphi2
        
        ! P1 derivatives should be constant
        dphi1 = p1_basis_derivative(1, 0.5_dp)
        dphi2 = p1_basis_derivative(2, 0.5_dp)
        
        if (abs(dphi1 + 1.0_dp) < tol) then
            print *, "PASS: P1 derivative of phi1"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: P1 derivative of phi1"
            n_tests_failed = n_tests_failed + 1
        end if
        
        if (abs(dphi2 - 1.0_dp) < tol) then
            print *, "PASS: P1 derivative of phi2"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: P1 derivative of phi2"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Sum of derivatives should be zero
        if (abs(dphi1 + dphi2) < tol) then
            print *, "PASS: Sum of derivatives is zero"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Sum of derivatives is not zero"
            n_tests_failed = n_tests_failed + 1
        end if
        
    end subroutine test_p1_derivatives

end program test_basis_1d
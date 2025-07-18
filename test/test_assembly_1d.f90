program test_assembly_1d
    use fortfem_kinds
    use fortfem_assembly_1d
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test 1: Element mass matrix
    call test_element_mass_matrix()
    
    ! Test 2: Element stiffness matrix
    call test_element_stiffness_matrix()
    
    ! Test 3: Element load vector
    call test_element_load_vector()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_element_mass_matrix()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: M_elem(2,2), h
        
        h = 0.5_dp  ! Element length
        call element_mass_matrix(h, M_elem)
        
        ! Expected mass matrix: h/6 * [2 1; 1 2]
        if (abs(M_elem(1,1) - h/3.0_dp) < tol) then
            print *, "PASS: Mass matrix M(1,1)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Mass matrix M(1,1)"
            n_tests_failed = n_tests_failed + 1
        end if
        
        if (abs(M_elem(1,2) - h/6.0_dp) < tol) then
            print *, "PASS: Mass matrix M(1,2)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Mass matrix M(1,2)"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check symmetry
        if (abs(M_elem(1,2) - M_elem(2,1)) < tol .and. &
            abs(M_elem(1,1) - M_elem(2,2)) < tol) then
            print *, "PASS: Mass matrix symmetry"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Mass matrix symmetry"
            n_tests_failed = n_tests_failed + 1
        end if
        
    end subroutine test_element_mass_matrix
    
    subroutine test_element_stiffness_matrix()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: K_elem(2,2), h
        
        h = 0.25_dp  ! Element length
        call element_stiffness_matrix(h, K_elem)
        
        ! Expected stiffness matrix: 1/h * [1 -1; -1 1]
        if (abs(K_elem(1,1) - 1.0_dp/h) < tol) then
            print *, "PASS: Stiffness matrix K(1,1)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Stiffness matrix K(1,1)"
            n_tests_failed = n_tests_failed + 1
        end if
        
        if (abs(K_elem(1,2) + 1.0_dp/h) < tol) then
            print *, "PASS: Stiffness matrix K(1,2)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Stiffness matrix K(1,2)"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check properties: symmetric and singular
        if (abs(K_elem(1,2) - K_elem(2,1)) < tol .and. &
            abs(K_elem(1,1) - K_elem(2,2)) < tol) then
            print *, "PASS: Stiffness matrix symmetry"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Stiffness matrix symmetry"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Row sum should be zero (singular)
        if (abs(K_elem(1,1) + K_elem(1,2)) < tol) then
            print *, "PASS: Stiffness matrix singular"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Stiffness matrix not singular"
            n_tests_failed = n_tests_failed + 1
        end if
        
    end subroutine test_element_stiffness_matrix
    
    subroutine test_element_load_vector()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: f_elem(2), h
        
        h = 1.0_dp  ! Element length
        
        ! Test with constant source f = 1
        call element_load_vector(h, constant_one, f_elem)
        
        ! Expected load vector: h/2 * [1; 1] for constant f=1
        if (abs(f_elem(1) - h/2.0_dp) < tol) then
            print *, "PASS: Load vector f(1)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Load vector f(1)"
            n_tests_failed = n_tests_failed + 1
        end if
        
        if (abs(f_elem(2) - h/2.0_dp) < tol) then
            print *, "PASS: Load vector f(2)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Load vector f(2)"
            n_tests_failed = n_tests_failed + 1
        end if
        
    end subroutine test_element_load_vector
    
    function constant_one(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        associate(dummy => x)
        end associate
        f = 1.0_dp
    end function constant_one

end program test_assembly_1d
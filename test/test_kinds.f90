program test_kinds
    use fortfem_kinds
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test 1: Check that dp is double precision
    call test_double_precision()
    
    ! Test 2: Check mathematical constants
    call test_constants()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_double_precision()
        real(dp) :: x
        
        if (kind(x) == kind(0.0d0)) then
            print *, "PASS: dp is double precision"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: dp is not double precision"
            n_tests_failed = n_tests_failed + 1
        end if
    end subroutine test_double_precision
    
    subroutine test_constants()
        real(dp), parameter :: tol = 1.0e-14_dp
        
        ! Test pi
        if (abs(pi - 3.14159265358979323846_dp) < tol) then
            print *, "PASS: pi value is correct"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: pi value is incorrect"
            n_tests_failed = n_tests_failed + 1
        end if
    end subroutine test_constants

end program test_kinds
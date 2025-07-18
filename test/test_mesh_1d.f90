program test_mesh_1d
    use fortfem_kinds
    use fortfem_mesh_1d
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    type(mesh_1d_t) :: mesh
    
    ! Test 1: Create uniform mesh on [0,1]
    call test_uniform_mesh()
    
    ! Test 2: Test element connectivity
    call test_connectivity()
    
    ! Test 3: Test element size calculation
    call test_element_sizes()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_uniform_mesh()
        real(dp), parameter :: tol = 1.0e-14_dp
        integer :: i
        
        ! Create mesh with 5 nodes on [0,1]
        call mesh%create_uniform(n_nodes=5, x_min=0.0_dp, x_max=1.0_dp)
        
        ! Check number of nodes
        if (mesh%n_nodes == 5) then
            print *, "PASS: Correct number of nodes"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Wrong number of nodes"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check number of elements
        if (mesh%n_elements == 4) then
            print *, "PASS: Correct number of elements"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Wrong number of elements"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check node positions
        do i = 1, mesh%n_nodes
            if (abs(mesh%nodes(i) - (i-1)*0.25_dp) < tol) then
                n_tests_passed = n_tests_passed + 1
            else
                print *, "FAIL: Wrong position for node", i
                n_tests_failed = n_tests_failed + 1
            end if
        end do
        
    end subroutine test_uniform_mesh
    
    subroutine test_connectivity()
        integer :: i
        
        ! Create simple mesh
        call mesh%create_uniform(n_nodes=3, x_min=0.0_dp, x_max=1.0_dp)
        
        ! Check connectivity: element 1 should connect nodes 1-2
        if (mesh%connectivity(1,1) == 1 .and. mesh%connectivity(2,1) == 2) then
            print *, "PASS: Element 1 connectivity"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Element 1 connectivity"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check connectivity: element 2 should connect nodes 2-3
        if (mesh%connectivity(1,2) == 2 .and. mesh%connectivity(2,2) == 3) then
            print *, "PASS: Element 2 connectivity"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Element 2 connectivity"
            n_tests_failed = n_tests_failed + 1
        end if
        
    end subroutine test_connectivity
    
    subroutine test_element_sizes()
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: h
        integer :: i
        
        ! Create uniform mesh
        call mesh%create_uniform(n_nodes=5, x_min=0.0_dp, x_max=2.0_dp)
        
        ! All elements should have size 0.5
        do i = 1, mesh%n_elements
            h = mesh%element_size(i)
            if (abs(h - 0.5_dp) < tol) then
                n_tests_passed = n_tests_passed + 1
            else
                print *, "FAIL: Wrong size for element", i
                n_tests_failed = n_tests_failed + 1
            end if
        end do
        
    end subroutine test_element_sizes

end program test_mesh_1d
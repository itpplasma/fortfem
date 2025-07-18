program test_mesh_2d
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test 1: Create rectangular mesh
    call test_rectangular_mesh()
    
    ! Test 2: Test mesh connectivity
    call test_mesh_connectivity()
    
    ! Test 3: Test boundary detection
    call test_boundary_detection()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_rectangular_mesh()
        type(mesh_2d_t) :: mesh
        real(dp), parameter :: tol = 1.0e-14_dp
        
        ! Create 2x2 rectangular mesh on [0,1]x[0,1]
        call mesh%create_rectangular(nx=3, ny=3, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        ! Check dimensions
        if (mesh%n_vertices == 9) then
            print *, "PASS: Correct number of vertices"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Wrong number of vertices"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! 2x2 grid has 8 triangles (2 per square)
        if (mesh%n_triangles == 8) then
            print *, "PASS: Correct number of triangles"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Wrong number of triangles"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check corner vertices
        if (abs(mesh%vertices(1,1) - 0.0_dp) < tol .and. &
            abs(mesh%vertices(2,1) - 0.0_dp) < tol) then
            print *, "PASS: Bottom-left corner correct"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Bottom-left corner incorrect"
            n_tests_failed = n_tests_failed + 1
        end if
        
        if (abs(mesh%vertices(1,9) - 1.0_dp) < tol .and. &
            abs(mesh%vertices(2,9) - 1.0_dp) < tol) then
            print *, "PASS: Top-right corner correct"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Top-right corner incorrect"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call mesh%destroy()
        
    end subroutine test_rectangular_mesh
    
    subroutine test_mesh_connectivity()
        type(mesh_2d_t) :: mesh
        integer :: t, v, count
        
        ! Create simple 2x2 mesh
        call mesh%create_rectangular(nx=3, ny=3, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        ! Build connectivity
        call mesh%build_connectivity()
        
        ! Check that connectivity was built
        if (allocated(mesh%vertex_to_triangles)) then
            print *, "PASS: Vertex-to-triangle connectivity allocated"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Vertex-to-triangle connectivity not allocated"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check center vertex (should be in 6 triangles)
        count = 0
        do t = 1, mesh%n_triangles
            do v = 1, 3
                if (mesh%triangles(v,t) == 5) count = count + 1
            end do
        end do
        
        if (count == 6) then
            print *, "PASS: Center vertex in correct number of triangles"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Center vertex triangle count incorrect"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call mesh%destroy()
        
    end subroutine test_mesh_connectivity
    
    subroutine test_boundary_detection()
        type(mesh_2d_t) :: mesh
        
        ! Create simple mesh
        call mesh%create_rectangular(nx=3, ny=3, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        call mesh%build_connectivity()
        call mesh%find_boundary()
        
        ! Check boundary edges
        if (mesh%n_boundary_edges == 8) then
            print *, "PASS: Correct number of boundary edges"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Wrong number of boundary edges"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check boundary vertices (should be 8 for square)
        if (count(mesh%is_boundary_vertex) == 8) then
            print *, "PASS: Correct number of boundary vertices"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Wrong number of boundary vertices"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call mesh%destroy()
        
    end subroutine test_boundary_detection

end program test_mesh_2d
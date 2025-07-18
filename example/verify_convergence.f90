program verify_convergence
    use fortfem
    implicit none
    
    type(poisson_2d_t) :: solver
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: solution(:), exact(:), error(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, L2_error, h, rate
    integer :: i, k, n, nx_values(5)
    real(dp) :: errors(5), h_values(5)
    
    print *, "FEM Convergence Rate Verification"
    print *, "================================="
    print *, ""
    print *, "Testing convergence for manufactured solution:"
    print *, "u = sin(π*x) * sin(π*y)"
    print *, "Expected convergence rate: O(h²) for P1 elements"
    print *, ""
    
    ! Test different mesh sizes
    nx_values = [5, 10, 20, 40, 80]
    
    do n = 1, size(nx_values)
        ! Create mesh
        call mesh%create_rectangular(nx=nx_values(n), ny=nx_values(n), &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        h = 1.0_dp / real(nx_values(n), dp)
        h_values(n) = h
        
        print '(a,i0,a,i0,a,f6.4)', "Mesh ", nx_values(n), "×", nx_values(n), &
              " (h=", h, ")"
        
        ! Initialize solver
        call solver%init("lapack")
        call solver%set_mesh(mesh)
        
        ! Set boundary conditions
        allocate(boundary_nodes(2*nx_values(n) + 2*nx_values(n) - 4))
        allocate(boundary_values(size(boundary_nodes)))
        
        k = 0
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
                abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
                k = k + 1
                boundary_nodes(k) = i
                boundary_values(k) = 0.0_dp
            end if
        end do
        
        call solver%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
        
        ! Solve
        call solver%solve(test_source_function)
        solution = solver%get_solution()
        
        ! Compute error
        allocate(exact(mesh%n_vertices))
        allocate(error(mesh%n_vertices))
        
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            exact(i) = exact_solution(x, y)
            error(i) = abs(solution(i) - exact(i))
        end do
        
        L2_error = sqrt(sum(error**2) / mesh%n_vertices)
        errors(n) = L2_error
        
        print '(a,i6,a,es10.3)', "  Vertices: ", mesh%n_vertices, &
              ", L2 error: ", L2_error
        
        if (n > 1) then
            rate = log(errors(n-1) / errors(n)) / log(h_values(n-1) / h_values(n))
            print '(a,f6.2)', "  Convergence rate: ", rate
        end if
        
        ! Clean up
        call solver%destroy()
        call mesh%destroy()
        deallocate(solution, exact, error)
        deallocate(boundary_nodes, boundary_values)
        
        print *, ""
    end do
    
    print *, "Summary of convergence rates:"
    print *, "h         L2 error     Rate"
    print *, "-----------------------------"
    do n = 1, size(nx_values)
        if (n > 1) then
            rate = log(errors(n-1) / errors(n)) / log(h_values(n-1) / h_values(n))
            print '(f6.4,2x,es10.3,2x,f6.2)', h_values(n), errors(n), rate
        else
            print '(f6.4,2x,es10.3,2x,a)', h_values(n), errors(n), "  --"
        end if
    end do
    
    print *, ""
    print *, "Expected rate for P1 elements: ~2.0"
    print *, "Verification complete!"
    
contains

    pure function test_source_function(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function test_source_function
    
    pure function exact_solution(x, y) result(u)
        real(dp), intent(in) :: x, y
        real(dp) :: u
        u = sin(pi * x) * sin(pi * y)
    end function exact_solution

end program verify_convergence
program verify_simple_problem
    use fortfem
    implicit none
    
    type(poisson_2d_t) :: solver
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: solution(:), exact(:), error(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, L2_error, max_error
    integer :: i, k
    
    print *, "Simple Problem Verification"
    print *, "==========================="
    print *, ""
    print *, "Testing with simple quadratic solution:"
    print *, "u = x*(1-x)*y*(1-y)"
    print *, "f = -2*x*(1-x) - 2*y*(1-y)"
    print *, "This satisfies u = 0 on boundary exactly"
    print *, ""
    
    ! Create small mesh first
    call mesh%create_rectangular(nx=10, ny=10, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Initialize solver
    call solver%init("lapack")
    call solver%set_mesh(mesh)
    
    ! Set boundary conditions
    allocate(boundary_nodes(4*10 - 4))
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
    print '(a,i0)', "Boundary nodes: ", k
    
    ! Solve
    call solver%solve(simple_source_function)
    solution = solver%get_solution()
    
    ! Compute error
    allocate(exact(mesh%n_vertices))
    allocate(error(mesh%n_vertices))
    
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        exact(i) = simple_exact_solution(x, y)
        error(i) = abs(solution(i) - exact(i))
    end do
    
    L2_error = sqrt(sum(error**2) / mesh%n_vertices)
    max_error = maxval(error)
    
    print *, ""
    print '(a,es12.5)', "L2 error: ", L2_error
    print '(a,es12.5)', "Max error: ", max_error
    
    ! Check some specific points
    print *, ""
    print *, "Solution at specific points:"
    print *, "Point      x       y     Computed   Exact     Error"
    print *, "-----------------------------------------------------"
    
    do i = 1, min(10, mesh%n_vertices)
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        print '(i5,2f8.3,3f10.6)', i, x, y, solution(i), exact(i), error(i)
    end do
    
    ! Test element matrices directly
    print *, ""
    print *, "Testing element matrices:"
    call test_element_matrices()
    
    ! Clean up
    call solver%destroy()
    call mesh%destroy()
    deallocate(solution, exact, error)
    deallocate(boundary_nodes, boundary_values)
    
    print *, ""
    print *, "Verification complete!"
    
contains

    pure function simple_source_function(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        ! f = -div(grad(u)) where u = x*(1-x)*y*(1-y)
        ! d²u/dx² = -2*y*(1-y), d²u/dy² = -2*x*(1-x)
        f = 2.0_dp * y * (1.0_dp - y) + 2.0_dp * x * (1.0_dp - x)
    end function simple_source_function
    
    pure function simple_exact_solution(x, y) result(u)
        real(dp), intent(in) :: x, y
        real(dp) :: u
        u = x * (1.0_dp - x) * y * (1.0_dp - y)
    end function simple_exact_solution
    
    subroutine test_element_matrices()
        type(assembly_2d_t) :: assembly
        real(dp) :: vertices(2,3), mass(3,3), stiff(3,3), load(3)
        integer :: i, j
        
        print *, "  Testing unit right triangle:"
        
        ! Unit right triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        ! Mass matrix
        call assembly%element_mass_matrix(vertices, mass)
        print '(a,f10.6)', "    Mass matrix sum: ", sum(mass)
        print '(a,f10.6)', "    Expected (area): ", 0.5_dp
        
        ! Stiffness matrix
        call assembly%element_stiffness_matrix(vertices, stiff)
        print '(a,f10.6)', "    Stiffness diagonal sum: ", stiff(1,1) + stiff(2,2) + stiff(3,3)
        print '(a,f10.6)', "    Max row sum: ", maxval(abs([sum(stiff(1,:)), sum(stiff(2,:)), sum(stiff(3,:))]))
        
        ! Load vector
        call assembly%element_load_vector(vertices, unit_source, load)
        print '(a,f10.6)', "    Load vector sum: ", sum(load)
        print '(a,f10.6)', "    Expected (area): ", 0.5_dp
        
    end subroutine test_element_matrices
    
    pure function unit_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 1.0_dp
    end function unit_source

end program verify_simple_problem
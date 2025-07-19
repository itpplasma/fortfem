program test_piola_debug
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_l2_projection
    use fortfem_basis_edge_2d_interface
    implicit none

    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    real(dp), allocatable :: coeff(:)
    integer :: n_dofs, i, j
    real(dp) :: xi, eta, triangle_area
    real(dp) :: basis_values(2, 3), E_exact(2), E_computed(2)
    integer :: triangle_dofs(3)
    
    ! Create simple 2x2 mesh
    call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()
    
    call space%init(mesh)
    n_dofs = space%get_n_dofs()
    allocate(coeff(n_dofs))
    
    print *, "=== Piola Transformation Debug ==="
    print *, "Triangles:", mesh%n_triangles
    print *, "Edges:", mesh%n_edges  
    print *, "DOFs:", n_dofs
    print *, ""
    
    ! Project constant field E = [1, 0]
    call project_edge_dof_analytical(mesh, space, constant_field_x, coeff)
    
    print *, "Projected coefficients for E = [1, 0]:"
    do i = 1, n_dofs
        print '("DOF", I2, ": ", F8.4)', i, coeff(i)
    end do
    print *, ""
    
    ! Test at triangle centers
    xi = 1.0_dp/3.0_dp
    eta = 1.0_dp/3.0_dp
    
    print *, "=== Basis Function Values at Triangle Centers ==="
    do i = 1, min(4, mesh%n_triangles)
        call space%get_triangle_dofs(i, triangle_dofs)
        triangle_area = compute_triangle_area(mesh, i)
        
        print '("Triangle", I1, " (area=", F6.3, "):")', i, triangle_area
        print '("  DOFs: ", 3I3)', triangle_dofs + 1  ! Convert to 1-based
        
        ! Evaluate basis functions (can't access internal Piola method)
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
        
        print *, "  Basis values at center:"
        do j = 1, 3
            print '("    φ", I1, " = [", F7.4, ", ", F7.4, "]")', j, basis_values(:, j)
        end do
        
        ! Evaluate field at this point
        call space%evaluate_at_point(i, xi, eta, coeff, E_computed)
        call constant_field_x(0.5_dp, 0.5_dp, E_exact)
        
        print '("  E_computed = [", F7.4, ", ", F7.4, "]")', E_computed
        print '("  E_exact    = [", F7.4, ", ", F7.4, "]")', E_exact
        print '("  Error      = ", F7.4)', sqrt(sum((E_computed - E_exact)**2))
        print *, ""
    end do
    
    deallocate(coeff)
    call space%destroy()
    call mesh%destroy()

contains

    subroutine constant_field_x(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        E(1) = 1.0_dp
        E(2) = 0.0_dp
    end subroutine

    function compute_triangle_area(mesh, triangle_idx) result(area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp) :: area
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        area = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
    end function

end program test_piola_debug
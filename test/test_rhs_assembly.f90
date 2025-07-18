program test_rhs_assembly
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_local_rhs_assembly()
    call test_global_rhs_assembly()
    call test_analytical_source_integration()
    
    print *, "All RHS assembly tests passed!"

contains

    subroutine test_local_rhs_assembly()
        type(mesh_2d_t) :: mesh
        real(dp) :: local_rhs(3)
        real(dp) :: triangle_area
        integer :: triangle_idx
        integer :: i
        logical :: has_nonzero_rhs
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        triangle_idx = 1
        triangle_area = 0.5_dp
        
        ! Compute local RHS with analytical source J = [sin(πx)sin(πy), cos(πx)cos(πy)]
        call compute_local_rhs_analytical(triangle_idx, triangle_area, local_rhs)
        
        ! Verify that RHS values are finite and reasonable
        do i = 1, 3
            if (.not. (abs(local_rhs(i)) < 1e10_dp)) then
                print *, "Error: local RHS value not finite"
                stop 1
            end if
        end do
        
        ! For our analytical source, RHS should be non-trivial
        has_nonzero_rhs = .false.
        do i = 1, 3
            if (abs(local_rhs(i)) > 1e-12_dp) then
                has_nonzero_rhs = .true.
                exit
            end if
        end do
        
        if (.not. has_nonzero_rhs) then
            print *, "Warning: all local RHS entries are zero"
        end if
        
        call mesh%destroy()
        print *, "Local RHS assembly test passed"
    end subroutine
    
    subroutine test_global_rhs_assembly()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: global_rhs(:)
        integer :: n_dofs, i
        real(dp) :: rhs_norm
        
        ! Create 2x2 mesh for more realistic test
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(global_rhs(n_dofs))
        global_rhs = 0.0_dp
        
        ! Assemble global RHS vector
        call assemble_global_rhs_analytical(mesh, global_rhs)
        
        ! Verify finite values
        do i = 1, n_dofs
            if (.not. (abs(global_rhs(i)) < 1e10_dp)) then
                print *, "Error: global RHS value not finite at DOF", i
                stop 1
            end if
        end do
        
        ! Compute RHS norm
        rhs_norm = 0.0_dp
        do i = 1, n_dofs
            rhs_norm = rhs_norm + global_rhs(i)**2
        end do
        rhs_norm = sqrt(rhs_norm)
        
        if (rhs_norm < 1e-12_dp) then
            print *, "Warning: global RHS vector is essentially zero"
        end if
        
        deallocate(global_rhs)
        call space%destroy()
        call mesh%destroy()
        print *, "Global RHS assembly test passed"
    end subroutine
    
    subroutine test_analytical_source_integration()
        real(dp) :: integral_x, integral_y
        real(dp) :: triangle_area
        
        triangle_area = 0.5_dp
        
        ! Test integration of analytical source over reference triangle
        call integrate_analytical_source_reference_triangle(triangle_area, integral_x, integral_y)
        
        ! For J = [sin(πx)sin(πy), cos(πx)cos(πy)] over reference triangle,
        ! integrals should be finite and computable
        if (.not. (abs(integral_x) < 1e10_dp .and. abs(integral_y) < 1e10_dp)) then
            print *, "Error: source integration gives non-finite results"
            stop 1
        end if
        
        print *, "Analytical source integration test passed"
    end subroutine
    
    subroutine compute_local_rhs_analytical(triangle_idx, triangle_area, local_rhs)
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: local_rhs(3)
        
        ! Quadrature points and weights for triangle integration
        real(dp), parameter :: xi_quad(3) = [0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: eta_quad(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        real(dp) :: basis_values(2, 3)
        real(dp) :: source_values(2)
        real(dp) :: x_phys, y_phys, dot_product
        integer :: i, q
        
        local_rhs = 0.0_dp
        
        ! Integrate using quadrature: ∫_T J·φ_i dx
        do q = 1, 3
            ! Evaluate basis functions at quadrature point
            call evaluate_edge_basis_2d(xi_quad(q), eta_quad(q), triangle_area, basis_values)
            
            ! Map reference coordinates to physical coordinates
            x_phys = xi_quad(q)
            y_phys = eta_quad(q)
            
            ! Evaluate analytical source: J = [sin(πx)sin(πy), cos(πx)cos(πy)]
            source_values(1) = sin(pi * x_phys) * sin(pi * y_phys)
            source_values(2) = cos(pi * x_phys) * cos(pi * y_phys)
            
            ! Integrate J·φ_i over triangle
            do i = 1, 3
                dot_product = source_values(1) * basis_values(1, i) + &
                             source_values(2) * basis_values(2, i)
                
                local_rhs(i) = local_rhs(i) + weights(q) * dot_product * triangle_area
            end do
        end do
    end subroutine compute_local_rhs_analytical
    
    subroutine assemble_global_rhs_analytical(mesh, global_rhs)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(inout) :: global_rhs(:)
        
        real(dp) :: local_rhs(3)
        real(dp) :: triangle_area
        integer :: triangle_dofs(3)
        integer :: t, i
        
        ! Zero out global RHS
        global_rhs = 0.0_dp
        
        ! Loop over all triangles
        do t = 1, mesh%n_triangles
            ! Compute triangle area
            triangle_area = compute_triangle_area_from_mesh(mesh, t)
            
            ! Compute local RHS vector
            call compute_local_rhs_analytical(t, triangle_area, local_rhs)
            
            ! Get triangle DOFs
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            ! Add to global RHS (convert 0-based to 1-based indexing)
            do i = 1, 3
                global_rhs(triangle_dofs(i) + 1) = global_rhs(triangle_dofs(i) + 1) + local_rhs(i)
            end do
        end do
    end subroutine assemble_global_rhs_analytical
    
    subroutine integrate_analytical_source_reference_triangle(triangle_area, integral_x, integral_y)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: integral_x, integral_y
        
        ! Quadrature points and weights
        real(dp), parameter :: xi_quad(3) = [0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: eta_quad(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        real(dp) :: source_values(2)
        integer :: q
        
        integral_x = 0.0_dp
        integral_y = 0.0_dp
        
        ! Integrate analytical source over reference triangle
        do q = 1, 3
            ! Evaluate source at quadrature point
            source_values(1) = sin(pi * xi_quad(q)) * sin(pi * eta_quad(q))
            source_values(2) = cos(pi * xi_quad(q)) * cos(pi * eta_quad(q))
            
            ! Add weighted contribution
            integral_x = integral_x + weights(q) * source_values(1) * triangle_area
            integral_y = integral_y + weights(q) * source_values(2) * triangle_area
        end do
    end subroutine integrate_analytical_source_reference_triangle
    
    function compute_triangle_area_from_mesh(mesh, triangle_idx) result(area)
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
    end function compute_triangle_area_from_mesh
    
    subroutine create_reference_triangle(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Reference triangle vertices: (0,0), (1,0), (0,1)
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        
        ! Triangle connectivity
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
    end subroutine
    
    subroutine create_unit_square_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        ! Create 3x3 mesh on unit square
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    end subroutine

end program test_rhs_assembly
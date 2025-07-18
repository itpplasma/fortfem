program test_vector_source_projection
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_vector_source_evaluation()
    call test_source_projection_onto_edges()
    call test_rhs_vector_assembly()
    
    print *, "All vector source projection tests passed!"

contains

    subroutine test_vector_source_evaluation()
        real(dp) :: x, y, source_values(2)
        
        ! Test analytical vector source: J = [sin(πx)sin(πy), cos(πx)cos(πy)]
        x = 0.5_dp
        y = 0.5_dp
        
        call evaluate_analytical_source(x, y, source_values)
        
        ! At (0.5, 0.5): sin(π/2) = 1, cos(π/2) = 0
        if (abs(source_values(1) - 1.0_dp) > 1e-12_dp) then
            print *, "Error: incorrect source evaluation at (0.5, 0.5)"
            print *, "Expected [1.0, 0.0], got [", source_values(1), ",", source_values(2), "]"
            stop 1
        end if
        
        if (abs(source_values(2)) > 1e-12_dp) then
            print *, "Error: incorrect source evaluation at (0.5, 0.5)"
            print *, "Expected [1.0, 0.0], got [", source_values(1), ",", source_values(2), "]"
            stop 1
        end if
        
        ! Test at origin
        x = 0.0_dp
        y = 0.0_dp
        call evaluate_analytical_source(x, y, source_values)
        
        ! At (0.0, 0.0): sin(0) = 0, cos(0) = 1
        if (abs(source_values(1)) > 1e-12_dp .or. abs(source_values(2) - 1.0_dp) > 1e-12_dp) then
            print *, "Error: incorrect source evaluation at origin"
            print *, "Expected [0.0, 1.0], got [", source_values(1), ",", source_values(2), "]"
            stop 1
        end if
        
        print *, "Vector source evaluation test passed"
    end subroutine
    
    subroutine test_source_projection_onto_edges()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp) :: local_rhs(3)
        real(dp) :: triangle_area
        integer :: triangle_idx
        integer :: i
        logical :: has_nonzero
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        triangle_idx = 1
        triangle_area = 0.5_dp
        
        ! Compute local RHS vector: ∫_T J·φ_i dx
        call compute_local_source_projection(triangle_idx, triangle_area, local_rhs)
        
        ! Verify that RHS values are reasonable
        do i = 1, 3
            if (abs(local_rhs(i)) > 1e10_dp) then
                print *, "Error: RHS value too large - numerical issues"
                stop 1
            end if
        end do
        
        ! For our analytical source on reference triangle, 
        ! RHS should be non-zero due to source term
        has_nonzero = .false.
        do i = 1, 3
            if (abs(local_rhs(i)) > 1e-12_dp) then
                has_nonzero = .true.
                exit
            end if
        end do
        
        if (.not. has_nonzero) then
            print *, "Warning: all RHS entries are zero - check source function"
        end if
        
        call space%destroy()
        call mesh%destroy()
        print *, "Source projection onto edges test passed"
    end subroutine
    
    subroutine test_rhs_vector_assembly()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: global_rhs(:)
        integer :: n_dofs
        integer :: i
        
        ! Create 2x2 mesh
        call create_simple_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(global_rhs(n_dofs))
        global_rhs = 0.0_dp
        
        ! Assemble global RHS vector from source term
        call assemble_global_rhs_vector(mesh, global_rhs)
        
        ! Verify assembly produced reasonable values
        do i = 1, n_dofs
            if (abs(global_rhs(i)) > 1e10_dp) then
                print *, "Error: global RHS value too large - numerical issues"
                stop 1
            end if
        end do
        
        deallocate(global_rhs)
        call space%destroy()
        call mesh%destroy()
        print *, "RHS vector assembly test passed"
    end subroutine
    
    subroutine evaluate_analytical_source(x, y, source_values)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: source_values(2)
        
        ! Analytical vector source: J = [sin(πx)sin(πy), cos(πx)cos(πy)]
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        source_values(1) = sin(pi * x) * sin(pi * y)
        source_values(2) = cos(pi * x) * cos(pi * y)
    end subroutine evaluate_analytical_source
    
    subroutine compute_local_source_projection(triangle_idx, triangle_area, local_rhs)
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: local_rhs(3)
        
        ! Quadrature points and weights for triangle integration
        real(dp), parameter :: xi_quad(3) = [0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: eta_quad(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
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
            ! For reference triangle: (x,y) = (xi, eta)
            x_phys = xi_quad(q)
            y_phys = eta_quad(q)
            
            ! Evaluate source at physical point
            call evaluate_analytical_source(x_phys, y_phys, source_values)
            
            ! Integrate J·φ_i over triangle
            do i = 1, 3
                dot_product = source_values(1) * basis_values(1, i) + &
                             source_values(2) * basis_values(2, i)
                
                local_rhs(i) = local_rhs(i) + weights(q) * dot_product * triangle_area
            end do
        end do
    end subroutine compute_local_source_projection
    
    subroutine assemble_global_rhs_vector(mesh, global_rhs)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(inout) :: global_rhs(:)
        
        real(dp) :: local_rhs(3)
        real(dp) :: triangle_area
        integer :: triangle_dofs(3)
        integer :: t, i
        
        ! Loop over all triangles
        do t = 1, mesh%n_triangles
            ! Compute triangle area
            triangle_area = compute_triangle_area_from_mesh(mesh, t)
            
            ! Compute local RHS vector
            call compute_local_source_projection(t, triangle_area, local_rhs)
            
            ! Get triangle DOFs
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            ! Add to global RHS (convert 0-based to 1-based indexing)
            do i = 1, 3
                global_rhs(triangle_dofs(i) + 1) = global_rhs(triangle_dofs(i) + 1) + local_rhs(i)
            end do
        end do
    end subroutine assemble_global_rhs_vector
    
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
    
    subroutine create_simple_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        ! Create 2x2 mesh
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    end subroutine

end program test_vector_source_projection
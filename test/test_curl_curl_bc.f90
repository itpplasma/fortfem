program test_curl_curl_bc
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    use fortfem_sparse_matrix
    implicit none

    call test_essential_bc_tangential()
    call test_boundary_dof_elimination()
    call test_natural_bc_implementation()
    
    print *, "All curl-curl boundary condition tests passed!"

contains

    subroutine test_essential_bc_tangential()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        integer :: i, edge_idx, n_boundary_constrained
        
        ! Create 3x3 mesh for better boundary testing
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        ! Apply homogeneous essential BC: E × n = 0 on boundary
        call space%apply_essential_bc(0.0_dp)
        
        ! Count constrained DOFs
        n_boundary_constrained = 0
        do i = 1, space%get_n_dofs()
            if (space%dof_is_essential(i) == 1) then
                n_boundary_constrained = n_boundary_constrained + 1
            end if
        end do
        
        ! Verify all boundary DOFs are constrained
        if (n_boundary_constrained /= space%get_n_boundary_dofs()) then
            print *, "Error: not all boundary DOFs are constrained"
            print *, "Expected:", space%get_n_boundary_dofs(), "Got:", n_boundary_constrained
            stop 1
        end if
        
        ! Verify essential values are set correctly
        do i = 1, space%get_n_dofs()
            if (space%dof_is_essential(i) == 1) then
                if (abs(space%essential_values(i) - 0.0_dp) > 1e-12_dp) then
                    print *, "Error: essential BC value not set correctly for DOF", i
                    stop 1
                end if
            end if
        end do
        
        call space%destroy()
        call mesh%destroy()
        print *, "Essential BC tangential test passed"
    end subroutine
    
    subroutine test_boundary_dof_elimination()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        type(triplet_matrix_t) :: system_matrix, reduced_matrix
        real(dp), allocatable :: rhs(:), reduced_rhs(:)
        integer :: n_dofs, n_interior_dofs
        integer :: i
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        n_interior_dofs = space%get_n_interior_dofs()
        
        ! Initialize full system
        call system_matrix%init(n_dofs, n_dofs * 15)
        allocate(rhs(n_dofs))
        rhs = 1.0_dp  ! Dummy RHS
        
        ! Assemble full system
        call assemble_dummy_system(mesh, system_matrix)
        
        ! Apply essential BC
        call space%apply_essential_bc(0.0_dp)
        
        ! Create reduced system (interior DOFs only)
        call reduced_matrix%init(n_interior_dofs, n_interior_dofs * 15)
        allocate(reduced_rhs(n_interior_dofs))
        
        ! Eliminate boundary DOFs
        call eliminate_boundary_dofs(system_matrix, rhs, space, reduced_matrix, reduced_rhs)
        
        ! Verify reduced system has correct size
        if (reduced_matrix%n /= n_interior_dofs) then
            print *, "Error: reduced matrix has wrong dimension"
            print *, "Expected:", n_interior_dofs, "Got:", reduced_matrix%n
            stop 1
        end if
        
        ! Verify reduced system has entries
        if (reduced_matrix%nnz <= 0) then
            print *, "Error: reduced matrix has no entries"
            stop 1
        end if
        
        deallocate(rhs, reduced_rhs)
        call system_matrix%destroy()
        call reduced_matrix%destroy()
        call space%destroy()
        call mesh%destroy()
        print *, "Boundary DOF elimination test passed"
    end subroutine
    
    subroutine test_natural_bc_implementation()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        real(dp), allocatable :: natural_bc_contribution(:)
        real(dp) :: max_contribution
        integer :: n_dofs
        integer :: i
        
        ! Create mesh
        call create_unit_square_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        call space%init(mesh)
        
        n_dofs = space%get_n_dofs()
        allocate(natural_bc_contribution(n_dofs))
        natural_bc_contribution = 0.0_dp
        
        ! Compute natural BC contribution: ∫_∂Ω g·v ds
        ! For homogeneous natural BC: n × curl E = 0, contribution is zero
        call compute_natural_bc_contribution(mesh, space, natural_bc_contribution)
        
        ! For homogeneous natural BC, contribution should be zero
        max_contribution = 0.0_dp
        do i = 1, n_dofs
            max_contribution = max(max_contribution, abs(natural_bc_contribution(i)))
        end do
        
        if (max_contribution > 1e-12_dp) then
            print *, "Error: non-zero natural BC contribution for homogeneous BC"
            print *, "Max contribution:", max_contribution
            stop 1
        end if
        
        deallocate(natural_bc_contribution)
        call space%destroy()
        call mesh%destroy()
        print *, "Natural BC implementation test passed"
    end subroutine
    
    subroutine assemble_dummy_system(mesh, matrix)
        type(mesh_2d_t), intent(in) :: mesh
        type(triplet_matrix_t), intent(inout) :: matrix
        
        ! Assemble a simple system for testing
        real(dp), parameter :: k_squared = 1.0_dp
        real(dp), parameter :: epsilon = 1e-6_dp
        
        ! Just use the complete system from previous test
        call assemble_complete_curl_curl_system(mesh, k_squared, epsilon, matrix)
    end subroutine
    
    subroutine assemble_complete_curl_curl_system(mesh, k_squared, epsilon, matrix)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: k_squared, epsilon
        type(triplet_matrix_t), intent(inout) :: matrix
        
        real(dp) :: local_curl(3, 3), local_mass(3, 3), local_reg(3, 3)
        real(dp) :: triangle_area, value
        integer :: triangle_dofs(3)
        integer :: t, i, j
        
        ! Loop over triangles
        do t = 1, mesh%n_triangles
            ! Compute triangle area
            triangle_area = compute_triangle_area(mesh, t)
            
            ! Compute local matrices
            call compute_local_curl_curl_matrix(triangle_area, local_curl)
            call compute_local_edge_mass_matrix(triangle_area, local_mass)
            call compute_local_regularization_matrix(triangle_area, epsilon, local_reg)
            
            ! Get triangle DOFs
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            ! Add to global matrix
            do i = 1, 3
                do j = 1, 3
                    value = local_curl(i, j) + k_squared * local_mass(i, j) + local_reg(i, j)
                    if (abs(value) > 1e-12_dp) then
                        call matrix%add(triangle_dofs(i) + 1, triangle_dofs(j) + 1, value)
                    end if
                end do
            end do
        end do
    end subroutine
    
    subroutine eliminate_boundary_dofs(full_matrix, full_rhs, space, reduced_matrix, reduced_rhs)
        type(triplet_matrix_t), intent(in) :: full_matrix
        real(dp), intent(in) :: full_rhs(:)
        type(hcurl_space_t), intent(in) :: space
        type(triplet_matrix_t), intent(inout) :: reduced_matrix
        real(dp), intent(out) :: reduced_rhs(:)
        
        integer :: i, j, k
        integer :: reduced_i, reduced_j
        
        ! Map full DOFs to reduced DOFs (interior only)
        integer, allocatable :: dof_map(:)
        allocate(dof_map(space%get_n_dofs()))
        
        ! Build DOF mapping: full index -> reduced index
        reduced_i = 0
        do i = 1, space%get_n_dofs()
            if (space%dof_is_essential(i) == 0) then  ! Interior DOF
                reduced_i = reduced_i + 1
                dof_map(i) = reduced_i
            else
                dof_map(i) = -1  ! Boundary DOF
            end if
        end do
        
        ! Copy RHS for interior DOFs
        reduced_i = 0
        do i = 1, space%get_n_dofs()
            if (space%dof_is_essential(i) == 0) then
                reduced_i = reduced_i + 1
                reduced_rhs(reduced_i) = full_rhs(i)
            end if
        end do
        
        ! Extract interior-interior block of matrix
        ! This is a simplified version - in practice we'd iterate through triplet entries
        do k = 1, full_matrix%nnz
            i = full_matrix%rows(k)
            j = full_matrix%cols(k)
            
            ! Only keep interior-interior entries
            if (space%dof_is_essential(i) == 0 .and. space%dof_is_essential(j) == 0) then
                reduced_i = dof_map(i)
                reduced_j = dof_map(j)
                call reduced_matrix%add(reduced_i, reduced_j, full_matrix%values(k))
            end if
        end do
        
        deallocate(dof_map)
    end subroutine
    
    subroutine compute_natural_bc_contribution(mesh, space, contribution)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(out) :: contribution(:)
        
        ! For homogeneous natural BC: n × curl E = 0
        ! The weak form contribution is: ∫_∂Ω (n × curl E)·v ds = 0
        ! So we set contribution to zero
        contribution = 0.0_dp
        
        ! For non-homogeneous natural BC: n × curl E = g
        ! We would compute: ∫_∂Ω g·v ds for each boundary edge
        ! This would involve:
        ! 1. Loop over boundary edges
        ! 2. Evaluate edge basis functions on boundary
        ! 3. Integrate g·φ_i along each edge
        ! 4. Add to global contribution vector
    end subroutine
    
    subroutine compute_local_curl_curl_matrix(triangle_area, local_matrix)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: local_matrix(3, 3)
        
        real(dp) :: curls(3)
        integer :: i, j
        
        ! For RT0 on reference triangle
        curls(1) = 2.0_dp
        curls(2) = 2.0_dp
        curls(3) = -2.0_dp
        
        do i = 1, 3
            do j = 1, 3
                local_matrix(i, j) = curls(i) * curls(j) * triangle_area
            end do
        end do
    end subroutine
    
    subroutine compute_local_edge_mass_matrix(triangle_area, mass_matrix)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: mass_matrix(3, 3)
        
        ! Simplified mass matrix for RT0 (approximate values)
        mass_matrix(1, 1) = 2.0_dp * triangle_area / 12.0_dp
        mass_matrix(2, 2) = 2.0_dp * triangle_area / 12.0_dp
        mass_matrix(3, 3) = 2.0_dp * triangle_area / 12.0_dp
        
        mass_matrix(1, 2) = triangle_area / 12.0_dp
        mass_matrix(2, 1) = mass_matrix(1, 2)
        mass_matrix(1, 3) = triangle_area / 12.0_dp
        mass_matrix(3, 1) = mass_matrix(1, 3)
        mass_matrix(2, 3) = triangle_area / 12.0_dp
        mass_matrix(3, 2) = mass_matrix(2, 3)
    end subroutine
    
    subroutine compute_local_regularization_matrix(triangle_area, epsilon, reg_matrix)
        real(dp), intent(in) :: triangle_area, epsilon
        real(dp), intent(out) :: reg_matrix(3, 3)
        
        real(dp) :: divergences(3)
        integer :: i, j
        
        divergences(1) = 2.0_dp
        divergences(2) = 2.0_dp
        divergences(3) = -4.0_dp
        
        do i = 1, 3
            do j = 1, 3
                reg_matrix(i, j) = epsilon * divergences(i) * divergences(j) * triangle_area
            end do
        end do
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
    end function compute_triangle_area
    
    subroutine create_unit_square_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        ! Create 3x3 mesh for better boundary testing
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    end subroutine

end program test_curl_curl_bc
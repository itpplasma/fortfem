program test_global_curl_curl_assembly
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    use fortfem_mesh_2d
    use fortfem_sparse_matrix
    implicit none

    call test_global_assembly_two_triangles()
    call test_global_assembly_consistency()
    call test_global_matrix_properties()
    
    print *, "All global curl-curl assembly tests passed!"

contains

    subroutine test_global_assembly_two_triangles()
        type(mesh_2d_t) :: mesh
        type(triplet_matrix_t) :: global_matrix
        real(dp) :: local_matrix_1(3, 3), local_matrix_2(3, 3)
        real(dp) :: triangle_area_1, triangle_area_2
        real(dp) :: curls_1(3), curls_2(3)
        integer :: triangle_dofs_1(3), triangle_dofs_2(3)
        integer :: i, j, row, col
        real(dp) :: value
        
        ! Create two triangles mesh
        call create_two_triangles_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        ! Get triangle areas
        triangle_area_1 = compute_triangle_area(mesh, 1)
        triangle_area_2 = compute_triangle_area(mesh, 2)
        
        ! Get curl values for each triangle
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area_1, curls_1)
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area_2, curls_2)
        
        ! Get DOF mappings for each triangle
        call mesh%get_triangle_edge_dofs(1, triangle_dofs_1)
        call mesh%get_triangle_edge_dofs(2, triangle_dofs_2)
        
        ! Compute local matrices
        do i = 1, 3
            do j = 1, 3
                local_matrix_1(i, j) = curls_1(i) * curls_1(j) * triangle_area_1
                local_matrix_2(i, j) = curls_2(i) * curls_2(j) * triangle_area_2
            end do
        end do
        
        ! Initialize global matrix
        call global_matrix%init(mesh%get_n_edge_dofs(), 100)  ! 100 max entries for small test
        
        ! Assemble global matrix
        call assemble_local_matrix(global_matrix, local_matrix_1, triangle_dofs_1)
        call assemble_local_matrix(global_matrix, local_matrix_2, triangle_dofs_2)
        
        ! Verify some properties
        if (global_matrix%n /= mesh%get_n_edge_dofs()) then
            print *, "Error: global matrix has wrong dimension"
            stop 1
        end if
        
        call global_matrix%destroy()
        call mesh%destroy()
        print *, "Global assembly two triangles test passed"
    end subroutine
    
    subroutine test_global_assembly_consistency()
        type(mesh_2d_t) :: mesh
        type(triplet_matrix_t) :: global_matrix
        real(dp) :: local_matrix(3, 3)
        real(dp) :: triangle_area
        real(dp) :: curls(3)
        integer :: triangle_dofs(3)
        integer :: i, j
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        ! Get triangle area and curl values
        triangle_area = compute_triangle_area(mesh, 1)
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! Get DOF mapping
        call mesh%get_triangle_edge_dofs(1, triangle_dofs)
        
        ! Compute local matrix
        do i = 1, 3
            do j = 1, 3
                local_matrix(i, j) = curls(i) * curls(j) * triangle_area
            end do
        end do
        
        ! Initialize and assemble global matrix
        call global_matrix%init(mesh%get_n_edge_dofs(), 50)  ! 50 max entries
        call assemble_local_matrix(global_matrix, local_matrix, triangle_dofs)
        
        ! For single triangle, verify that entries were added correctly
        if (global_matrix%nnz /= 9) then  ! Should have 3x3 = 9 entries
            print *, "Error: expected 9 entries, got", global_matrix%nnz
            stop 1
        end if
        
        call global_matrix%destroy()
        call mesh%destroy()
        print *, "Global assembly consistency test passed"
    end subroutine
    
    subroutine test_global_matrix_properties()
        type(mesh_2d_t) :: mesh
        type(triplet_matrix_t) :: global_matrix
        real(dp) :: local_matrix(3, 3)
        real(dp) :: triangle_area
        real(dp) :: curls(3)
        integer :: triangle_dofs(3)
        integer :: i, j
        
        ! Create single triangle mesh
        call create_single_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        ! Get triangle area and curl values
        triangle_area = compute_triangle_area(mesh, 1)
        call evaluate_edge_basis_curl_2d(0.33_dp, 0.33_dp, triangle_area, curls)
        
        ! Get DOF mapping
        call mesh%get_triangle_edge_dofs(1, triangle_dofs)
        
        ! Compute local matrix
        do i = 1, 3
            do j = 1, 3
                local_matrix(i, j) = curls(i) * curls(j) * triangle_area
            end do
        end do
        
        ! Initialize and assemble global matrix
        call global_matrix%init(mesh%get_n_edge_dofs(), 50)  ! 50 max entries
        call assemble_local_matrix(global_matrix, local_matrix, triangle_dofs)
        
        ! Test that assembly completed successfully
        if (global_matrix%nnz == 0) then
            print *, "Error: no entries added to global matrix"
            stop 1
        end if
        
        call global_matrix%destroy()
        call mesh%destroy()
        print *, "Global matrix properties test passed"
    end subroutine
    
    subroutine assemble_local_matrix(global_matrix, local_matrix, dofs)
        type(triplet_matrix_t), intent(inout) :: global_matrix
        real(dp), intent(in) :: local_matrix(3, 3)
        integer, intent(in) :: dofs(3)
        integer :: i, j
        
        do i = 1, 3
            do j = 1, 3
                call global_matrix%add(dofs(i), dofs(j), local_matrix(i, j))
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
    end function
    
    subroutine create_single_triangle_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Triangle vertices
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
    
    subroutine create_two_triangles_mesh(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 4
        mesh%n_triangles = 2
        
        allocate(mesh%vertices(2, 4))
        allocate(mesh%triangles(3, 2))
        
        ! Four vertices
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 1.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 1.0_dp
        mesh%vertices(1, 4) = 1.0_dp
        mesh%vertices(2, 4) = 1.0_dp
        
        ! Two triangles sharing edge 2-3
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
        mesh%triangles(1, 2) = 2
        mesh%triangles(2, 2) = 4
        mesh%triangles(3, 2) = 3
    end subroutine

end program test_global_curl_curl_assembly
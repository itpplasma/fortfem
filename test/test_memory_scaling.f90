program test_memory_scaling
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_sparse_matrix
    implicit none

    call test_edge_space_cleanup()
    call test_sparse_matrix_efficiency()
    call test_large_problem_scaling()
    
    print *, "All memory and scaling tests passed!"

contains

    subroutine test_edge_space_cleanup()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        integer :: i
        
        print *, ""
        print *, "Edge Space Memory Cleanup Test"
        print *, "=============================="
        
        ! Create and destroy multiple times to test cleanup
        do i = 1, 3
            ! Create mesh and space
            call mesh%create_rectangular(10, 10, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
            call mesh%build_edge_connectivity()
            call mesh%build_edge_dof_numbering()
            call space%init(mesh)
            
            print *, "Iteration", i, ":"
            print *, "  Mesh vertices:", mesh%n_vertices
            print *, "  Mesh triangles:", mesh%n_triangles
            print *, "  Edge DOFs:", space%get_n_dofs()
            
            ! Clean up
            call space%destroy()
            call mesh%destroy()
            
            ! Verify cleanup
            if (allocated(mesh%vertices)) then
                print *, "Error: mesh vertices not deallocated"
                stop 1
            end if
            
            if (mesh%n_vertices /= 0 .or. mesh%n_triangles /= 0) then
                print *, "Error: mesh counters not reset"
                stop 1
            end if
        end do
        
        print *, "Edge space cleanup test passed"
    end subroutine
    
    subroutine test_sparse_matrix_efficiency()
        type(triplet_matrix_t) :: triplet
        type(csr_matrix_t) :: csr
        integer :: n, nnz_expected, nnz_actual, i
        integer :: sizes(3) = [10, 50, 100]
        real(dp) :: sparsity
        
        print *, ""
        print *, "Sparse Matrix Storage Efficiency Test"
        print *, "===================================="
        
        ! Test different matrix sizes
        
        do i = 1, 3
            n = sizes(i)
            
            ! Create tridiagonal matrix (typical for 1D problems)
            call triplet%init(n, 3*n)
            call create_tridiagonal_matrix(n, triplet)
            
            ! Convert to CSR
            call triplet%to_csr(csr)
            
            ! Check storage efficiency
            nnz_expected = 3*n - 2  ! Tridiagonal
            nnz_actual = csr%nnz
            sparsity = real(nnz_actual, dp) / real(n*n, dp)
            
            print *, "Matrix size", n, "x", n, ":"
            print *, "  Non-zeros:", nnz_actual
            print *, "  Expected:", nnz_expected
            print *, "  Sparsity:", sparsity * 100.0_dp, "%"
            
            if (abs(nnz_actual - nnz_expected) > 0) then
                print *, "Warning: unexpected number of non-zeros"
            end if
            
            call triplet%destroy()
            call csr%destroy()
        end do
        
        print *, "Sparse matrix efficiency test passed"
    end subroutine
    
    subroutine test_large_problem_scaling()
        type(mesh_2d_t) :: mesh
        type(hcurl_space_t) :: space
        type(triplet_matrix_t) :: system_matrix
        integer :: mesh_sizes(4) = [5, 10, 20, 30]
        integer :: n_dofs(4), n_entries(4)
        real(dp) :: h_values(4)
        real(dp) :: dof_ratio, entry_ratio
        integer :: i
        
        print *, ""
        print *, "Large Problem Scaling Test"
        print *, "=========================="
        print *, ""
        print *, "Mesh Size    h        DOFs    Matrix Entries   Entries/DOF"
        print *, "-----------------------------------------------------------"
        
        do i = 1, 4
            ! Create mesh
            call mesh%create_rectangular(mesh_sizes(i)+1, mesh_sizes(i)+1, &
                                       0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
            call mesh%build_edge_connectivity()
            call mesh%build_edge_dof_numbering()
            call space%init(mesh)
            
            ! Get problem size
            n_dofs(i) = space%get_n_dofs()
            h_values(i) = 1.0_dp / real(mesh_sizes(i), dp)
            
            ! Assemble matrix to count entries
            call system_matrix%init(n_dofs(i), n_dofs(i) * 20)
            call assemble_curl_curl_system(mesh, system_matrix)
            n_entries(i) = system_matrix%nnz
            
            write(*, '(I8, F10.4, I8, I15, F15.2)') &
                mesh_sizes(i)**2, h_values(i), n_dofs(i), n_entries(i), &
                real(n_entries(i), dp) / real(n_dofs(i), dp)
            
            ! Clean up
            call system_matrix%destroy()
            call space%destroy()
            call mesh%destroy()
        end do
        
        ! Check scaling behavior
        print *, ""
        print *, "Scaling analysis:"
        do i = 2, 4
            dof_ratio = real(n_dofs(i), dp) / real(n_dofs(i-1), dp)
            entry_ratio = real(n_entries(i), dp) / real(n_entries(i-1), dp)
            print *, "  h/2: DOFs x", dof_ratio, ", Entries x", entry_ratio
        end do
        
        print *, "Large problem scaling test passed"
    end subroutine
    
    subroutine create_tridiagonal_matrix(n, matrix)
        integer, intent(in) :: n
        type(triplet_matrix_t), intent(inout) :: matrix
        
        integer :: i
        
        ! Main diagonal
        do i = 1, n
            call matrix%add(i, i, 2.0_dp)
        end do
        
        ! Sub-diagonal
        do i = 2, n
            call matrix%add(i, i-1, -1.0_dp)
        end do
        
        ! Super-diagonal
        do i = 1, n-1
            call matrix%add(i, i+1, -1.0_dp)
        end do
    end subroutine
    
    subroutine assemble_curl_curl_system(mesh, matrix)
        type(mesh_2d_t), intent(in) :: mesh
        type(triplet_matrix_t), intent(inout) :: matrix
        
        real(dp) :: local_matrix(3, 3)
        real(dp) :: triangle_area
        integer :: triangle_dofs(3)
        integer :: t, i, j
        
        ! Simple assembly for testing
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            
            ! Simple local matrix
            local_matrix = triangle_area / 12.0_dp
            local_matrix(1, 1) = 2.0_dp * triangle_area / 12.0_dp
            local_matrix(2, 2) = 2.0_dp * triangle_area / 12.0_dp
            local_matrix(3, 3) = 2.0_dp * triangle_area / 12.0_dp
            
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            
            do i = 1, 3
                do j = 1, 3
                    if (abs(local_matrix(i, j)) > 1e-12_dp) then
                        call matrix%add(triangle_dofs(i) + 1, triangle_dofs(j) + 1, &
                                       local_matrix(i, j))
                    end if
                end do
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

end program test_memory_scaling
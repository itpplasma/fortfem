program verify_assembly
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(assembly_2d_t) :: assembly
    type(triplet_matrix_t) :: triplet
    type(csr_matrix_t) :: csr
    real(dp), allocatable :: rhs(:), solution(:)
    integer :: i, j
    
    print *, "Global Assembly Verification"
    print *, "============================"
    print *, ""
    
    ! Create minimal 2x2 mesh (2 triangles)
    call mesh%create_rectangular(nx=2, ny=2, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Print mesh details
    print *, ""
    print *, "Mesh vertices:"
    do i = 1, mesh%n_vertices
        print '(a,i0,a,2f8.3)', "  V", i, ": ", mesh%vertices(:,i)
    end do
    
    print *, ""
    print *, "Mesh triangles:"
    do i = 1, mesh%n_triangles
        print '(a,i0,a,3i3)', "  T", i, ": ", mesh%triangles(:,i)
    end do
    
    ! Initialize assembly
    call assembly%init(mesh%n_vertices, 9 * mesh%n_vertices)
    
    ! Assemble global system
    allocate(rhs(mesh%n_vertices))
    call assembly%assemble_global(mesh, unit_source, triplet, rhs)
    
    print *, ""
    print '(a,i0)', "Matrix entries: ", triplet%nnz
    print '(a,f10.6)', "RHS sum: ", sum(rhs)
    print '(a,f10.6)', "Expected RHS (area): ", 1.0_dp
    
    ! Show matrix entries
    print *, ""
    print *, "Matrix entries (row, col, value):"
    do i = 1, min(20, triplet%nnz)
        print '(2i5,f12.6)', triplet%rows(i), triplet%cols(i), triplet%values(i)
    end do
    
    ! Convert to CSR and verify
    call triplet%to_csr(csr)
    print *, ""
    print '(a,i0)', "CSR matrix size: ", csr%n
    print '(a,i0)', "CSR nnz: ", csr%nnz
    
    ! Test matrix-vector multiplication
    allocate(solution(mesh%n_vertices))
    solution = 1.0_dp
    call csr%matvec(solution, rhs)
    
    print *, ""
    print *, "Matrix-vector test (x=1):"
    print '(a,f10.6)', "  Result sum: ", sum(rhs)
    print '(a)', "  (Should be sum of matrix rows)"
    
    ! Test simple solve problem
    print *, ""
    print *, "Testing simple problem with known solution..."
    
    ! For this simple mesh, manually set up a problem
    ! Let's solve with RHS = [1,1,1,1] and see if we get reasonable results
    rhs = 1.0_dp
    
    ! Apply boundary conditions (set boundary nodes to 0)
    ! For 2x2 mesh, boundary nodes are typically edges
    do i = 1, mesh%n_vertices
        if (abs(mesh%vertices(1,i)) < 1e-10 .or. abs(mesh%vertices(1,i) - 1.0_dp) < 1e-10 .or. &
            abs(mesh%vertices(2,i)) < 1e-10 .or. abs(mesh%vertices(2,i) - 1.0_dp) < 1e-10) then
            
            ! Set diagonal to 1, off-diagonal to 0, RHS to 0
            do j = 1, csr%nnz
                if (csr%col_idx(j) == i) then
                    ! This is column i
                    if (csr%row_ptr(i) <= j .and. j < csr%row_ptr(i+1)) then
                        ! This is also row i (diagonal)
                        csr%values(j) = 1.0_dp
                    else
                        ! Off-diagonal
                        csr%values(j) = 0.0_dp
                    end if
                end if
            end do
            
            ! Zero row i
            do j = csr%row_ptr(i), csr%row_ptr(i+1)-1
                if (csr%col_idx(j) == i) then
                    csr%values(j) = 1.0_dp
                else
                    csr%values(j) = 0.0_dp
                end if
            end do
            
            rhs(i) = 0.0_dp
        end if
    end do
    
    print *, "Matrix after boundary conditions:"
    do i = 1, csr%n
        print '(a,i0,a)', "  Row ", i, ":"
        do j = csr%row_ptr(i), csr%row_ptr(i+1)-1
            print '(a,i0,a,f10.6)', "    (", csr%col_idx(j), ") = ", csr%values(j)
        end do
    end do
    
    print *, ""
    print *, "RHS after boundary conditions:"
    do i = 1, csr%n
        print '(a,i0,a,f10.6)', "  b[", i, "] = ", rhs(i)
    end do
    
    ! Clean up
    call assembly%destroy()
    call triplet%destroy()
    call csr%destroy()
    call mesh%destroy()
    deallocate(rhs, solution)
    
    print *, ""
    print *, "Assembly verification complete!"
    
contains

    pure function unit_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 1.0_dp
    end function unit_source

end program verify_assembly
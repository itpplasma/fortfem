program test_assembly_2d_example
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(assembly_2d_t) :: assembly
    type(triplet_matrix_t) :: triplet
    real(dp), allocatable :: rhs(:)
    real(dp) :: vertices(2,3), mass(3,3), stiff(3,3), load(3)
    integer :: i, j
    
    print *, "2D Assembly Test"
    print *, "================"
    
    ! Test 1: Element matrices for unit triangle
    print *, ""
    print *, "Test 1: Element matrices for unit right triangle"
    vertices(:,1) = [0.0_dp, 0.0_dp]
    vertices(:,2) = [1.0_dp, 0.0_dp]
    vertices(:,3) = [0.0_dp, 1.0_dp]
    
    ! Mass matrix
    call assembly%element_mass_matrix(vertices, mass)
    print *, "Mass matrix:"
    do i = 1, 3
        print '(3f10.6)', mass(i,:)
    end do
    print '(a,f10.6)', "Total mass (sum of entries): ", sum(mass)
    
    ! Stiffness matrix
    call assembly%element_stiffness_matrix(vertices, stiff)
    print *, ""
    print *, "Stiffness matrix:"
    do i = 1, 3
        print '(3f10.6)', stiff(i,:)
    end do
    print '(a,f10.6)', "Row sums (should be ~0): ", &
        maxval(abs([sum(stiff(1,:)), sum(stiff(2,:)), sum(stiff(3,:))]))
    
    ! Load vector
    call assembly%element_load_vector(vertices, const_one, load)
    print *, ""
    print *, "Load vector (f=1):"
    print '(3f10.6)', load
    print '(a,f10.6)', "Total load: ", sum(load)
    
    ! Test 2: Small mesh assembly
    print *, ""
    print *, "Test 2: Assembly on 2x2 mesh"
    
    call mesh%create_rectangular(nx=2, ny=2, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Number of vertices: ", mesh%n_vertices
    print '(a,i0)', "Number of triangles: ", mesh%n_triangles
    
    ! Initialize assembly
    call assembly%init(mesh%n_vertices, 9 * mesh%n_vertices)
    
    ! Assemble global system
    allocate(rhs(mesh%n_vertices))
    call assembly%assemble_global(mesh, const_one, triplet, rhs)
    
    print '(a,i0)', "Number of matrix entries: ", triplet%nnz
    print '(a,f10.6)', "Sum of RHS: ", sum(rhs)
    print '(a,f10.6)', "Max RHS value: ", maxval(rhs)
    print '(a,f10.6)', "Min RHS value: ", minval(rhs)
    
    ! Show first few matrix entries
    print *, ""
    print *, "First 10 matrix entries (row, col, value):"
    do i = 1, min(10, triplet%nnz)
        print '(2i5,f12.6)', triplet%rows(i), triplet%cols(i), triplet%values(i)
    end do
    
    ! Clean up
    call assembly%destroy()
    call triplet%destroy()
    deallocate(rhs)
    call mesh%destroy()
    
    print *, ""
    print *, "All tests completed!"
    
contains

    pure function const_one(x, y) result(val)
        real(dp), intent(in) :: x, y
        real(dp) :: val
        val = 1.0_dp
    end function const_one

end program test_assembly_2d_example
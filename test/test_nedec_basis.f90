program test_nedec_basis
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    use fortfem_mesh_2d
    implicit none

    call test_nedec_rt0_basis_properties()
    call test_nedec_rt0_basis_values()
    
    print *, "All Nédélec RT0 basis tests passed!"

contains

    subroutine test_nedec_rt0_basis_properties()
        type(edge_basis_2d_t) :: basis
        type(mesh_2d_t) :: mesh
        integer :: n_dofs
        
        ! Create reference triangle mesh
        call create_reference_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Initialize edge basis
        call basis%init(mesh)
        
        ! Test number of DOFs (should be 3 for RT0)
        n_dofs = basis%n_dofs()
        if (n_dofs /= 3) then
            print *, "Error: expected 3 DOFs, got", n_dofs
            stop 1
        end if
        
        call basis%destroy()
        call mesh%destroy()
        
        print *, "Nédélec RT0 basis properties test passed"
    end subroutine
    
    subroutine test_nedec_rt0_basis_values()
        type(edge_basis_2d_t) :: basis
        type(mesh_2d_t) :: mesh
        real(dp) :: xi, eta, triangle_area
        real(dp) :: values(2, 3)  ! 2 components, 3 basis functions
        real(dp) :: curls(3)      ! Curl of each basis function
        real(dp) :: divs(3)       ! Divergence of each basis function
        integer :: i
        
        ! Create reference triangle mesh
        call create_reference_triangle_mesh(mesh)
        call mesh%build_edge_connectivity()
        
        ! Initialize edge basis
        call basis%init(mesh)
        
        ! Test at triangle center (barycentric coordinates (1/3, 1/3))
        xi = 1.0_dp / 3.0_dp
        eta = 1.0_dp / 3.0_dp
        triangle_area = 0.5_dp  ! Area of reference triangle
        
        ! Evaluate basis functions
        call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
        call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
        call evaluate_edge_basis_div_2d(xi, eta, triangle_area, divs)
        
        ! Check that values are reasonable (non-zero for RT0)
        do i = 1, 3
            if (abs(values(1, i)) < 1e-16_dp .and. abs(values(2, i)) < 1e-16_dp) then
                print *, "Warning: basis function", i, "has zero value at center"
            end if
        end do
        
        ! Check that curls are constant (RT0 property)
        do i = 1, 3
            if (abs(curls(i)) < 1e-16_dp) then
                print *, "Warning: basis function", i, "has zero curl"
            end if
        end do
        
        ! For RT0, div should be zero for all basis functions
        do i = 1, 3
            if (abs(divs(i)) > 1e-12_dp) then
                print *, "Error: RT0 basis function", i, "has non-zero divergence:", divs(i)
                stop 1
            end if
        end do
        
        call basis%destroy()
        call mesh%destroy()
        
        print *, "Nédélec RT0 basis values test passed"
    end subroutine
    
    subroutine create_reference_triangle_mesh(mesh)
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

end program test_nedec_basis
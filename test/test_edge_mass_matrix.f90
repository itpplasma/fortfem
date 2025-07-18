program test_edge_mass_matrix
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none

    call test_local_mass_matrix()
    call test_mass_matrix_properties()
    call test_vector_valued_integration()
    
    print *, "All edge element mass matrix tests passed!"

contains

    subroutine test_local_mass_matrix()
        type(mesh_2d_t) :: mesh
        real(dp) :: local_mass_matrix(3, 3)
        real(dp) :: triangle_area
        integer :: i, j
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        triangle_area = 0.5_dp  ! Reference triangle area
        
        ! Compute local mass matrix: M_ij = ∫_T φ_i · φ_j dx
        ! For RT0, this requires integration over the triangle
        call compute_local_mass_matrix(triangle_area, local_mass_matrix)
        
        ! Verify symmetry
        do i = 1, 3
            do j = 1, 3
                if (abs(local_mass_matrix(i, j) - local_mass_matrix(j, i)) > 1e-12_dp) then
                    print *, "Error: mass matrix not symmetric at (", i, ",", j, ")"
                    stop 1
                end if
            end do
        end do
        
        ! Verify positive definiteness (diagonal entries should be positive)
        do i = 1, 3
            if (local_mass_matrix(i, i) <= 0.0_dp) then
                print *, "Error: mass matrix diagonal entry", i, "not positive:", local_mass_matrix(i, i)
                stop 1
            end if
        end do
        
        call mesh%destroy()
        print *, "Local mass matrix test passed"
    end subroutine
    
    subroutine test_mass_matrix_properties()
        type(mesh_2d_t) :: mesh
        real(dp) :: local_mass_matrix(3, 3)
        real(dp) :: triangle_area
        real(dp) :: test_vector(3), quadratic_form
        integer :: i, j
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        triangle_area = 0.5_dp
        call compute_local_mass_matrix(triangle_area, local_mass_matrix)
        
        ! Test positive definiteness with several test vectors
        
        ! Test vector 1: [1, 0, 0]
        test_vector = [1.0_dp, 0.0_dp, 0.0_dp]
        quadratic_form = 0.0_dp
        do i = 1, 3
            do j = 1, 3
                quadratic_form = quadratic_form + test_vector(i) * local_mass_matrix(i, j) * test_vector(j)
            end do
        end do
        
        if (quadratic_form <= 0.0_dp) then
            print *, "Error: mass matrix not positive definite for test vector [1,0,0]"
            stop 1
        end if
        
        ! Test vector 2: [1, 1, 1]
        test_vector = [1.0_dp, 1.0_dp, 1.0_dp]
        quadratic_form = 0.0_dp
        do i = 1, 3
            do j = 1, 3
                quadratic_form = quadratic_form + test_vector(i) * local_mass_matrix(i, j) * test_vector(j)
            end do
        end do
        
        if (quadratic_form <= 0.0_dp) then
            print *, "Error: mass matrix not positive definite for test vector [1,1,1]"
            stop 1
        end if
        
        call mesh%destroy()
        print *, "Mass matrix properties test passed"
    end subroutine
    
    subroutine test_vector_valued_integration()
        type(mesh_2d_t) :: mesh
        real(dp) :: triangle_area
        real(dp) :: xi, eta
        real(dp) :: basis_values_1(2, 3), basis_values_2(2, 3)
        real(dp) :: dot_product
        integer :: i, j
        
        ! Create reference triangle
        call create_reference_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        triangle_area = 0.5_dp
        
        ! Test vector-valued integration at different points
        xi = 0.25_dp
        eta = 0.25_dp
        
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values_1)
        
        xi = 0.5_dp
        eta = 0.25_dp
        
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values_2)
        
        ! Compute dot products between basis functions
        do i = 1, 3
            do j = 1, 3
                dot_product = basis_values_1(1, i) * basis_values_2(1, j) + &
                             basis_values_1(2, i) * basis_values_2(2, j)
                
                ! Verify that dot products are well-defined
                if (abs(dot_product) > 1e10_dp) then
                    print *, "Error: dot product too large - numerical issues"
                    stop 1
                end if
            end do
        end do
        
        call mesh%destroy()
        print *, "Vector-valued integration test passed"
    end subroutine
    
    subroutine compute_local_mass_matrix(triangle_area, mass_matrix)
        real(dp), intent(in) :: triangle_area
        real(dp), intent(out) :: mass_matrix(3, 3)
        
        ! Quadrature points and weights for triangle integration
        real(dp), parameter :: xi_quad(3) = [0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: eta_quad(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
        real(dp) :: basis_values(2, 3)
        real(dp) :: dot_product
        integer :: i, j, q
        
        mass_matrix = 0.0_dp
        
        ! Integrate using quadrature
        do q = 1, 3
            call evaluate_edge_basis_2d(xi_quad(q), eta_quad(q), triangle_area, basis_values)
            
            do i = 1, 3
                do j = 1, 3
                    ! Compute φ_i · φ_j
                    dot_product = basis_values(1, i) * basis_values(1, j) + &
                                 basis_values(2, i) * basis_values(2, j)
                    
                    ! Add weighted contribution
                    mass_matrix(i, j) = mass_matrix(i, j) + weights(q) * dot_product * triangle_area
                end do
            end do
        end do
    end subroutine compute_local_mass_matrix
    
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

end program test_edge_mass_matrix
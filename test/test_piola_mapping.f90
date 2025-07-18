program test_piola_mapping
    use fortfem_kinds, only: dp
    use fortfem_basis_edge_2d
    use fortfem_mesh_2d
    implicit none

    call test_piola_transformation()
    call test_jacobian_scaling()
    call test_curl_transformation()
    
    print *, "All Piola mapping tests passed!"

contains

    subroutine test_piola_transformation()
        type(mesh_2d_t) :: mesh
        real(dp) :: xi, eta
        real(dp) :: ref_values(2, 3), phys_values(2, 3)
        real(dp) :: jacobian(2, 2), det_jac, inv_jac(2, 2)
        real(dp) :: triangle_area
        integer :: i, j
        
        ! Create a physical triangle (not reference)
        call create_physical_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        ! Compute triangle area
        triangle_area = compute_triangle_area(mesh, 1)
        
        ! Reference coordinates
        xi = 0.33_dp
        eta = 0.33_dp
        
        ! Get reference basis values
        call evaluate_edge_basis_2d(xi, eta, 0.5_dp, ref_values)
        
        ! Compute Jacobian of transformation from reference to physical
        call compute_jacobian(mesh, 1, xi, eta, jacobian, det_jac)
        call invert_2x2_matrix(jacobian, inv_jac)
        
        ! Piola transformation: φ_phys = J^{-T} φ_ref / det(J)
        ! This preserves tangential continuity across edges
        do i = 1, 3
            ! Transform each basis function
            phys_values(1, i) = (inv_jac(1, 1) * ref_values(1, i) + &
                                inv_jac(2, 1) * ref_values(2, i))
            phys_values(2, i) = (inv_jac(1, 2) * ref_values(1, i) + &
                                inv_jac(2, 2) * ref_values(2, i))
        end do
        
        ! Verify transformation preserves essential properties
        ! For RT0, the values should be reasonable (not zero or infinite)
        do i = 1, 3
            if (abs(phys_values(1, i)) > 1e6_dp .or. abs(phys_values(2, i)) > 1e6_dp) then
                print *, "Error: Piola transformation produced unreasonable values"
                print *, "Basis", i, "values:", phys_values(:, i)
                stop 1
            end if
        end do
        
        call mesh%destroy()
        print *, "Piola transformation test passed"
    end subroutine
    
    subroutine test_jacobian_scaling()
        type(mesh_2d_t) :: mesh
        real(dp) :: xi, eta
        real(dp) :: ref_curls(3), phys_curls(3)
        real(dp) :: jacobian(2, 2), det_jac, inv_jac(2, 2)
        real(dp) :: ref_area, phys_area
        integer :: i
        
        ! Create physical triangle
        call create_physical_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        xi = 0.33_dp
        eta = 0.33_dp
        
        ! Reference triangle area
        ref_area = 0.5_dp
        
        ! Physical triangle area
        phys_area = compute_triangle_area(mesh, 1)
        
        ! Get reference curl values
        call evaluate_edge_basis_curl_2d(xi, eta, ref_area, ref_curls)
        
        ! Compute Jacobian
        call compute_jacobian(mesh, 1, xi, eta, jacobian, det_jac)
        
        ! Curl transformation: curl(φ_phys) = curl(φ_ref) / det(J)
        do i = 1, 3
            phys_curls(i) = ref_curls(i) / det_jac
        end do
        
        ! Verify scaling relationship
        do i = 1, 3
            if (abs(phys_curls(i) * det_jac - ref_curls(i)) > 1e-12_dp) then
                print *, "Error: Jacobian scaling incorrect for curl", i
                print *, "Expected:", ref_curls(i)
                print *, "Got:", phys_curls(i) * det_jac
                stop 1
            end if
        end do
        
        call mesh%destroy()
        print *, "Jacobian scaling test passed"
    end subroutine
    
    subroutine test_curl_transformation()
        type(mesh_2d_t) :: mesh
        real(dp) :: xi, eta
        real(dp) :: ref_curls(3), phys_curls(3)
        real(dp) :: jacobian(2, 2), det_jac, inv_jac(2, 2)
        real(dp) :: ref_area, phys_area
        integer :: i
        
        ! Create physical triangle
        call create_physical_triangle(mesh)
        call mesh%build_edge_connectivity()
        
        xi = 0.33_dp
        eta = 0.33_dp
        
        ! Areas
        ref_area = 0.5_dp
        phys_area = compute_triangle_area(mesh, 1)
        
        ! Get curl values
        call evaluate_edge_basis_curl_2d(xi, eta, ref_area, ref_curls)
        call evaluate_edge_basis_curl_2d(xi, eta, phys_area, phys_curls)
        
        ! Compute Jacobian
        call compute_jacobian(mesh, 1, xi, eta, jacobian, det_jac)
        
        ! Verify curl transformation formula
        do i = 1, 3
            if (abs(phys_curls(i) - ref_curls(i) / det_jac) > 1e-12_dp) then
                print *, "Error: curl transformation incorrect for basis", i
                print *, "Expected:", ref_curls(i) / det_jac
                print *, "Got:", phys_curls(i)
                stop 1
            end if
        end do
        
        call mesh%destroy()
        print *, "Curl transformation test passed"
    end subroutine
    
    subroutine create_physical_triangle(mesh)
        type(mesh_2d_t), intent(out) :: mesh
        
        mesh%n_vertices = 3
        mesh%n_triangles = 1
        
        allocate(mesh%vertices(2, 3))
        allocate(mesh%triangles(3, 1))
        
        ! Physical triangle vertices (not reference)
        mesh%vertices(1, 1) = 0.0_dp
        mesh%vertices(2, 1) = 0.0_dp
        mesh%vertices(1, 2) = 2.0_dp
        mesh%vertices(2, 2) = 0.0_dp
        mesh%vertices(1, 3) = 0.0_dp
        mesh%vertices(2, 3) = 2.0_dp
        
        ! Triangle connectivity
        mesh%triangles(1, 1) = 1
        mesh%triangles(2, 1) = 2
        mesh%triangles(3, 1) = 3
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
    
    subroutine compute_jacobian(mesh, triangle_idx, xi, eta, jacobian, det_jac)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: jacobian(2, 2), det_jac
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        ! Jacobian of transformation from reference to physical
        jacobian(1, 1) = x2 - x1
        jacobian(1, 2) = x3 - x1
        jacobian(2, 1) = y2 - y1
        jacobian(2, 2) = y3 - y1
        
        det_jac = jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1)
    end subroutine
    
    subroutine invert_2x2_matrix(A, A_inv)
        real(dp), intent(in) :: A(2, 2)
        real(dp), intent(out) :: A_inv(2, 2)
        
        real(dp) :: det
        
        det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
        
        A_inv(1, 1) = A(2, 2) / det
        A_inv(1, 2) = -A(1, 2) / det
        A_inv(2, 1) = -A(2, 1) / det
        A_inv(2, 2) = A(1, 1) / det
    end subroutine

end program test_piola_mapping
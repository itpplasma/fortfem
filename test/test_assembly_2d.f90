module test_assembly_2d
    use fortfem
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private
    
    public :: collect_assembly_2d
    
contains

    subroutine collect_assembly_2d(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        
        testsuite = [ &
            new_unittest("element_mass_matrix", test_element_mass_matrix), &
            new_unittest("element_stiffness_matrix", test_element_stiffness_matrix), &
            new_unittest("element_load_vector", test_element_load_vector), &
            new_unittest("mass_matrix_properties", test_mass_matrix_properties), &
            new_unittest("stiffness_matrix_properties", test_stiffness_matrix_properties), &
            new_unittest("assembly_small_mesh", test_assembly_small_mesh) &
        ]
        
    end subroutine collect_assembly_2d
    
    subroutine test_element_mass_matrix(error)
        type(error_type), allocatable, intent(out) :: error
        type(assembly_2d_t) :: assembly
        real(dp) :: vertices(2,3), mass(3,3)
        real(dp) :: area, expected_diag, expected_off
        
        ! Unit right triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        ! Compute element mass matrix
        call assembly%element_mass_matrix(vertices, mass)
        
        ! For P1 elements: M_ij = area/12 for i != j, area/6 for i = j
        area = 0.5_dp
        expected_diag = area / 6.0_dp
        expected_off = area / 12.0_dp
        
        ! Check diagonal elements
        call check(error, abs(mass(1,1) - expected_diag) < 1e-14, &
                   "Mass matrix diagonal incorrect")
        call check(error, abs(mass(2,2) - expected_diag) < 1e-14, &
                   "Mass matrix diagonal incorrect")
        call check(error, abs(mass(3,3) - expected_diag) < 1e-14, &
                   "Mass matrix diagonal incorrect")
        
        ! Check off-diagonal elements
        call check(error, abs(mass(1,2) - expected_off) < 1e-14, &
                   "Mass matrix off-diagonal incorrect")
        call check(error, abs(mass(2,3) - expected_off) < 1e-14, &
                   "Mass matrix off-diagonal incorrect")
        
        ! Check symmetry
        call check(error, abs(mass(1,2) - mass(2,1)) < 1e-14, &
                   "Mass matrix not symmetric")
        
    end subroutine test_element_mass_matrix
    
    subroutine test_element_stiffness_matrix(error)
        type(error_type), allocatable, intent(out) :: error
        type(assembly_2d_t) :: assembly
        real(dp) :: vertices(2,3), stiff(3,3)
        real(dp) :: sum_row
        integer :: i
        
        ! Unit right triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        ! Compute element stiffness matrix
        call assembly%element_stiffness_matrix(vertices, stiff)
        
        ! Check symmetry
        do i = 1, 3
            call check(error, abs(stiff(i,1) - stiff(1,i)) < 1e-14, &
                       "Stiffness matrix not symmetric")
            call check(error, abs(stiff(i,2) - stiff(2,i)) < 1e-14, &
                       "Stiffness matrix not symmetric")
            call check(error, abs(stiff(i,3) - stiff(3,i)) < 1e-14, &
                       "Stiffness matrix not symmetric")
            if (allocated(error)) return
        end do
        
        ! Check that rows sum to zero (constant function has zero energy)
        do i = 1, 3
            sum_row = stiff(i,1) + stiff(i,2) + stiff(i,3)
            call check(error, abs(sum_row) < 1e-14, &
                       "Stiffness matrix rows should sum to zero")
            if (allocated(error)) return
        end do
        
    end subroutine test_element_stiffness_matrix
    
    subroutine test_element_load_vector(error)
        type(error_type), allocatable, intent(out) :: error
        type(assembly_2d_t) :: assembly
        real(dp) :: vertices(2,3), load(3)
        real(dp) :: area, expected_val
        
        ! Unit right triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [1.0_dp, 0.0_dp]
        vertices(:,3) = [0.0_dp, 1.0_dp]
        
        ! Compute load vector for f = 1
        call assembly%element_load_vector(vertices, const_one, load)
        
        ! For constant source, each entry should be area/3
        area = 0.5_dp
        expected_val = area / 3.0_dp
        
        call check(error, abs(load(1) - expected_val) < 1e-14, &
                   "Load vector incorrect for constant source")
        call check(error, abs(load(2) - expected_val) < 1e-14, &
                   "Load vector incorrect for constant source")
        call check(error, abs(load(3) - expected_val) < 1e-14, &
                   "Load vector incorrect for constant source")
        
    end subroutine test_element_load_vector
    
    subroutine test_mass_matrix_properties(error)
        type(error_type), allocatable, intent(out) :: error
        type(assembly_2d_t) :: assembly
        real(dp) :: vertices(2,3), mass(3,3)
        real(dp) :: eigenvalues(3), work(9)
        integer :: info
        
        ! General triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [1.0_dp, 1.5_dp]
        
        call assembly%element_mass_matrix(vertices, mass)
        
        ! Check positive definiteness by computing eigenvalues
        ! For symmetric matrix, we can use DSYEV
        call dsyev('N', 'U', 3, mass, 3, eigenvalues, work, 9, info)
        
        call check(error, info == 0, "Eigenvalue computation failed")
        if (allocated(error)) return
        
        ! All eigenvalues should be positive
        call check(error, all(eigenvalues > 0.0_dp), &
                   "Mass matrix not positive definite")
        
    end subroutine test_mass_matrix_properties
    
    subroutine test_stiffness_matrix_properties(error)
        type(error_type), allocatable, intent(out) :: error
        type(assembly_2d_t) :: assembly
        real(dp) :: vertices(2,3), stiff(3,3)
        real(dp) :: eigenvalues(3), work(9)
        integer :: info, i
        
        ! General triangle
        vertices(:,1) = [0.0_dp, 0.0_dp]
        vertices(:,2) = [2.0_dp, 0.0_dp]
        vertices(:,3) = [1.0_dp, 1.5_dp]
        
        call assembly%element_stiffness_matrix(vertices, stiff)
        
        ! Check eigenvalues - should have one zero eigenvalue (constants)
        call dsyev('N', 'U', 3, stiff, 3, eigenvalues, work, 9, info)
        
        call check(error, info == 0, "Eigenvalue computation failed")
        if (allocated(error)) return
        
        ! One eigenvalue should be ~0, others positive
        call check(error, abs(eigenvalues(1)) < 1e-12, &
                   "Stiffness matrix should have one zero eigenvalue")
        call check(error, eigenvalues(2) > 1e-12, &
                   "Stiffness matrix should have positive eigenvalues")
        call check(error, eigenvalues(3) > 1e-12, &
                   "Stiffness matrix should have positive eigenvalues")
        
    end subroutine test_stiffness_matrix_properties
    
    subroutine test_assembly_small_mesh(error)
        type(error_type), allocatable, intent(out) :: error
        type(mesh_2d_t) :: mesh
        type(assembly_2d_t) :: assembly
        type(triplet_matrix_t) :: triplet
        real(dp), allocatable :: rhs(:)
        integer :: nnz_estimate
        
        ! Create small 2x2 mesh
        call mesh%create_rectangular(nx=2, ny=2, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        
        ! Estimate non-zeros (upper bound)
        nnz_estimate = 9 * mesh%n_vertices
        
        ! Initialize assembly
        call assembly%init(mesh%n_vertices, nnz_estimate)
        
        ! Assemble global system
        allocate(rhs(mesh%n_vertices))
        call assembly%assemble_global(mesh, const_one, triplet, rhs)
        
        ! Check that we got some entries
        call check(error, triplet%nnz > 0, "No entries assembled")
        call check(error, triplet%nnz <= nnz_estimate, "Too many entries")
        
        ! Check RHS is positive (for positive source)
        call check(error, all(rhs > 0.0_dp), "RHS should be positive")
        
        ! Clean up
        call assembly%destroy()
        call triplet%destroy()
        deallocate(rhs)
        call mesh%destroy()
        
    end subroutine test_assembly_small_mesh
    
    ! Helper function for constant source
    pure function const_one(x, y) result(val)
        real(dp), intent(in) :: x, y
        real(dp) :: val
        val = 1.0_dp
    end function const_one

end module test_assembly_2d
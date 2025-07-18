module fortfem_assembly_1d
    use fortfem_kinds
    use fortfem_basis_1d
    implicit none
    private
    
    public :: element_mass_matrix
    public :: element_stiffness_matrix
    public :: element_load_vector
    public :: source_function
    
    ! Interface for source term functions
    abstract interface
        function source_function(x) result(f)
            import :: dp
            real(dp), intent(in) :: x
            real(dp) :: f
        end function source_function
    end interface
    
contains

    subroutine element_mass_matrix(h, M_elem)
        real(dp), intent(in) :: h          ! Element length
        real(dp), intent(out) :: M_elem(2,2)
        
        ! For P1 elements, exact integration gives:
        ! M = h/6 * [2 1; 1 2]
        M_elem(1,1) = h / 3.0_dp
        M_elem(1,2) = h / 6.0_dp
        M_elem(2,1) = h / 6.0_dp
        M_elem(2,2) = h / 3.0_dp
        
    end subroutine element_mass_matrix
    
    subroutine element_stiffness_matrix(h, K_elem)
        real(dp), intent(in) :: h          ! Element length
        real(dp), intent(out) :: K_elem(2,2)
        
        ! For P1 elements with constant Jacobian:
        ! K = 1/h * [1 -1; -1 1]
        K_elem(1,1) =  1.0_dp / h
        K_elem(1,2) = -1.0_dp / h
        K_elem(2,1) = -1.0_dp / h
        K_elem(2,2) =  1.0_dp / h
        
    end subroutine element_stiffness_matrix
    
    subroutine element_load_vector(h, source, f_elem)
        real(dp), intent(in) :: h          ! Element length
        procedure(source_function) :: source
        real(dp), intent(out) :: f_elem(2)
        
        ! 2-point Gauss quadrature for P1 elements
        real(dp), parameter :: xi1 = 0.5_dp - 1.0_dp/(2.0_dp*sqrt(3.0_dp))
        real(dp), parameter :: xi2 = 0.5_dp + 1.0_dp/(2.0_dp*sqrt(3.0_dp))
        real(dp), parameter :: w = 0.5_dp
        
        real(dp) :: x1, x2
        
        ! Map quadrature points to physical element
        ! Assuming element from x to x+h
        x1 = xi1 * h
        x2 = xi2 * h
        
        ! Integrate source * basis functions
        f_elem(1) = h * w * (source(x1) * p1_basis(1, xi1) + &
                            source(x2) * p1_basis(1, xi2))
        f_elem(2) = h * w * (source(x1) * p1_basis(2, xi1) + &
                            source(x2) * p1_basis(2, xi2))
        
    end subroutine element_load_vector

end module fortfem_assembly_1d
module fortfem_basis_1d
    use fortfem_kinds
    implicit none
    private
    
    public :: p1_basis, p1_basis_derivative
    
contains

    function p1_basis(i, xi) result(phi)
        integer, intent(in) :: i      ! Basis function index (1 or 2)
        real(dp), intent(in) :: xi    ! Local coordinate in [0,1]
        real(dp) :: phi
        
        select case(i)
        case(1)
            phi = 1.0_dp - xi
        case(2)
            phi = xi
        case default
            error stop "Invalid basis function index"
        end select
        
    end function p1_basis
    
    function p1_basis_derivative(i, xi) result(dphi)
        integer, intent(in) :: i      ! Basis function index (1 or 2)
        real(dp), intent(in) :: xi    ! Local coordinate (not used for P1)
        real(dp) :: dphi
        
        associate(dummy => xi)
        end associate
        
        select case(i)
        case(1)
            dphi = -1.0_dp
        case(2)
            dphi = 1.0_dp
        case default
            error stop "Invalid basis function index"
        end select
        
    end function p1_basis_derivative

end module fortfem_basis_1d
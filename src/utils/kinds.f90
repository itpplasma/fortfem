module fortfem_kinds
    implicit none
    private
    
    ! Precision kinds
    integer, parameter, public :: dp = kind(0.0d0)
    
    ! Mathematical constants
    real(dp), parameter, public :: pi = 3.14159265358979323846_dp
    
end module fortfem_kinds
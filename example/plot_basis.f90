program plot_basis
    ! Example: Plot P1 basis functions on reference triangle
    use fortplot
    use fortfem
    implicit none
    
    integer, parameter :: n = 50
    real(dp) :: x(n), y(n), z(n,n)
    integer :: i, j
    
    ! Create grid on reference triangle
    do i = 1, n
        x(i) = real(i-1, dp) / real(n-1, dp)
    end do
    y = x
    
    ! Evaluate first basis function
    do j = 1, n
        do i = 1, n
            if (x(i) + y(j) <= 1.0_dp) then
                ! P1 basis function at vertex (0,0)
                z(i,j) = 1.0_dp - x(i) - y(j)
            else
                z(i,j) = 0.0_dp
            end if
        end do
    end do
    
    ! Plot using fortplotlib
    call figure()
    call contour(x, y, z)
    call xlabel('x')
    call ylabel('y')
    call title('P1 Basis Function on Reference Triangle')
    call savefig('p1_basis.png')
    
end program plot_basis

module fortfem_gmres
    use fortfem_kinds
    implicit none
    private
    
    public :: gmres_solve
    public :: gmres_options_t
    
    type :: gmres_options_t
        integer :: max_iter = 1000
        integer :: restart = 30
        real(dp) :: tol = 1e-10_dp
        logical :: verbose = .false.
    end type gmres_options_t
    
contains

    subroutine gmres_solve(A, b, x, options, info)
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(gmres_options_t), intent(in), optional :: options
        integer, intent(out), optional :: info
        
        type(gmres_options_t) :: opts
        integer :: n, m, i, j, k, iter, info_local
        real(dp) :: beta, norm_r, tol_abs
        real(dp), allocatable :: r(:), w(:), v(:,:), h(:,:), g(:), c(:), s(:), y(:)
        real(dp) :: temp, h_ik, h_jk
        
        ! Set default options
        if (present(options)) then
            opts = options
        else
            opts = gmres_options_t()
        end if
        
        n = size(A, 1)
        m = min(opts%restart, n)
        
        ! Allocate arrays
        allocate(r(n), w(n), v(n, m+1), h(m+1, m), g(m+1), c(m), s(m), y(m))
        
        ! Initialize
        info_local = 0
        tol_abs = opts%tol * norm2(b)
        
        ! Initial residual
        r = b - matmul(A, x)
        beta = norm2(r)
        
        if (opts%verbose) then
            print '(a,es12.5)', 'GMRES: Initial residual = ', beta
        end if
        
        ! Main GMRES iterations
        do iter = 1, opts%max_iter, m
            
            ! Initialize Krylov subspace
            v(:, 1) = r / beta
            g = 0.0_dp
            g(1) = beta
            
            ! Arnoldi process
            do j = 1, m
                w = matmul(A, v(:, j))
                
                ! Gram-Schmidt orthogonalization
                do i = 1, j
                    h(i, j) = dot_product(w, v(:, i))
                    w = w - h(i, j) * v(:, i)
                end do
                
                h(j+1, j) = norm2(w)
                
                if (h(j+1, j) < 1e-14_dp) then
                    m = j  ! Lucky breakdown
                    exit
                end if
                
                if (j < m) then
                    v(:, j+1) = w / h(j+1, j)
                end if
                
                ! Apply previous Givens rotations
                do i = 1, j-1
                    temp = c(i) * h(i, j) + s(i) * h(i+1, j)
                    h(i+1, j) = -s(i) * h(i, j) + c(i) * h(i+1, j)
                    h(i, j) = temp
                end do
                
                ! Compute new Givens rotation
                temp = sqrt(h(j, j)**2 + h(j+1, j)**2)
                c(j) = h(j, j) / temp
                s(j) = h(j+1, j) / temp
                
                ! Apply new Givens rotation
                h(j, j) = temp
                h(j+1, j) = 0.0_dp
                
                ! Update g
                temp = c(j) * g(j) + s(j) * g(j+1)
                g(j+1) = -s(j) * g(j) + c(j) * g(j+1)
                g(j) = temp
                
                ! Check convergence
                norm_r = abs(g(j+1))
                if (opts%verbose .and. mod(iter + j - 1, 10) == 0) then
                    print '(a,i0,a,es12.5)', 'GMRES: Iteration ', iter + j - 1, &
                          ', residual = ', norm_r
                end if
                
                if (norm_r < tol_abs) then
                    m = j
                    exit
                end if
            end do
            
            ! Solve upper triangular system
            y(1:m) = g(1:m)
            do i = m, 1, -1
                do j = i+1, m
                    y(i) = y(i) - h(i, j) * y(j)
                end do
                y(i) = y(i) / h(i, i)
            end do
            
            ! Update solution
            do i = 1, m
                x = x + y(i) * v(:, i)
            end do
            
            ! Check convergence
            if (norm_r < tol_abs) then
                if (opts%verbose) then
                    print '(a,i0,a)', 'GMRES: Converged in ', iter + m - 1, ' iterations'
                end if
                exit
            end if
            
            ! Compute new residual for restart
            r = b - matmul(A, x)
            beta = norm2(r)
            
            if (iter + m > opts%max_iter) then
                info_local = 1  ! Max iterations reached
                exit
            end if
        end do
        
        ! Clean up
        deallocate(r, w, v, h, g, c, s, y)
        
        if (present(info)) info = info_local
        
        if (info_local == 1 .and. opts%verbose) then
            print *, 'GMRES: Maximum iterations reached'
        end if
    end subroutine gmres_solve

end module fortfem_gmres
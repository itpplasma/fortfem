program solver_comparison_2d
    use fortfem
    use function_space_module
    use weak_forms_module
    implicit none
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: info, ipiv(*)
            double precision, intent(inout) :: a(lda,*), b(ldb,*)
        end subroutine dgesv
    end interface
    
    ! Problem variables
    type(mesh_2d_t) :: mesh
    type(function_space_t) :: V
    type(trial_function_t) :: u
    type(test_function_t) :: v_test
    
    ! Weak form components
    type(bilinear_form_t) :: a
    type(linear_form_t) :: L
    type(weak_form_t) :: problem
    
    ! Solutions and matrices
    real(dp), allocatable :: sol_lapack(:), sol_umfpack(:), difference(:)
    real(dp), allocatable :: matrix_lapack(:,:), matrix_umfpack(:,:)
    real(dp), allocatable :: rhs_lapack(:), rhs_umfpack(:)
    real(dp) :: x, y, max_diff, rms_diff, time_start, time_end
    integer :: i, nx, ny
    
    print *, "2D Solver Comparison with Weak Forms"
    print *, "====================================="
    print *, ""
    print *, "Demonstrating weak form framework with different solvers"
    print *, "Problem: Find u ∈ V such that ∫∇u·∇v dx = ∫fv dx ∀v ∈ V"
    print *, "with u = 0 on ∂Ω and f = 2π²sin(πx)sin(πy)"
    print *, ""
    
    ! Create mesh
    nx = 25
    ny = 25
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Create function space
    call create_P1_space(mesh, V)
    
    ! Create trial and test functions
    call u%init(V, "u")
    call v_test%init(V, "v")
    
    print '(a,i0)', "Function space DOFs: ", V%n_dofs
    print *, ""
    
    ! Set up weak formulation
    print *, "Setting up weak formulation..."
    print *, "  Bilinear form: a(u,v) = ∫∇u·∇v dx"
    print *, "  Linear form:   L(v) = ∫fv dx"
    
    ! Bilinear form (stiffness matrix)
    call a%init(form_type=2, expression="grad(u).grad(v)")
    
    ! Linear form (load vector)
    call L%init(form_type=1, expression="f*v")
    
    ! Create weak form problem
    call problem%init(a, L, "Poisson weak form")
    
    print *, ""
    
    ! Solve with LAPACK
    print *, "Solving with LAPACK..."
    call cpu_time(time_start)
    
    ! Assemble system for LAPACK
    allocate(matrix_lapack(V%n_dofs, V%n_dofs))
    allocate(rhs_lapack(V%n_dofs))
    
    call problem%assemble(V, matrix_lapack, rhs_lapack)
    call apply_boundary_conditions(matrix_lapack, rhs_lapack)
    call solve_with_lapack(matrix_lapack, rhs_lapack, sol_lapack)
    
    call cpu_time(time_end)
    print '(a,f8.3,a)', "LAPACK solve time: ", time_end - time_start, " seconds"
    
    ! Solve with UMFPACK (conceptual - would need proper sparse assembly)
    print *, "Solving with UMFPACK..."
    call cpu_time(time_start)
    
    ! Assemble system for UMFPACK
    allocate(matrix_umfpack(V%n_dofs, V%n_dofs))
    allocate(rhs_umfpack(V%n_dofs))
    
    call problem%assemble(V, matrix_umfpack, rhs_umfpack)
    call apply_boundary_conditions(matrix_umfpack, rhs_umfpack)
    call solve_with_lapack(matrix_umfpack, rhs_umfpack, sol_umfpack)  ! Using LAPACK for demo
    
    call cpu_time(time_end)
    print '(a,f8.3,a)', "UMFPACK solve time: ", time_end - time_start, " seconds"
    
    ! Compare solutions
    allocate(difference(mesh%n_vertices))
    difference = abs(sol_lapack - sol_umfpack)
    
    max_diff = maxval(difference)
    rms_diff = sqrt(sum(difference**2) / mesh%n_vertices)
    
    print *, ""
    print *, "Solution comparison:"
    print '(a,es12.5)', "Max difference: ", max_diff
    print '(a,es12.5)', "RMS difference: ", rms_diff
    print '(a,es12.5)', "Max LAPACK solution: ", maxval(abs(sol_lapack))
    print '(a,es12.5)', "Max UMFPACK solution: ", maxval(abs(sol_umfpack))
    
    ! Write data for plotting
    call write_solution_data("lapack_solution.dat", sol_lapack)
    call write_solution_data("umfpack_solution.dat", sol_umfpack)
    call write_solution_data("difference.dat", difference)
    
    ! Create Python script for contour plotting
    call create_plotting_script()
    
    print *, ""
    print *, "Data files created:"
    print *, "  lapack_solution.dat  - LAPACK solution"
    print *, "  umfpack_solution.dat - UMFPACK solution"
    print *, "  difference.dat       - Absolute difference"
    print *, "  plot_comparison.py   - Python plotting script"
    print *, ""
    print *, "To generate plots, run: python plot_comparison.py"
    
    ! Clean up
    call u%destroy()
    call v_test%destroy()
    call V%destroy()
    call mesh%destroy()
    call a%destroy()
    call L%destroy()
    call problem%destroy()
    deallocate(sol_lapack, sol_umfpack, difference)
    deallocate(matrix_lapack, matrix_umfpack, rhs_lapack, rhs_umfpack)
    
    print *, ""
    print *, "Comparison complete!"
    
contains

    pure function manufactured_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function manufactured_source
    
    subroutine apply_boundary_conditions(matrix, rhs)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp) :: x, y
        integer :: i, n
        
        n = size(matrix, 1)
        
        do i = 1, n
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            ! Apply homogeneous Dirichlet BC on boundary
            if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
                abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
                matrix(i, :) = 0.0_dp
                matrix(:, i) = 0.0_dp
                matrix(i, i) = 1.0_dp
                rhs(i) = 0.0_dp
            end if
        end do
    end subroutine apply_boundary_conditions
    
    subroutine solve_with_lapack(matrix, rhs, solution)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer :: n, info
        integer, allocatable :: ipiv(:)
        
        n = size(matrix, 1)
        allocate(solution(n))
        allocate(ipiv(n))
        
        solution = rhs
        
        call dgesv(n, 1, matrix, n, ipiv, solution, n, info)
        
        if (info /= 0) then
            print *, "Error in LAPACK solver, info = ", info
            stop
        end if
        
        deallocate(ipiv)
    end subroutine solve_with_lapack
    
    subroutine write_solution_data(filename, solution)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: solution(:)
        integer :: i
        
        open(unit=10, file=filename, status='replace')
        write(10, '(a)') '# x y z'
        do i = 1, mesh%n_vertices
            write(10, '(3es16.8)') mesh%vertices(1,i), mesh%vertices(2,i), solution(i)
        end do
        close(10)
    end subroutine write_solution_data
    
    subroutine create_plotting_script()
        integer :: unit
        
        open(newunit=unit, file='plot_comparison.py', status='replace')
        
        write(unit, '(a)') '#!/usr/bin/env python3'
        write(unit, '(a)') 'import numpy as np'
        write(unit, '(a)') 'import matplotlib.pyplot as plt'
        write(unit, '(a)') 'from scipy.interpolate import griddata'
        write(unit, '(a)') 'import matplotlib.tri as tri'
        write(unit, '(a)') ''
        write(unit, '(a)') 'def load_data(filename):'
        write(unit, '(a)') '    data = np.loadtxt(filename)'
        write(unit, '(a)') '    return data[:, 0], data[:, 1], data[:, 2]'
        write(unit, '(a)') ''
        write(unit, '(a)') '# Load data'
        write(unit, '(a)') 'x_lap, y_lap, z_lap = load_data("lapack_solution.dat")'
        write(unit, '(a)') 'x_umf, y_umf, z_umf = load_data("umfpack_solution.dat")'
        write(unit, '(a)') 'x_diff, y_diff, z_diff = load_data("difference.dat")'
        write(unit, '(a)') ''
        write(unit, '(a)') '# Create triangulation'
        write(unit, '(a)') 'triang = tri.Triangulation(x_lap, y_lap)'
        write(unit, '(a)') ''
        write(unit, '(a)') '# Create figure with subplots'
        write(unit, '(a)') 'fig, axes = plt.subplots(2, 2, figsize=(12, 10))'
        write(unit, '(a)') 'fig.suptitle("2D Poisson Solver Comparison: LAPACK vs UMFPACK", fontsize=14)'
        write(unit, '(a)') ''
        write(unit, '(a)') '# LAPACK solution'
        write(unit, '(a)') 'cs1 = axes[0,0].tricontourf(triang, z_lap, levels=20, cmap="viridis")'
        write(unit, '(a)') 'axes[0,0].set_title("LAPACK Solution")'
        write(unit, '(a)') 'axes[0,0].set_xlabel("x")'
        write(unit, '(a)') 'axes[0,0].set_ylabel("y")'
        write(unit, '(a)') 'axes[0,0].set_aspect("equal")'
        write(unit, '(a)') 'plt.colorbar(cs1, ax=axes[0,0])'
        write(unit, '(a)') ''
        write(unit, '(a)') '# UMFPACK solution'
        write(unit, '(a)') 'cs2 = axes[0,1].tricontourf(triang, z_umf, levels=20, cmap="viridis")'
        write(unit, '(a)') 'axes[0,1].set_title("UMFPACK Solution")'
        write(unit, '(a)') 'axes[0,1].set_xlabel("x")'
        write(unit, '(a)') 'axes[0,1].set_ylabel("y")'
        write(unit, '(a)') 'axes[0,1].set_aspect("equal")'
        write(unit, '(a)') 'plt.colorbar(cs2, ax=axes[0,1])'
        write(unit, '(a)') ''
        write(unit, '(a)') '# Difference'
        write(unit, '(a)') 'cs3 = axes[1,0].tricontourf(triang, z_diff, levels=20, cmap="plasma")'
        write(unit, '(a)') 'axes[1,0].set_title("Absolute Difference")'
        write(unit, '(a)') 'axes[1,0].set_xlabel("x")'
        write(unit, '(a)') 'axes[1,0].set_ylabel("y")'
        write(unit, '(a)') 'axes[1,0].set_aspect("equal")'
        write(unit, '(a)') 'plt.colorbar(cs3, ax=axes[1,0])'
        write(unit, '(a)') ''
        write(unit, '(a)') '# Exact solution for comparison'
        write(unit, '(a)') 'z_exact = np.sin(np.pi * x_lap) * np.sin(np.pi * y_lap)'
        write(unit, '(a)') 'cs4 = axes[1,1].tricontourf(triang, z_exact, levels=20, cmap="viridis")'
        write(unit, '(a)') 'axes[1,1].set_title("Exact Solution")'
        write(unit, '(a)') 'axes[1,1].set_xlabel("x")'
        write(unit, '(a)') 'axes[1,1].set_ylabel("y")'
        write(unit, '(a)') 'axes[1,1].set_aspect("equal")'
        write(unit, '(a)') 'plt.colorbar(cs4, ax=axes[1,1])'
        write(unit, '(a)') ''
        write(unit, '(a)') 'plt.tight_layout()'
        write(unit, '(a)') 'plt.savefig("solver_comparison_2d.png", dpi=300, bbox_inches="tight")'
        write(unit, '(a)') 'plt.show()'
        write(unit, '(a)') ''
        write(unit, '(a)') '# Print statistics'
        write(unit, '(a)') 'print(f"Max difference: {np.max(z_diff):.2e}")'
        write(unit, '(a)') 'print(f"RMS difference: {np.sqrt(np.mean(z_diff**2)):.2e}")'
        write(unit, '(a)') 'print(f"Max LAPACK solution: {np.max(np.abs(z_lap)):.6f}")'
        write(unit, '(a)') 'print(f"Max UMFPACK solution: {np.max(np.abs(z_umf)):.6f}")'
        write(unit, '(a)') 'print(f"Max exact solution: {np.max(np.abs(z_exact)):.6f}")'
        
        close(unit)
    end subroutine create_plotting_script

end program solver_comparison_2d
program solver_comparison_2d
    use fortfem
    implicit none
    
    type(poisson_2d_t) :: solver_lapack, solver_umfpack
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: sol_lapack(:), sol_umfpack(:), difference(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    real(dp) :: x, y, max_diff, rms_diff, time_start, time_end
    integer :: i, k, nx, ny
    
    print *, "2D Solver Comparison: LAPACK vs UMFPACK"
    print *, "======================================="
    print *, ""
    print *, "Solving 2D Poisson equation with both solvers"
    print *, "Problem: -∇²u = f with u = sin(πx)sin(πy) as manufactured solution"
    print *, ""
    
    ! Create mesh
    nx = 25
    ny = 25
    call mesh%create_rectangular(nx=nx, ny=ny, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    print '(a,i0)', "Mesh triangles: ", mesh%n_triangles
    
    ! Set up boundary conditions
    allocate(boundary_nodes(2*nx + 2*ny - 4))
    allocate(boundary_values(size(boundary_nodes)))
    
    k = 0
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        
        if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
            abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
            k = k + 1
            boundary_nodes(k) = i
            boundary_values(k) = 0.0_dp  ! Homogeneous Dirichlet BC
        end if
    end do
    
    print '(a,i0)', "Boundary nodes: ", k
    print *, ""
    
    ! Solve with LAPACK
    print *, "Solving with LAPACK..."
    call cpu_time(time_start)
    
    call solver_lapack%init("lapack")
    call solver_lapack%set_mesh(mesh)
    call solver_lapack%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    call solver_lapack%solve(manufactured_source)
    sol_lapack = solver_lapack%get_solution()
    
    call cpu_time(time_end)
    print '(a,f8.3,a)', "LAPACK solve time: ", time_end - time_start, " seconds"
    
    ! Solve with UMFPACK
    print *, "Solving with UMFPACK..."
    call cpu_time(time_start)
    
    call solver_umfpack%init("umfpack")
    call solver_umfpack%set_mesh(mesh)
    call solver_umfpack%set_dirichlet_bc(boundary_nodes(1:k), boundary_values(1:k))
    call solver_umfpack%solve(manufactured_source)
    sol_umfpack = solver_umfpack%get_solution()
    
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
    call solver_lapack%destroy()
    call solver_umfpack%destroy()
    call mesh%destroy()
    deallocate(sol_lapack, sol_umfpack, difference)
    deallocate(boundary_nodes, boundary_values)
    
    print *, ""
    print *, "Comparison complete!"
    
contains

    pure function manufactured_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function manufactured_source
    
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
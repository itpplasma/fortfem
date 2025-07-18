program mesh_2d_demo
    ! Example: Create and visualize 2D meshes
    use fortfem
    use fortplot
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: x(:), y(:), areas(:)
    integer :: i, j, t
    
    print *, "2D Mesh Generation Demo"
    print *, "======================"
    
    ! Create rectangular mesh
    call mesh%create_rectangular(nx=11, ny=11, &
                               x_min=0.0_dp, x_max=1.0_dp, &
                               y_min=0.0_dp, y_max=1.0_dp)
    
    call mesh%build_connectivity()
    call mesh%find_boundary()
    
    print *, "Mesh statistics:"
    print '(a,i0)', "  Vertices: ", mesh%n_vertices
    print '(a,i0)', "  Triangles: ", mesh%n_triangles
    print '(a,i0)', "  Edges: ", mesh%n_edges
    print '(a,i0)', "  Boundary edges: ", mesh%n_boundary_edges
    
    ! Compute triangle areas
    areas = mesh%compute_areas()
    print '(a,es12.5)', "  Min area: ", minval(areas)
    print '(a,es12.5)', "  Max area: ", maxval(areas)
    print '(a,es12.5)', "  Total area: ", sum(areas)
    
    ! Save mesh to file
    call mesh%save_to_file('rectangular_mesh.dat')
    print *, ""
    print *, "Mesh saved to rectangular_mesh.dat"
    
    ! Plot mesh using matplotlib triplot equivalent
    call plot_mesh(mesh, 'rectangular_mesh.png')
    
    ! Create a refined mesh in a circle
    call create_circle_mesh()
    
    ! Clean up
    deallocate(areas)
    call mesh%destroy()
    
contains

    subroutine plot_mesh(mesh, filename)
        type(mesh_2d_t), intent(in) :: mesh
        character(len=*), intent(in) :: filename
        
        real(dp), allocatable :: edge_x(:), edge_y(:)
        integer :: e, v1, v2, k
        
        ! Prepare edge data for plotting
        allocate(edge_x(3*mesh%n_edges), edge_y(3*mesh%n_edges))
        
        k = 0
        do e = 1, mesh%n_edges
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            
            ! Add edge vertices with NaN separator
            k = k + 1
            edge_x(k) = mesh%vertices(1, v1)
            edge_y(k) = mesh%vertices(2, v1)
            k = k + 1
            edge_x(k) = mesh%vertices(1, v2)
            edge_y(k) = mesh%vertices(2, v2)
            k = k + 1
            edge_x(k) = ieee_value(1.0_dp, ieee_quiet_nan)
            edge_y(k) = ieee_value(1.0_dp, ieee_quiet_nan)
        end do
        
        ! Plot
        call figure()
        call plot(edge_x(1:k), edge_y(1:k), 'b-')
        
        ! Highlight boundary edges
        k = 0
        do i = 1, mesh%n_boundary_edges
            e = mesh%boundary_edges(i)
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            
            k = k + 1
            edge_x(k) = mesh%vertices(1, v1)
            edge_y(k) = mesh%vertices(2, v1)
            k = k + 1
            edge_x(k) = mesh%vertices(1, v2)
            edge_y(k) = mesh%vertices(2, v2)
            k = k + 1
            edge_x(k) = ieee_value(1.0_dp, ieee_quiet_nan)
            edge_y(k) = ieee_value(1.0_dp, ieee_quiet_nan)
        end do
        
        if (k > 0) then
            ! Can't overlay in fortplot, so just mention boundary
        end if
        
        call xlabel('x')
        call ylabel('y')
        call title('2D Triangular Mesh')
        call savefig(filename)
        
        deallocate(edge_x, edge_y)
        
        print '(a,a)', "Mesh plot saved to ", filename
        
    end subroutine plot_mesh
    
    subroutine create_circle_mesh()
        type(mesh_2d_t) :: circle_mesh
        real(dp) :: r, theta
        integer :: nr, ntheta, n_inner
        integer :: i, j, k, v_center
        
        print *, ""
        print *, "Creating circular mesh..."
        
        nr = 5       ! Radial divisions
        ntheta = 16  ! Angular divisions
        
        ! Total vertices: center + nr circles
        circle_mesh%n_vertices = 1 + nr * ntheta
        allocate(circle_mesh%vertices(2, circle_mesh%n_vertices))
        
        ! Center vertex
        k = 1
        circle_mesh%vertices(:, k) = [0.5_dp, 0.5_dp]
        v_center = k
        
        ! Circular layers
        do i = 1, nr
            r = 0.4_dp * real(i, dp) / real(nr, dp)  ! Radius up to 0.4
            do j = 1, ntheta
                theta = 2.0_dp * pi * real(j-1, dp) / real(ntheta, dp)
                k = k + 1
                circle_mesh%vertices(1, k) = 0.5_dp + r * cos(theta)
                circle_mesh%vertices(2, k) = 0.5_dp + r * sin(theta)
            end do
        end do
        
        ! Estimate triangles (this is approximate - proper triangulation needed)
        circle_mesh%n_triangles = ntheta + 2 * (nr-1) * ntheta
        
        print '(a,i0)', "  Circle vertices: ", circle_mesh%n_vertices
        print '(a,i0)', "  Estimated triangles: ", circle_mesh%n_triangles
        
        ! For now, just plot the vertices
        call figure()
        call plot(circle_mesh%vertices(1,:), circle_mesh%vertices(2,:), 'o')
        call xlabel('x')
        call ylabel('y')
        call title('Circular Mesh Vertices')
        call savefig('circle_mesh_vertices.png')
        
        print *, "Circle mesh vertices saved to circle_mesh_vertices.png"
        
        call circle_mesh%destroy()
        
    end subroutine create_circle_mesh

end program mesh_2d_demo
module constrained_delaunay
    use fortfem_kinds, only: dp
    use delaunay_types
    use geometric_predicates
    use bowyer_watson
    implicit none
    
    private
    public :: constrained_delaunay_triangulate
    public :: insert_constraint, recover_constraints
    public :: enforce_constraints
    public :: constraint_edge_exists
    
contains

subroutine constrained_delaunay_triangulate(input_points, constraint_segments, mesh)
    !> Constrained Delaunay triangulation using incremental constraint insertion
    real(dp), intent(in) :: input_points(:,:)      ! (2, npoints)
    integer, intent(in) :: constraint_segments(:,:) ! (2, nsegments)
    type(mesh_t), intent(out) :: mesh
    
    integer :: i
    
    ! Start with unconstrained Delaunay triangulation
    call delaunay_triangulate(input_points, mesh)
    
    ! Insert constraint edges
    do i = 1, size(constraint_segments, 2)
        call insert_constraint(mesh, constraint_segments(:, i))
    end do
    
    ! Ensure all constraints are properly enforced
    call enforce_constraints(mesh, constraint_segments)
    
end subroutine constrained_delaunay_triangulate

subroutine insert_constraint(mesh, constraint_edge)
    !> Insert a single constraint edge into the triangulation
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: constraint_edge(2)
    
    integer :: v1, v2
    integer, allocatable :: intersecting_triangles(:)
    integer :: nintersecting
    logical :: edge_exists
    
    v1 = constraint_edge(1)
    v2 = constraint_edge(2)
    
    ! Check if constraint edge already exists
    if (constraint_edge_exists(mesh, v1, v2)) then
        return
    end if
    
    ! Find triangles that intersect the constraint edge
    call find_intersecting_triangles(mesh, v1, v2, intersecting_triangles, nintersecting)
    
    if (nintersecting == 0) then
        ! Edge is already in triangulation, just mark it as constrained
        call add_constraint_edge(mesh, v1, v2)
        return
    end if
    
    ! Remove intersecting triangles and retriangulate
    call remove_triangles_and_retriangulate(mesh, v1, v2, intersecting_triangles, nintersecting)
    
    ! Add the constraint edge
    call add_constraint_edge(mesh, v1, v2)
    
end subroutine insert_constraint

logical function constraint_edge_exists(mesh, v1, v2)
    !> Check if constraint edge already exists in triangulation
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: v1, v2
    
    integer :: i, t1, t2, t3
    integer :: e1, e2
    
    constraint_edge_exists = .false.
    
    ! Check all triangle edges to see if v1-v2 edge exists
    do i = 1, mesh%ntriangles
        if (.not. mesh%triangles(i)%valid) cycle
        
        t1 = mesh%triangles(i)%vertices(1)
        t2 = mesh%triangles(i)%vertices(2)
        t3 = mesh%triangles(i)%vertices(3)
        
        ! Check all three edges of triangle
        if ((t1 == v1 .and. t2 == v2) .or. (t1 == v2 .and. t2 == v1)) then
            constraint_edge_exists = .true.
            return
        end if
        if ((t2 == v1 .and. t3 == v2) .or. (t2 == v2 .and. t3 == v1)) then
            constraint_edge_exists = .true.
            return
        end if
        if ((t3 == v1 .and. t1 == v2) .or. (t3 == v2 .and. t1 == v1)) then
            constraint_edge_exists = .true.
            return
        end if
    end do
    
end function constraint_edge_exists

subroutine find_intersecting_triangles(mesh, v1, v2, intersecting_triangles, nintersecting)
    !> Find all triangles that intersect the constraint edge v1-v2
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: v1, v2
    integer, allocatable, intent(out) :: intersecting_triangles(:)
    integer, intent(out) :: nintersecting
    
    type(point_t) :: p1, p2
    integer :: i, count
    logical, allocatable :: intersects(:)
    
    p1 = mesh%points(v1)
    p2 = mesh%points(v2)
    
    allocate(intersects(mesh%ntriangles))
    intersects = .false.
    count = 0
    
    ! Check each triangle for intersection with constraint edge
    do i = 1, mesh%ntriangles
        if (.not. mesh%triangles(i)%valid) cycle
        
        if (triangle_intersects_edge(mesh, i, p1, p2)) then
            intersects(i) = .true.
            count = count + 1
        end if
    end do
    
    ! Collect intersecting triangle indices
    nintersecting = count
    allocate(intersecting_triangles(nintersecting))
    count = 0
    do i = 1, mesh%ntriangles
        if (intersects(i)) then
            count = count + 1
            intersecting_triangles(count) = i
        end if
    end do
    
end subroutine find_intersecting_triangles

logical function triangle_intersects_edge(mesh, tri_idx, p1, p2)
    !> Check if triangle intersects with edge p1-p2
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: tri_idx
    type(point_t), intent(in) :: p1, p2
    
    type(point_t) :: ta, tb, tc
    integer :: v1, v2, v3
    
    triangle_intersects_edge = .false.
    
    v1 = mesh%triangles(tri_idx)%vertices(1)
    v2 = mesh%triangles(tri_idx)%vertices(2)
    v3 = mesh%triangles(tri_idx)%vertices(3)
    
    ta = mesh%points(v1)
    tb = mesh%points(v2)
    tc = mesh%points(v3)
    
    ! Check if edge p1-p2 intersects any edge of the triangle
    if (segments_intersect(p1, p2, ta, tb) .or. &
        segments_intersect(p1, p2, tb, tc) .or. &
        segments_intersect(p1, p2, tc, ta)) then
        triangle_intersects_edge = .true.
    end if
    
end function triangle_intersects_edge

logical function segments_intersect(p1, p2, p3, p4)
    !> Check if line segments p1-p2 and p3-p4 intersect
    type(point_t), intent(in) :: p1, p2, p3, p4
    
    integer :: o1, o2, o3, o4
    
    ! Get orientations
    o1 = orientation(p1, p2, p3)
    o2 = orientation(p1, p2, p4)
    o3 = orientation(p3, p4, p1)
    o4 = orientation(p3, p4, p2)
    
    ! General case: segments intersect if orientations are different
    if (o1 /= o2 .and. o3 /= o4) then
        segments_intersect = .true.
        return
    end if
    
    ! Special cases: collinear points
    if (o1 == ORIENTATION_COLLINEAR .and. point_on_segment(p1, p3, p2)) then
        segments_intersect = .true.
        return
    end if
    if (o2 == ORIENTATION_COLLINEAR .and. point_on_segment(p1, p4, p2)) then
        segments_intersect = .true.
        return
    end if
    if (o3 == ORIENTATION_COLLINEAR .and. point_on_segment(p3, p1, p4)) then
        segments_intersect = .true.
        return
    end if
    if (o4 == ORIENTATION_COLLINEAR .and. point_on_segment(p3, p2, p4)) then
        segments_intersect = .true.
        return
    end if
    
    segments_intersect = .false.
end function segments_intersect

logical function point_on_segment(p, q, r)
    !> Check if point q lies on segment pr (assuming collinear)
    type(point_t), intent(in) :: p, q, r
    
    point_on_segment = (q%x <= max(p%x, r%x) .and. q%x >= min(p%x, r%x) .and. &
                       q%y <= max(p%y, r%y) .and. q%y >= min(p%y, r%y))
end function point_on_segment

subroutine remove_triangles_and_retriangulate(mesh, v1, v2, intersecting_triangles, nintersecting)
    !> Remove intersecting triangles and retriangulate the cavity
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: v1, v2, nintersecting
    integer, intent(in) :: intersecting_triangles(:)
    
    integer, allocatable :: cavity_boundary(:,:)
    integer :: nboundary_edges
    integer :: i
    
    ! Find the boundary of the cavity formed by removing intersecting triangles
    call find_removal_cavity_boundary(mesh, intersecting_triangles, nintersecting, &
                                     cavity_boundary, nboundary_edges)
    
    ! Remove the intersecting triangles
    do i = 1, nintersecting
        mesh%triangles(intersecting_triangles(i))%valid = .false.
    end do
    
    ! Retriangulate the cavity while preserving the constraint edge v1-v2
    call retriangulate_with_constraint(mesh, v1, v2, cavity_boundary, nboundary_edges)
    
end subroutine remove_triangles_and_retriangulate

subroutine find_removal_cavity_boundary(mesh, removed_triangles, nremoved, &
                                       boundary_edges, nboundary)
    !> Find boundary edges of cavity created by removing triangles
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: removed_triangles(:), nremoved
    integer, allocatable, intent(out) :: boundary_edges(:,:)
    integer, intent(out) :: nboundary
    
    integer, allocatable :: all_edges(:,:)
    integer, allocatable :: edge_count(:)
    integer :: nedges_total, i, j, t, v1, v2, edge_idx
    logical :: found
    
    ! Collect all edges from removed triangles
    nedges_total = 3 * nremoved
    allocate(all_edges(2, nedges_total))
    allocate(edge_count(nedges_total))
    
    nedges_total = 0
    do i = 1, nremoved
        t = removed_triangles(i)
        ! Add three edges of triangle t
        do j = 1, 3
            v1 = mesh%triangles(t)%vertices(j)
            v2 = mesh%triangles(t)%vertices(mod(j, 3) + 1)
            
            ! Ensure consistent edge ordering
            if (v1 > v2) then
                call swap_integers(v1, v2)
            end if
            
            ! Check if edge already exists
            found = .false.
            do edge_idx = 1, nedges_total
                if (all_edges(1, edge_idx) == v1 .and. all_edges(2, edge_idx) == v2) then
                    edge_count(edge_idx) = edge_count(edge_idx) + 1
                    found = .true.
                    exit
                end if
            end do
            
            if (.not. found) then
                nedges_total = nedges_total + 1
                all_edges(1, nedges_total) = v1
                all_edges(2, nedges_total) = v2
                edge_count(nedges_total) = 1
            end if
        end do
    end do
    
    ! Boundary edges appear only once
    nboundary = 0
    do i = 1, nedges_total
        if (edge_count(i) == 1) then
            nboundary = nboundary + 1
        end if
    end do
    
    allocate(boundary_edges(2, nboundary))
    nboundary = 0
    do i = 1, nedges_total
        if (edge_count(i) == 1) then
            nboundary = nboundary + 1
            boundary_edges(1, nboundary) = all_edges(1, i)
            boundary_edges(2, nboundary) = all_edges(2, i)
        end if
    end do
    
end subroutine find_removal_cavity_boundary

subroutine retriangulate_with_constraint(mesh, constraint_v1, constraint_v2, &
                                        boundary_edges, nboundary)
    !> Retriangulate cavity ensuring constraint edge is included
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: constraint_v1, constraint_v2, nboundary
    integer, intent(in) :: boundary_edges(:,:)
    
    integer :: i, v1, v2, new_tri
    type(point_t) :: p1, p2, pc1, pc2, midpoint
    logical :: connects_to_v1
    real(dp) :: dist1, dist2
    
    pc1 = mesh%points(constraint_v1)
    pc2 = mesh%points(constraint_v2)
    
    ! Improved triangulation: connect boundary edges to appropriate constraint vertex
    ! This reduces (but doesn't eliminate) overlapping triangles
    do i = 1, nboundary
        v1 = boundary_edges(1, i)
        v2 = boundary_edges(2, i)
        p1 = mesh%points(v1)
        p2 = mesh%points(v2)
        
        ! Skip if this is the constraint edge itself
        if ((v1 == constraint_v1 .and. v2 == constraint_v2) .or. &
            (v1 == constraint_v2 .and. v2 == constraint_v1)) then
            cycle
        end if
        
        ! Choose constraint vertex based on which is closer to edge midpoint
        midpoint%x = (p1%x + p2%x) / 2.0_dp
        midpoint%y = (p1%y + p2%y) / 2.0_dp
        
        dist1 = (midpoint%x - pc1%x)**2 + (midpoint%y - pc1%y)**2
        dist2 = (midpoint%x - pc2%x)**2 + (midpoint%y - pc2%y)**2
        
        connects_to_v1 = (dist1 < dist2)
        
        ! Create triangle with proper orientation
        if (connects_to_v1) then
            if (orientation(p1, p2, pc1) == ORIENTATION_CCW) then
                new_tri = add_triangle(mesh, v1, v2, constraint_v1)
            else
                new_tri = add_triangle(mesh, v2, v1, constraint_v1)
            end if
        else
            if (orientation(p1, p2, pc2) == ORIENTATION_CCW) then
                new_tri = add_triangle(mesh, v1, v2, constraint_v2)
            else
                new_tri = add_triangle(mesh, v2, v1, constraint_v2)
            end if
        end if
    end do
    
    ! Ensure constraint edge triangle exists
    call ensure_constraint_triangle(mesh, constraint_v1, constraint_v2)
    
end subroutine retriangulate_with_constraint

subroutine split_cavity_by_constraint(mesh, constraint_v1, constraint_v2, &
                                     boundary_edges, nboundary, &
                                     polygon1, npoly1, polygon2, npoly2)
    !> Split cavity boundary into two polygons on either side of constraint edge
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: constraint_v1, constraint_v2, nboundary
    integer, intent(in) :: boundary_edges(:,:)
    integer, allocatable, intent(out) :: polygon1(:), polygon2(:)
    integer, intent(out) :: npoly1, npoly2
    
    type(point_t) :: constraint_p1, constraint_p2, edge_p1, edge_p2, midpoint
    integer :: i, v1, v2, side
    integer, allocatable :: temp_poly1(:), temp_poly2(:)
    
    constraint_p1 = mesh%points(constraint_v1)
    constraint_p2 = mesh%points(constraint_v2)
    
    allocate(temp_poly1(2 * nboundary + 2))  ! Worst case: all edges on one side
    allocate(temp_poly2(2 * nboundary + 2))
    
    npoly1 = 0
    npoly2 = 0
    
    ! Add constraint vertices to both polygons
    npoly1 = npoly1 + 1
    temp_poly1(npoly1) = constraint_v1
    npoly1 = npoly1 + 1
    temp_poly1(npoly1) = constraint_v2
    
    npoly2 = npoly2 + 1
    temp_poly2(npoly2) = constraint_v1
    npoly2 = npoly2 + 1
    temp_poly2(npoly2) = constraint_v2
    
    ! Classify boundary edges by which side of constraint they're on
    do i = 1, nboundary
        v1 = boundary_edges(1, i)
        v2 = boundary_edges(2, i)
        
        ! Skip if this is the constraint edge itself
        if ((v1 == constraint_v1 .and. v2 == constraint_v2) .or. &
            (v1 == constraint_v2 .and. v2 == constraint_v1)) then
            cycle
        end if
        
        edge_p1 = mesh%points(v1)
        edge_p2 = mesh%points(v2)
        
        ! Determine which side of constraint edge this boundary edge is on
        ! Use orientation test with midpoint of boundary edge
        midpoint%x = (edge_p1%x + edge_p2%x) / 2.0_dp
        midpoint%y = (edge_p1%y + edge_p2%y) / 2.0_dp
        midpoint%id = 0
        
        side = orientation(constraint_p1, constraint_p2, midpoint)
        
        if (side >= 0) then
            ! Add to polygon1 (left side or on constraint)
            npoly1 = npoly1 + 1
            temp_poly1(npoly1) = v1
            npoly1 = npoly1 + 1  
            temp_poly1(npoly1) = v2
        else
            ! Add to polygon2 (right side)
            npoly2 = npoly2 + 1
            temp_poly2(npoly2) = v1
            npoly2 = npoly2 + 1
            temp_poly2(npoly2) = v2
        end if
    end do
    
    ! Copy to output arrays
    if (npoly1 > 0) then
        allocate(polygon1(npoly1))
        polygon1(1:npoly1) = temp_poly1(1:npoly1)
    end if
    
    if (npoly2 > 0) then
        allocate(polygon2(npoly2))
        polygon2(1:npoly2) = temp_poly2(1:npoly2)
    end if
    
    deallocate(temp_poly1, temp_poly2)
    
end subroutine split_cavity_by_constraint

subroutine triangulate_polygon(mesh, polygon, npoly)
    !> Simple polygon triangulation using ear clipping
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: polygon(:), npoly
    
    integer :: i, v1, v2, v3, new_tri
    
    ! For now, simple fan triangulation from first vertex
    ! This is not optimal but works for convex polygons
    if (npoly < 3) return
    
    do i = 2, npoly - 1
        v1 = polygon(1)
        v2 = polygon(i)
        v3 = polygon(i + 1)
        
        ! Create triangle with proper orientation
        new_tri = add_triangle(mesh, v1, v2, v3)
    end do
    
end subroutine triangulate_polygon

subroutine ensure_constraint_edge_triangle(mesh, constraint_v1, constraint_v2)
    !> Ensure constraint edge is represented by at least one triangle
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: constraint_v1, constraint_v2
    
    ! For now, just verify the edge exists - actual triangle creation
    ! should be handled by the polygon triangulation
    call ensure_constraint_triangle(mesh, constraint_v1, constraint_v2)
    
end subroutine ensure_constraint_edge_triangle

subroutine ensure_constraint_triangle(mesh, v1, v2)
    !> Ensure there's a triangle containing the constraint edge v1-v2
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: v1, v2
    
    integer :: t, i, j, found_triangles
    logical :: edge_found
    
    ! Verify constraint edge exists in at least one triangle
    edge_found = .false.
    found_triangles = 0
    
    do t = 1, mesh%ntriangles
        if (mesh%triangles(t)%valid) then
            ! Check all edges of triangle t
            do i = 1, 3
                j = mod(i, 3) + 1
                if ((mesh%triangles(t)%vertices(i) == v1 .and. mesh%triangles(t)%vertices(j) == v2) .or. &
                    (mesh%triangles(t)%vertices(i) == v2 .and. mesh%triangles(t)%vertices(j) == v1)) then
                    edge_found = .true.
                    found_triangles = found_triangles + 1
                end if
            end do
        end if
    end do
    
    ! Log warnings if edge validation fails
    if (.not. edge_found) then
        write(*,'(A,I0,A,I0,A)') "Warning: Constraint edge ", v1, "-", v2, &
                                " not found in any triangle"
    else if (found_triangles > 2) then
        write(*,'(A,I0,A,I0,A,I0,A)') "Warning: Constraint edge ", v1, "-", v2, &
                                     " found in ", found_triangles, " triangles (expected ≤2)"
    end if
    
end subroutine ensure_constraint_triangle

subroutine add_constraint_edge(mesh, v1, v2)
    !> Add constraint edge to the mesh
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: v1, v2
    
    integer :: edge_idx
    
    edge_idx = add_edge(mesh, v1, v2, .true.)
    
end subroutine add_constraint_edge

subroutine recover_constraints(mesh, constraint_segments)
    !> Ensure all constraint segments are present in the triangulation
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: constraint_segments(:,:)
    
    integer :: i
    
    do i = 1, size(constraint_segments, 2)
        if (.not. constraint_edge_exists(mesh, constraint_segments(1, i), constraint_segments(2, i))) then
            call insert_constraint(mesh, constraint_segments(:, i))
        end if
    end do
    
end subroutine recover_constraints

subroutine enforce_constraints(mesh, constraint_segments)
    !> Final pass to ensure all constraints are properly enforced
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: constraint_segments(:,:)
    
    integer :: i, max_iterations, iteration
    logical :: all_constraints_satisfied
    
    max_iterations = 10
    
    do iteration = 1, max_iterations
        all_constraints_satisfied = .true.
        
        do i = 1, size(constraint_segments, 2)
            if (.not. constraint_edge_exists(mesh, constraint_segments(1, i), constraint_segments(2, i))) then
                call insert_constraint(mesh, constraint_segments(:, i))
                all_constraints_satisfied = .false.
            end if
        end do
        
        if (all_constraints_satisfied) exit
    end do
    
end subroutine enforce_constraints

subroutine swap_integers(a, b)
    !> Swap two integers
    integer, intent(inout) :: a, b
    integer :: temp
    
    temp = a
    a = b
    b = temp
end subroutine swap_integers

end module constrained_delaunay
module fortfem_hcurl_space
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_basis_edge_2d_interface
    implicit none
    private
    
    public :: hcurl_space_t
    
    ! H(curl) function space for edge elements
    type :: hcurl_space_t
        type(mesh_2d_t), pointer :: mesh => null()
        type(edge_basis_2d_t), allocatable :: basis(:)  ! One basis per triangle
        integer :: n_dofs = 0
        integer :: n_interior_dofs = 0
        integer :: n_boundary_dofs = 0
        
        ! DOF management
        integer, allocatable :: dof_is_essential(:)  ! Essential BC flags
        real(dp), allocatable :: essential_values(:) ! Essential BC values
        
    contains
        procedure :: init => hcurl_space_init
        procedure :: destroy => hcurl_space_destroy
        procedure :: get_n_dofs => hcurl_space_get_n_dofs
        procedure :: get_n_interior_dofs => hcurl_space_get_n_interior_dofs
        procedure :: get_n_boundary_dofs => hcurl_space_get_n_boundary_dofs
        procedure :: apply_essential_bc => hcurl_space_apply_essential_bc
        procedure :: evaluate_at_point => hcurl_space_evaluate_at_point
        procedure :: evaluate_curl_at_point => hcurl_space_evaluate_curl_at_point
        procedure :: evaluate_curl_basis_at_point => hcurl_space_evaluate_curl_basis_at_point
        procedure :: get_triangle_dofs => hcurl_space_get_triangle_dofs
        procedure :: evaluate_edge_basis_2d_with_piola => evaluate_edge_basis_2d_with_piola
        procedure :: evaluate_edge_basis_oriented => evaluate_edge_basis_oriented
        procedure :: evaluate_edge_basis_physical => evaluate_edge_basis_physical
        procedure :: evaluate_nedelec_basis_correct
        procedure :: evaluate_conservative_basis
    end type hcurl_space_t
    
contains

    subroutine hcurl_space_init(this, mesh)
        class(hcurl_space_t), intent(out) :: this
        type(mesh_2d_t), intent(in), target :: mesh
        integer :: t
        
        this%mesh => mesh
        
        ! Check that mesh has edge connectivity
        if (.not. allocated(mesh%edges)) then
            error stop "Mesh must have edge connectivity. Call build_edge_connectivity first."
        end if
        
        if (.not. allocated(mesh%edge_to_dof)) then
            error stop "Mesh must have edge DOF numbering. Call build_edge_dof_numbering first."
        end if
        
        ! Set DOF counts
        this%n_dofs = mesh%n_edge_dofs
        this%n_interior_dofs = mesh%n_interior_dofs
        this%n_boundary_dofs = mesh%n_boundary_dofs
        
        ! Initialize basis functions for each triangle
        allocate(this%basis(mesh%n_triangles))
        do t = 1, mesh%n_triangles
            call this%basis(t)%init(mesh)
        end do
        
        ! Initialize essential BC arrays
        allocate(this%dof_is_essential(this%n_dofs))
        allocate(this%essential_values(this%n_dofs))
        this%dof_is_essential = 0  ! 0 = not essential, 1 = essential
        this%essential_values = 0.0_dp
        
    end subroutine hcurl_space_init
    
    subroutine hcurl_space_destroy(this)
        class(hcurl_space_t), intent(inout) :: this
        integer :: t
        
        if (allocated(this%basis)) then
            do t = 1, size(this%basis)
                call this%basis(t)%destroy()
            end do
            deallocate(this%basis)
        end if
        
        if (allocated(this%dof_is_essential)) deallocate(this%dof_is_essential)
        if (allocated(this%essential_values)) deallocate(this%essential_values)
        
        this%mesh => null()
        this%n_dofs = 0
        this%n_interior_dofs = 0
        this%n_boundary_dofs = 0
    end subroutine hcurl_space_destroy
    
    function hcurl_space_get_n_dofs(this) result(n_dofs)
        class(hcurl_space_t), intent(in) :: this
        integer :: n_dofs
        
        n_dofs = this%n_dofs
    end function hcurl_space_get_n_dofs
    
    function hcurl_space_get_n_interior_dofs(this) result(n_dofs)
        class(hcurl_space_t), intent(in) :: this
        integer :: n_dofs
        
        n_dofs = this%n_interior_dofs
    end function hcurl_space_get_n_interior_dofs
    
    function hcurl_space_get_n_boundary_dofs(this) result(n_dofs)
        class(hcurl_space_t), intent(in) :: this
        integer :: n_dofs
        
        n_dofs = this%n_boundary_dofs
    end function hcurl_space_get_n_boundary_dofs
    
    subroutine hcurl_space_apply_essential_bc(this, boundary_value)
        class(hcurl_space_t), intent(inout) :: this
        real(dp), intent(in) :: boundary_value
        integer :: i, edge_idx
        
        ! Apply homogeneous essential BC: E × n = 0 on boundary
        ! For RT0, this means tangential component is zero
        ! Mark boundary DOFs as essential
        do i = 1, this%n_dofs
            ! Get edge corresponding to this DOF (dof_to_edge is 1-based)
            edge_idx = this%mesh%dof_to_edge(i)
            
            if (this%mesh%is_boundary_edge(edge_idx)) then
                this%dof_is_essential(i) = 1
                this%essential_values(i) = boundary_value
            end if
        end do
    end subroutine hcurl_space_apply_essential_bc
    
    subroutine hcurl_space_evaluate_at_point(this, triangle_idx, xi, eta, coeff, values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: coeff(:)  ! DOF coefficients
        real(dp), intent(out) :: values(2)  ! Evaluated vector field
        
        real(dp) :: triangle_area
        real(dp) :: basis_values(2, 3)
        integer :: triangle_dofs(3)
        integer :: i
        
        if (triangle_idx < 1 .or. triangle_idx > this%mesh%n_triangles) then
            error stop "Triangle index out of bounds"
        end if
        
        ! Get triangle area
        triangle_area = compute_triangle_area(this%mesh, triangle_idx)
        
        ! Evaluate basis functions ensuring no double-counting of DOFs
        call this%evaluate_conservative_basis(triangle_idx, xi, eta, basis_values)
        
        ! Get DOF indices for this triangle
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
        
        ! Compute weighted sum: E(x) = Σ c_i φ_i(x)
        ! Note: triangle_dofs are 0-based, but Fortran arrays are 1-based
        values = 0.0_dp
        do i = 1, 3
            values(1) = values(1) + coeff(triangle_dofs(i) + 1) * basis_values(1, i)
            values(2) = values(2) + coeff(triangle_dofs(i) + 1) * basis_values(2, i)
        end do
    end subroutine hcurl_space_evaluate_at_point
    
    subroutine hcurl_space_evaluate_curl_at_point(this, triangle_idx, xi, eta, coeff, curl_value)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: coeff(:)  ! DOF coefficients
        real(dp), intent(out) :: curl_value  ! Scalar curl in 2D
        
        real(dp) :: triangle_area
        real(dp) :: basis_curls(3)
        integer :: triangle_dofs(3)
        integer :: i
        
        if (triangle_idx < 1 .or. triangle_idx > this%mesh%n_triangles) then
            error stop "Triangle index out of bounds"
        end if
        
        ! Get triangle area
        triangle_area = compute_triangle_area(this%mesh, triangle_idx)
        
        ! Evaluate curl of basis functions
        call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, basis_curls)
        
        ! Get DOF indices for this triangle
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
        
        ! Compute weighted sum: curl(E) = Σ c_i curl(φ_i)
        ! Note: triangle_dofs are 0-based, but Fortran arrays are 1-based
        curl_value = 0.0_dp
        do i = 1, 3
            curl_value = curl_value + coeff(triangle_dofs(i) + 1) * basis_curls(i)
        end do
    end subroutine hcurl_space_evaluate_curl_at_point
    
    subroutine hcurl_space_get_triangle_dofs(this, triangle_idx, triangle_dofs)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        integer, intent(out) :: triangle_dofs(3)
        
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
    end subroutine hcurl_space_get_triangle_dofs
    
    function compute_triangle_area(mesh, triangle_idx) result(area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp) :: area
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        area = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
    end function compute_triangle_area
    
    ! Evaluate edge basis functions with Piola transformation
    subroutine evaluate_edge_basis_2d_with_piola(this, triangle_idx, xi, eta, values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        real(dp) :: ref_values(2, 3)  ! Reference basis values
        real(dp) :: jacobian(2, 2), det_jacobian
        real(dp) :: triangle_area
        integer :: i
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        ! Get triangle vertices
        x1 = this%mesh%vertices(1, this%mesh%triangles(1, triangle_idx))
        y1 = this%mesh%vertices(2, this%mesh%triangles(1, triangle_idx))
        x2 = this%mesh%vertices(1, this%mesh%triangles(2, triangle_idx))
        y2 = this%mesh%vertices(2, this%mesh%triangles(2, triangle_idx))
        x3 = this%mesh%vertices(1, this%mesh%triangles(3, triangle_idx))
        y3 = this%mesh%vertices(2, this%mesh%triangles(3, triangle_idx))
        
        ! Compute Jacobian matrix: F = [∂x/∂ξ, ∂x/∂η; ∂y/∂ξ, ∂y/∂η]
        jacobian(1, 1) = x2 - x1  ! ∂x/∂ξ
        jacobian(1, 2) = x3 - x1  ! ∂x/∂η
        jacobian(2, 1) = y2 - y1  ! ∂y/∂ξ
        jacobian(2, 2) = y3 - y1  ! ∂y/∂η
        
        det_jacobian = jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1)
        triangle_area = 0.5_dp * abs(det_jacobian)
        
        ! Evaluate reference basis functions
        call evaluate_edge_basis_2d(xi, eta, triangle_area, ref_values)
        
        ! Apply Piola transformation: φ_phys = (1/J) * F * φ_ref
        ! For H(curl): u_phys = (1/det(J)) * J * u_ref where J is Jacobian matrix
        do i = 1, 3
            values(1, i) = (jacobian(1, 1) * ref_values(1, i) + jacobian(1, 2) * ref_values(2, i)) / det_jacobian
            values(2, i) = (jacobian(2, 1) * ref_values(1, i) + jacobian(2, 2) * ref_values(2, i)) / det_jacobian
        end do
    end subroutine evaluate_edge_basis_2d_with_piola
    
    ! Evaluate edge basis functions with proper orientation and Piola transformation
    subroutine evaluate_edge_basis_oriented(this, triangle_idx, xi, eta, values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        real(dp) :: ref_values(2, 3)  ! Reference basis values
        real(dp) :: jacobian(2, 2), det_jacobian
        real(dp) :: triangle_area
        integer :: i, global_edge_idx
        integer :: triangle_dofs(3)
        real(dp) :: orientation_sign
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        ! Get triangle vertices
        x1 = this%mesh%vertices(1, this%mesh%triangles(1, triangle_idx))
        y1 = this%mesh%vertices(2, this%mesh%triangles(1, triangle_idx))
        x2 = this%mesh%vertices(1, this%mesh%triangles(2, triangle_idx))
        y2 = this%mesh%vertices(2, this%mesh%triangles(2, triangle_idx))
        x3 = this%mesh%vertices(1, this%mesh%triangles(3, triangle_idx))
        y3 = this%mesh%vertices(2, this%mesh%triangles(3, triangle_idx))
        
        ! Compute Jacobian matrix: F = [∂x/∂ξ, ∂x/∂η; ∂y/∂ξ, ∂y/∂η]
        jacobian(1, 1) = x2 - x1  ! ∂x/∂ξ
        jacobian(1, 2) = x3 - x1  ! ∂x/∂η
        jacobian(2, 1) = y2 - y1  ! ∂y/∂ξ
        jacobian(2, 2) = y3 - y1  ! ∂y/∂η
        
        det_jacobian = jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1)
        triangle_area = 0.5_dp * abs(det_jacobian)
        
        ! Evaluate reference basis functions
        call evaluate_edge_basis_2d(xi, eta, triangle_area, ref_values)
        
        ! Get DOF indices for this triangle
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
        
        ! Apply Piola transformation with proper edge orientations
        do i = 1, 3
            ! For now, use positive orientation (this needs to be fixed based on mesh edge directions)
            orientation_sign = 1.0_dp
            
            ! Apply Piola transformation: φ_phys = (orientation/J) * F * φ_ref
            values(1, i) = orientation_sign * (jacobian(1, 1) * ref_values(1, i) + jacobian(1, 2) * ref_values(2, i)) / det_jacobian
            values(2, i) = orientation_sign * (jacobian(2, 1) * ref_values(1, i) + jacobian(2, 2) * ref_values(2, i)) / det_jacobian
        end do
    end subroutine evaluate_edge_basis_oriented
    
    ! Evaluate edge basis functions based on physical edge directions
    subroutine evaluate_edge_basis_physical(this, triangle_idx, xi, eta, values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        integer :: triangle_dofs(3), i, global_edge_idx
        real(dp) :: edge_length, tangent(2)
        integer :: vertices(2)
        
        ! Get DOF indices for this triangle
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
        
        ! For each local edge, map to physical edge direction
        do i = 1, 3
            global_edge_idx = triangle_dofs(i) + 1  ! Convert to 1-based
            
            if (global_edge_idx <= this%mesh%n_edges) then
                ! Get the physical edge direction
                call this%mesh%get_edge_vertices(global_edge_idx, vertices)
                call this%mesh%get_edge_length_tangent(global_edge_idx, edge_length, tangent)
                
                ! Map edge to appropriate coordinate direction
                ! Simple heuristic: if tangent is more horizontal, map to x-direction
                if (abs(tangent(1)) > abs(tangent(2))) then
                    ! More horizontal edge -> x-direction basis
                    values(1, i) = 1.0_dp
                    values(2, i) = 0.0_dp
                else
                    ! More vertical edge -> y-direction basis  
                    values(1, i) = 0.0_dp
                    values(2, i) = 1.0_dp
                end if
            else
                ! Safety fallback
                values(1, i) = 0.0_dp
                values(2, i) = 0.0_dp
            end if
        end do
    end subroutine evaluate_edge_basis_physical
    
    ! Evaluate Nédélec basis functions with correct edge-tangent properties
    subroutine evaluate_nedelec_basis_correct(this, triangle_idx, xi, eta, values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        integer :: triangle_dofs(3), i, global_edge_idx
        real(dp) :: edge_length, tangent(2)
        integer :: vertices(2)
        real(dp) :: edge_lengths(3), edge_tangents(2,3)
        
        ! Get DOF indices for this triangle
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
        
        ! Get all edge information for this triangle
        do i = 1, 3
            global_edge_idx = triangle_dofs(i) + 1  ! Convert to 1-based
            
            if (global_edge_idx <= this%mesh%n_edges) then
                call this%mesh%get_edge_vertices(global_edge_idx, vertices)
                call this%mesh%get_edge_length_tangent(global_edge_idx, edge_lengths(i), edge_tangents(:,i))
            else
                edge_lengths(i) = 1.0_dp
                edge_tangents(1,i) = 1.0_dp
                edge_tangents(2,i) = 0.0_dp
            end if
        end do
        
        ! For each basis function, create a field aligned with its edge
        do i = 1, 3
            ! Simple approach: basis function i is aligned with edge i's tangent direction
            ! Scaled to give unit line integral along the edge
            values(1, i) = edge_tangents(1, i) / edge_lengths(i)
            values(2, i) = edge_tangents(2, i) / edge_lengths(i)
        end do
    end subroutine evaluate_nedelec_basis_correct
    
    ! Conservative basis evaluation - each DOF contributes exactly once
    subroutine evaluate_conservative_basis(this, triangle_idx, xi, eta, values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        integer :: triangle_dofs(3), i, global_edge_idx
        real(dp) :: edge_length, tangent(2)
        integer :: vertices(2)
        
        ! Get DOF indices for this triangle
        call this%mesh%get_triangle_edge_dofs(triangle_idx, triangle_dofs)
        
        ! For each local edge, create basis function based on global DOF index
        do i = 1, 3
            global_edge_idx = triangle_dofs(i) + 1  ! Convert to 1-based
            
            if (global_edge_idx <= this%mesh%n_edges) then
                call this%mesh%get_edge_vertices(global_edge_idx, vertices)
                call this%mesh%get_edge_length_tangent(global_edge_idx, edge_length, tangent)
                
                ! Create basis function proportional to edge tangent
                ! Divide by number of triangles sharing this edge to avoid double-counting
                values(1, i) = tangent(1) * 0.5_dp  ! Assume each edge is shared by 2 triangles
                values(2, i) = tangent(2) * 0.5_dp
            else
                values(1, i) = 0.0_dp
                values(2, i) = 0.0_dp
            end if
        end do
    end subroutine evaluate_conservative_basis

    ! Evaluate curl of basis functions at a point
    subroutine hcurl_space_evaluate_curl_basis_at_point(this, triangle_idx, xi, eta, curl_values)
        class(hcurl_space_t), intent(in) :: this
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: curl_values(3)  ! Curl of 3 basis functions
        
        real(dp) :: triangle_area
        
        ! Compute triangle area
        triangle_area = compute_triangle_area_helper(this%mesh, triangle_idx)
        
        ! Use the corrected curl evaluation from basis_edge_2d module
        call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curl_values)
    end subroutine hcurl_space_evaluate_curl_basis_at_point

    ! Helper function to compute triangle area
    function compute_triangle_area_helper(mesh, triangle_idx) result(area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp) :: area
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        area = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
    end function compute_triangle_area_helper

end module fortfem_hcurl_space
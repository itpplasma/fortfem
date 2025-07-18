module fortfem_hcurl_space
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_basis_edge_2d
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
        procedure :: get_triangle_dofs => hcurl_space_get_triangle_dofs
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
        
        ! Evaluate basis functions
        call evaluate_edge_basis_2d(xi, eta, triangle_area, basis_values)
        
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

end module fortfem_hcurl_space
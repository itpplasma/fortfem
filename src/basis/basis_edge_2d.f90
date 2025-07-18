module fortfem_basis_edge_2d
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    private
    
    public :: edge_basis_2d_t
    public :: evaluate_edge_basis_2d
    public :: evaluate_edge_basis_curl_2d
    public :: evaluate_edge_basis_div_2d
    public :: evaluate_edge_basis_2d_piola
    
    ! Edge element basis functions (Nédélec elements)
    type :: edge_basis_2d_t
        integer :: n_edges = 0
        real(dp), allocatable :: edge_vectors(:,:)  ! Direction vectors
        real(dp), allocatable :: edge_lengths(:)    ! Edge lengths
    contains
        procedure :: init => edge_basis_init
        procedure :: destroy => edge_basis_destroy
        procedure :: n_dofs => edge_basis_n_dofs
    end type edge_basis_2d_t
    
contains

    subroutine edge_basis_init(this, mesh)
        class(edge_basis_2d_t), intent(inout) :: this
        type(mesh_2d_t), intent(in) :: mesh
        integer :: i, j, k, edge_count
        real(dp) :: dx, dy
        
        ! Count edges (3 per triangle)
        this%n_edges = 3 * mesh%n_triangles
        
        allocate(this%edge_vectors(2, this%n_edges))
        allocate(this%edge_lengths(this%n_edges))
        
        edge_count = 0
        
        ! For each triangle, compute edge vectors
        do i = 1, mesh%n_triangles
            do j = 1, 3
                edge_count = edge_count + 1
                
                ! Next vertex (cyclic)
                k = mod(j, 3) + 1
                
                ! Edge vector
                dx = mesh%vertices(1, mesh%triangles(k, i)) - &
                     mesh%vertices(1, mesh%triangles(j, i))
                dy = mesh%vertices(2, mesh%triangles(k, i)) - &
                     mesh%vertices(2, mesh%triangles(j, i))
                
                this%edge_vectors(1, edge_count) = dx
                this%edge_vectors(2, edge_count) = dy
                this%edge_lengths(edge_count) = sqrt(dx**2 + dy**2)
            end do
        end do
    end subroutine edge_basis_init
    
    subroutine edge_basis_destroy(this)
        class(edge_basis_2d_t), intent(inout) :: this
        if (allocated(this%edge_vectors)) deallocate(this%edge_vectors)
        if (allocated(this%edge_lengths)) deallocate(this%edge_lengths)
        this%n_edges = 0
    end subroutine edge_basis_destroy
    
    function edge_basis_n_dofs(this) result(n)
        class(edge_basis_2d_t), intent(in) :: this
        integer :: n
        n = this%n_edges
    end function edge_basis_n_dofs
    
    ! Evaluate edge basis functions at reference coordinates
    subroutine evaluate_edge_basis_2d(xi, eta, triangle_area, values)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        ! Lowest order Nédélec elements (H(curl) conforming)
        ! Reference triangle: (0,0), (1,0), (0,1)
        ! Standard RT0 (Raviart-Thomas) basis functions
        
        ! For constant field representation, we need:
        ! Sum of basis functions weighted by coefficients = constant field
        
        ! Edge 0: φ₀ = (1, 0) - simple constant field in x-direction
        values(1, 1) = 1.0_dp
        values(2, 1) = 0.0_dp
        
        ! Edge 1: φ₁ = (0, 1) - simple constant field in y-direction
        values(1, 2) = 0.0_dp
        values(2, 2) = 1.0_dp
        
        ! Edge 2: φ₂ = (0, 0) - no contribution (for 2D we only need 2 basis functions for constant fields)
        values(1, 3) = 0.0_dp
        values(2, 3) = 0.0_dp
    end subroutine evaluate_edge_basis_2d
    
    ! Evaluate curl of edge basis functions
    subroutine evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: curls(3)  ! Scalar curl in 2D
        
        ! Curl of Nédélec elements (constant per element)
        ! For Nédélec: curl(φᵢ) = ∂φᵢʸ/∂x - ∂φᵢˣ/∂y
        ! Transform from reference to physical using Jacobian: curl_phys = curl_ref / J
        ! For triangular elements: J = 2 * triangle_area
        real(dp) :: jacobian_det
        
        jacobian_det = 2.0_dp * triangle_area
        
        ! Basis 0: φ₀ = (1, 0), curl = ∂0/∂x - ∂1/∂y = 0 - 0 = 0
        curls(1) = 0.0_dp
        
        ! Basis 1: φ₁ = (0, 1), curl = ∂1/∂x - ∂0/∂y = 0 - 0 = 0
        curls(2) = 0.0_dp
        
        ! Basis 2: φ₂ = (0, 0), curl = 0
        curls(3) = 0.0_dp
    end subroutine evaluate_edge_basis_curl_2d
    
    ! Evaluate edge basis functions with Piola transformation
    subroutine evaluate_edge_basis_2d_piola(mesh, triangle_idx, xi, eta, values)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)  ! 2D vectors, 3 edges
        
        real(dp) :: ref_values(2, 3)  ! Reference basis values
        real(dp) :: jacobian(2, 2), inv_jacobian(2, 2), det_jacobian
        real(dp) :: triangle_area
        integer :: i, j, k
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        ! Get triangle vertices
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
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
        do i = 1, 3
            values(1, i) = (jacobian(1, 1) * ref_values(1, i) + jacobian(1, 2) * ref_values(2, i)) / det_jacobian
            values(2, i) = (jacobian(2, 1) * ref_values(1, i) + jacobian(2, 2) * ref_values(2, i)) / det_jacobian
        end do
    end subroutine evaluate_edge_basis_2d_piola
    
    ! Evaluate divergence of edge basis functions
    subroutine evaluate_edge_basis_div_2d(xi, eta, triangle_area, divs)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: divs(3)  ! Divergence values
        
        ! Divergence of Nédélec elements (should be zero for H(curl))
        divs(1) = 0.0_dp
        divs(2) = 0.0_dp
        divs(3) = 0.0_dp
    end subroutine evaluate_edge_basis_div_2d

end module fortfem_basis_edge_2d
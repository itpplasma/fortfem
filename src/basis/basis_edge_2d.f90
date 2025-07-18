module fortfem_basis_edge_2d
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    private
    
    public :: edge_basis_2d_t
    public :: evaluate_edge_basis_2d
    public :: evaluate_edge_basis_curl_2d
    public :: evaluate_edge_basis_div_2d
    
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
        
        ! Lowest order Nédélec elements (RT0)
        ! Edge 1: (1-eta, 0)
        values(1, 1) = 1.0_dp - eta
        values(2, 1) = 0.0_dp
        
        ! Edge 2: (0, xi)
        values(1, 2) = 0.0_dp
        values(2, 2) = xi
        
        ! Edge 3: (eta, 1-xi)
        values(1, 3) = eta
        values(2, 3) = 1.0_dp - xi
    end subroutine evaluate_edge_basis_2d
    
    ! Evaluate curl of edge basis functions
    subroutine evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: curls(3)  ! Scalar curl in 2D
        
        ! Curl of Nédélec elements (constant per element)
        curls(1) = 1.0_dp / triangle_area
        curls(2) = 1.0_dp / triangle_area
        curls(3) = -1.0_dp / triangle_area
    end subroutine evaluate_edge_basis_curl_2d
    
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
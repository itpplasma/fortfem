module fortfem_mesh_1d
    use fortfem_kinds
    implicit none
    private
    
    type, public :: mesh_1d_t
        integer :: n_nodes = 0
        integer :: n_elements = 0
        real(dp), allocatable :: nodes(:)
        integer, allocatable :: connectivity(:,:)  ! (2, n_elements)
    contains
        procedure :: create_uniform
        procedure :: element_size
        procedure :: deallocate => mesh_deallocate
    end type mesh_1d_t
    
contains

    subroutine create_uniform(this, n_nodes, x_min, x_max)
        class(mesh_1d_t), intent(out) :: this
        integer, intent(in) :: n_nodes
        real(dp), intent(in) :: x_min, x_max
        
        integer :: i
        real(dp) :: dx
        
        this%n_nodes = n_nodes
        this%n_elements = n_nodes - 1
        
        ! Allocate arrays
        allocate(this%nodes(n_nodes))
        allocate(this%connectivity(2, this%n_elements))
        
        ! Create uniform node positions
        dx = (x_max - x_min) / real(n_nodes - 1, dp)
        do i = 1, n_nodes
            this%nodes(i) = x_min + (i - 1) * dx
        end do
        
        ! Create connectivity
        do i = 1, this%n_elements
            this%connectivity(1, i) = i
            this%connectivity(2, i) = i + 1
        end do
        
    end subroutine create_uniform
    
    function element_size(this, elem_id) result(h)
        class(mesh_1d_t), intent(in) :: this
        integer, intent(in) :: elem_id
        real(dp) :: h
        
        integer :: n1, n2
        
        n1 = this%connectivity(1, elem_id)
        n2 = this%connectivity(2, elem_id)
        h = abs(this%nodes(n2) - this%nodes(n1))
        
    end function element_size
    
    subroutine mesh_deallocate(this)
        class(mesh_1d_t), intent(inout) :: this
        
        if (allocated(this%nodes)) deallocate(this%nodes)
        if (allocated(this%connectivity)) deallocate(this%connectivity)
        this%n_nodes = 0
        this%n_elements = 0
        
    end subroutine mesh_deallocate

end module fortfem_mesh_1d
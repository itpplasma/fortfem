module function_space_module
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    private
    
    public :: function_space_t, P1, P2, create_P1_space, create_P2_space
    public :: trial_function_t, test_function_t
    public :: vector_function_space_t, vector_trial_function_t, vector_test_function_t
    public :: create_edge_space  ! For Nédélec edge elements
    
    ! Element types
    integer, parameter :: ELEMENT_P1 = 1
    integer, parameter :: ELEMENT_P2 = 2
    integer, parameter :: ELEMENT_RT0 = 3
    integer, parameter :: ELEMENT_EDGE = 4  ! Nédélec edge elements
    
    ! Function space definition
    type :: function_space_t
        type(mesh_2d_t) :: mesh
        integer :: element_type = ELEMENT_P1
        integer :: n_dofs = 0
        integer :: n_components = 1  ! 1 for scalar, 2+ for vector
    contains
        procedure :: init => function_space_init
        procedure :: destroy => function_space_destroy
        procedure :: n_dofs_per_element
        procedure :: get_dof_indices
    end type function_space_t
    
    ! Scalar trial function (unknown)
    type :: trial_function_t
        type(function_space_t), allocatable :: space
        character(len=32) :: name = ""
        real(dp), allocatable :: coefficients(:)
    contains
        procedure :: init => trial_function_init
        procedure :: destroy => trial_function_destroy
        procedure :: assign_values
        generic :: assignment(=) => trial_function_assign
        procedure :: trial_function_assign
    end type trial_function_t
    
    ! Scalar test function
    type :: test_function_t
        type(function_space_t), allocatable :: space
        character(len=32) :: name = ""
    contains
        procedure :: init => test_function_init
        procedure :: destroy => test_function_destroy
        generic :: assignment(=) => test_function_assign
        procedure :: test_function_assign
    end type test_function_t
    
    ! Vector function space
    type :: vector_function_space_t
        type(function_space_t) :: scalar_space
        integer :: dim = 2  ! spatial dimension
    contains
        procedure :: init => vector_function_space_init
        procedure :: destroy => vector_function_space_destroy
        procedure :: total_dofs
    end type vector_function_space_t
    
    ! Vector trial function
    type :: vector_trial_function_t
        type(vector_function_space_t), allocatable :: space
        character(len=32) :: name = ""
        real(dp), allocatable :: coefficients(:)
    contains
        procedure :: init => vector_trial_function_init
        procedure :: destroy => vector_trial_function_destroy
        procedure :: component
        generic :: assignment(=) => vector_trial_function_assign
        procedure :: vector_trial_function_assign
    end type vector_trial_function_t
    
    ! Vector test function
    type :: vector_test_function_t
        type(vector_function_space_t), allocatable :: space
        character(len=32) :: name = ""
    contains
        procedure :: init => vector_test_function_init
        procedure :: destroy => vector_test_function_destroy
        procedure :: component => vector_test_component
        generic :: assignment(=) => vector_test_function_assign
        procedure :: vector_test_function_assign
    end type vector_test_function_t
    
contains

    ! Function space constructor functions
    function P1(mesh) result(space)
        type(mesh_2d_t), intent(in) :: mesh
        type(function_space_t) :: space
        
        call space%init(mesh, ELEMENT_P1)
    end function P1
    
    function P2(mesh) result(space)
        type(mesh_2d_t), intent(in) :: mesh
        type(function_space_t) :: space
        
        call space%init(mesh, ELEMENT_P2)
    end function P2

    ! Subroutine versions to avoid returning allocatables
    subroutine create_P1_space(mesh, space)
        type(mesh_2d_t), intent(in) :: mesh
        type(function_space_t), intent(out) :: space
        
        call space%init(mesh, ELEMENT_P1)
    end subroutine create_P1_space

    subroutine create_P2_space(mesh, space)
        type(mesh_2d_t), intent(in) :: mesh
        type(function_space_t), intent(out) :: space
        
        call space%init(mesh, ELEMENT_P2)
    end subroutine create_P2_space
    
    ! Function space methods
    subroutine function_space_init(this, mesh, element_type)
        class(function_space_t), intent(inout) :: this
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: element_type
        
        this%mesh = mesh
        this%element_type = element_type
        this%n_components = 1
        
        select case (element_type)
        case (ELEMENT_P1)
            this%n_dofs = mesh%n_vertices
        case (ELEMENT_P2)
            ! P2 has vertices + edge midpoints
            this%n_dofs = mesh%n_vertices + mesh%n_edges
        case default
            error stop "Unknown element type"
        end select
    end subroutine function_space_init
    
    subroutine function_space_destroy(this)
        class(function_space_t), intent(inout) :: this
        
        call this%mesh%destroy()
        this%n_dofs = 0
        this%element_type = 0
    end subroutine function_space_destroy
    
    function n_dofs_per_element(this) result(n)
        class(function_space_t), intent(in) :: this
        integer :: n
        
        select case (this%element_type)
        case (ELEMENT_P1)
            n = 3  ! 3 vertices per triangle
        case (ELEMENT_P2)
            n = 6  ! 3 vertices + 3 edges per triangle
        case default
            n = 0
        end select
    end function n_dofs_per_element
    
    subroutine get_dof_indices(this, element, dof_indices)
        class(function_space_t), intent(in) :: this
        integer, intent(in) :: element
        integer, intent(out) :: dof_indices(:)
        
        select case (this%element_type)
        case (ELEMENT_P1)
            dof_indices(1:3) = this%mesh%triangles(:, element)
        case (ELEMENT_P2)
            ! P2 elements have 6 DOFs per triangle: 3 vertices + 3 edge midpoints
            dof_indices(1:3) = this%mesh%triangles(:, element)  ! Vertices
            ! Edge midpoints (simplified mapping)
            dof_indices(4) = this%mesh%n_vertices + 3*(element-1) + 1  ! Edge 1-2
            dof_indices(5) = this%mesh%n_vertices + 3*(element-1) + 2  ! Edge 2-3
            dof_indices(6) = this%mesh%n_vertices + 3*(element-1) + 3  ! Edge 3-1
        end select
    end subroutine get_dof_indices
    
    ! Trial function methods
    subroutine trial_function_init(this, space, name)
        class(trial_function_t), intent(inout) :: this
        type(function_space_t), intent(in) :: space
        character(len=*), intent(in), optional :: name
        
        allocate(this%space)
        this%space = space
        if (present(name)) this%name = name
        
        allocate(this%coefficients(space%n_dofs))
        this%coefficients = 0.0_dp
    end subroutine trial_function_init
    
    subroutine trial_function_destroy(this)
        class(trial_function_t), intent(inout) :: this
        
        if (allocated(this%coefficients)) deallocate(this%coefficients)
        if (allocated(this%space)) deallocate(this%space)
    end subroutine trial_function_destroy
    
    subroutine assign_values(this, values)
        class(trial_function_t), intent(inout) :: this
        real(dp), intent(in) :: values(:)
        
        if (size(values) /= size(this%coefficients)) then
            error stop "Size mismatch in trial function assignment"
        end if
        
        this%coefficients = values
    end subroutine assign_values
    
    ! Test function methods
    subroutine test_function_init(this, space, name)
        class(test_function_t), intent(inout) :: this
        type(function_space_t), intent(in) :: space
        character(len=*), intent(in), optional :: name
        
        allocate(this%space)
        this%space = space
        if (present(name)) this%name = name
    end subroutine test_function_init
    
    subroutine test_function_destroy(this)
        class(test_function_t), intent(inout) :: this
        
        if (allocated(this%space)) deallocate(this%space)
    end subroutine test_function_destroy
    
    ! Vector function space methods
    subroutine vector_function_space_init(this, scalar_space, dim)
        class(vector_function_space_t), intent(inout) :: this
        type(function_space_t), intent(in) :: scalar_space
        integer, intent(in) :: dim
        
        this%scalar_space = scalar_space
        this%dim = dim
    end subroutine vector_function_space_init
    
    subroutine vector_function_space_destroy(this)
        class(vector_function_space_t), intent(inout) :: this
        
        call this%scalar_space%destroy()
        this%dim = 0
    end subroutine vector_function_space_destroy
    
    function total_dofs(this) result(n)
        class(vector_function_space_t), intent(in) :: this
        integer :: n
        
        n = this%dim * this%scalar_space%n_dofs
    end function total_dofs
    
    ! Vector trial function methods
    subroutine vector_trial_function_init(this, space, name)
        class(vector_trial_function_t), intent(inout) :: this
        type(vector_function_space_t), intent(in) :: space
        character(len=*), intent(in), optional :: name
        
        allocate(this%space)
        this%space = space
        if (present(name)) this%name = name
        
        allocate(this%coefficients(space%total_dofs()))
        this%coefficients = 0.0_dp
    end subroutine vector_trial_function_init
    
    subroutine vector_trial_function_destroy(this)
        class(vector_trial_function_t), intent(inout) :: this
        
        if (allocated(this%coefficients)) deallocate(this%coefficients)
        if (allocated(this%space)) deallocate(this%space)
    end subroutine vector_trial_function_destroy
    
    subroutine component(this, comp, u_comp)
        class(vector_trial_function_t), intent(in) :: this
        integer, intent(in) :: comp
        type(trial_function_t), intent(out) :: u_comp
        
        ! This would extract the component - simplified for now
        call u_comp%init(this%space%scalar_space, &
                         trim(this%name) // "_" // char(comp + 48))
    end subroutine component
    
    ! Vector test function methods
    subroutine vector_test_function_init(this, space, name)
        class(vector_test_function_t), intent(inout) :: this
        type(vector_function_space_t), intent(in) :: space
        character(len=*), intent(in), optional :: name
        
        allocate(this%space)
        this%space = space
        if (present(name)) this%name = name
    end subroutine vector_test_function_init
    
    subroutine vector_test_function_destroy(this)
        class(vector_test_function_t), intent(inout) :: this
        
        if (allocated(this%space)) deallocate(this%space)
    end subroutine vector_test_function_destroy
    
    subroutine vector_test_component(this, comp, v_comp)
        class(vector_test_function_t), intent(in) :: this
        integer, intent(in) :: comp
        type(test_function_t), intent(out) :: v_comp
        
        call v_comp%init(this%space%scalar_space, &
                         trim(this%name) // "_" // char(comp + 48))
    end subroutine vector_test_component

    ! Deep-copy assignment operators
    subroutine trial_function_assign(this, other)
        class(trial_function_t), intent(inout) :: this
        type(trial_function_t), intent(in) :: other
        
        call this%destroy()
        
        if (allocated(other%space)) then
            allocate(this%space)
            this%space = other%space
        end if
        
        this%name = other%name
        
        if (allocated(other%coefficients)) then
            allocate(this%coefficients(size(other%coefficients)))
            this%coefficients = other%coefficients
        end if
    end subroutine trial_function_assign
    
    subroutine test_function_assign(this, other)
        class(test_function_t), intent(inout) :: this
        type(test_function_t), intent(in) :: other
        
        call this%destroy()
        
        if (allocated(other%space)) then
            allocate(this%space)
            this%space = other%space
        end if
        
        this%name = other%name
    end subroutine test_function_assign
    
    subroutine vector_trial_function_assign(this, other)
        class(vector_trial_function_t), intent(inout) :: this
        type(vector_trial_function_t), intent(in) :: other
        
        call this%destroy()
        
        if (allocated(other%space)) then
            allocate(this%space)
            this%space = other%space
        end if
        
        this%name = other%name
        
        if (allocated(other%coefficients)) then
            allocate(this%coefficients(size(other%coefficients)))
            this%coefficients = other%coefficients
        end if
    end subroutine vector_trial_function_assign
    
    subroutine vector_test_function_assign(this, other)
        class(vector_test_function_t), intent(inout) :: this
        type(vector_test_function_t), intent(in) :: other
        
        call this%destroy()
        
        if (allocated(other%space)) then
            allocate(this%space)
            this%space = other%space
        end if
        
        this%name = other%name
    end subroutine vector_test_function_assign
    
    ! Create edge element function space (Nédélec H(curl))
    subroutine create_edge_space(mesh, space)
        type(mesh_2d_t), intent(in) :: mesh
        type(function_space_t), intent(out) :: space
        
        space%mesh = mesh
        space%element_type = ELEMENT_EDGE
        space%n_components = 2  ! Vector field in 2D
        
        ! For lowest order Nédélec: 1 DOF per edge
        ! Each triangle has 3 edges, but edges are shared
        ! Approximate: 3 * n_triangles (will be refined in actual implementation)
        space%n_dofs = 3 * mesh%n_triangles
    end subroutine create_edge_space

end module function_space_module
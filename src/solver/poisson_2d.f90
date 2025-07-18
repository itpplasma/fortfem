module poisson_2d_module
    use fortfem_kinds
    use fortfem_mesh_2d
    use assembly_2d_module
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    implicit none
    private
    
    public :: poisson_2d_t
    
    type :: poisson_2d_t
        type(mesh_2d_t) :: mesh
        type(assembly_2d_t) :: assembly
        class(linear_solver_t), allocatable :: solver
        real(dp), allocatable :: solution(:)
        real(dp), allocatable :: rhs(:)
        integer, allocatable :: dirichlet_nodes(:)
        real(dp), allocatable :: dirichlet_values(:)
        integer :: n_dirichlet = 0
    contains
        procedure :: init
        procedure :: destroy
        procedure :: set_mesh
        procedure :: set_dirichlet_bc
        procedure :: solve
        procedure :: get_solution
        procedure :: apply_boundary_conditions
    end type poisson_2d_t
    
contains

    subroutine init(this, solver_type)
        class(poisson_2d_t), intent(inout) :: this
        character(len=*), intent(in), optional :: solver_type
        character(len=20) :: solver_name
        
        solver_name = "lapack"
        if (present(solver_type)) solver_name = solver_type
        
        ! Allocate solver
        select case (trim(solver_name))
        case ("lapack")
            allocate(lapack_dense_solver_t :: this%solver)
        case ("suitesparse", "umfpack")
            allocate(suitesparse_solver_t :: this%solver)
        case default
            error stop "Unknown solver type"
        end select
        
    end subroutine init
    
    subroutine destroy(this)
        class(poisson_2d_t), intent(inout) :: this
        
        call this%mesh%destroy()
        call this%assembly%destroy()
        if (allocated(this%solver)) then
            call this%solver%destroy()
            deallocate(this%solver)
        end if
        if (allocated(this%solution)) deallocate(this%solution)
        if (allocated(this%rhs)) deallocate(this%rhs)
        if (allocated(this%dirichlet_nodes)) deallocate(this%dirichlet_nodes)
        if (allocated(this%dirichlet_values)) deallocate(this%dirichlet_values)
        
    end subroutine destroy
    
    subroutine set_mesh(this, mesh)
        class(poisson_2d_t), intent(inout) :: this
        type(mesh_2d_t), intent(in) :: mesh
        
        this%mesh = mesh
        
        ! Allocate solution and RHS
        if (allocated(this%solution)) deallocate(this%solution)
        if (allocated(this%rhs)) deallocate(this%rhs)
        allocate(this%solution(mesh%n_vertices))
        allocate(this%rhs(mesh%n_vertices))
        
        ! Initialize assembly with estimated non-zeros
        call this%assembly%init(mesh%n_vertices, 9 * mesh%n_vertices)
        
    end subroutine set_mesh
    
    subroutine set_dirichlet_bc(this, nodes, values)
        class(poisson_2d_t), intent(inout) :: this
        integer, intent(in) :: nodes(:)
        real(dp), intent(in) :: values(:)
        
        this%n_dirichlet = size(nodes)
        
        if (allocated(this%dirichlet_nodes)) deallocate(this%dirichlet_nodes)
        if (allocated(this%dirichlet_values)) deallocate(this%dirichlet_values)
        
        allocate(this%dirichlet_nodes(this%n_dirichlet))
        allocate(this%dirichlet_values(this%n_dirichlet))
        
        this%dirichlet_nodes = nodes
        this%dirichlet_values = values
        
    end subroutine set_dirichlet_bc
    
    subroutine solve(this, source_func)
        class(poisson_2d_t), intent(inout) :: this
        interface
            pure function source_func(x, y) result(val)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: val
            end function source_func
        end interface
        
        type(triplet_matrix_t) :: triplet
        type(csr_matrix_t) :: csr
        real(dp), allocatable :: x(:)
        integer :: i, j, info
        
        ! Assemble system
        call this%assembly%assemble_global(this%mesh, source_func, triplet, this%rhs)
        
        ! Apply boundary conditions
        call this%apply_boundary_conditions(triplet)
        
        ! Convert to CSR format
        call triplet%to_csr(csr)
        
        ! Allocate solution vector
        allocate(x(this%mesh%n_vertices))
        
        ! Initialize solver and solve
        call this%solver%init()
        call this%solver%solve(csr, this%rhs, x, info)
        
        if (info /= 0) then
            error stop "Linear solver failed"
        end if
        
        ! Copy solution
        this%solution = x
        
        ! Clean up
        deallocate(x)
        call csr%destroy()
        
        call triplet%destroy()
        
    end subroutine solve
    
    subroutine apply_boundary_conditions(this, matrix)
        class(poisson_2d_t), intent(inout) :: this
        type(triplet_matrix_t), intent(inout) :: matrix
        
        integer :: i, j, k, node
        real(dp) :: diag_val, bc_val
        logical :: is_bc_node
        
        ! For each Dirichlet node, modify matrix and RHS
        do i = 1, this%n_dirichlet
            node = this%dirichlet_nodes(i)
            bc_val = this%dirichlet_values(i)
            
            ! Find diagonal entry
            diag_val = 0.0_dp
            do j = 1, matrix%nnz
                if (matrix%rows(j) == node .and. matrix%cols(j) == node) then
                    diag_val = abs(matrix%values(j))
                    if (diag_val < 1e-10) diag_val = 1.0_dp
                    exit
                end if
            end do
            
            ! Zero out row and column, set diagonal
            do j = 1, matrix%nnz
                if (matrix%rows(j) == node) then
                    if (matrix%cols(j) == node) then
                        matrix%values(j) = diag_val
                    else
                        matrix%values(j) = 0.0_dp
                    end if
                else if (matrix%cols(j) == node) then
                    ! Modify RHS for non-BC nodes
                    is_bc_node = .false.
                    do k = 1, this%n_dirichlet
                        if (matrix%rows(j) == this%dirichlet_nodes(k)) then
                            is_bc_node = .true.
                            exit
                        end if
                    end do
                    
                    if (.not. is_bc_node) then
                        this%rhs(matrix%rows(j)) = this%rhs(matrix%rows(j)) - &
                                                  matrix%values(j) * bc_val
                    end if
                    matrix%values(j) = 0.0_dp
                end if
            end do
            
            ! Set RHS
            this%rhs(node) = diag_val * bc_val
        end do
        
    end subroutine apply_boundary_conditions
    
    function get_solution(this) result(sol)
        class(poisson_2d_t), intent(in) :: this
        real(dp), allocatable :: sol(:)
        
        if (allocated(this%solution)) then
            allocate(sol(size(this%solution)))
            sol = this%solution
        else
            allocate(sol(0))
        end if
        
    end function get_solution

end module poisson_2d_module
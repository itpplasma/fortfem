module fortfem_poisson_1d_sparse
    use fortfem_kinds
    use fortfem_mesh_1d
    use fortfem_assembly_1d
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    implicit none
    private
    
    ! Re-export the interface
    public :: source_function
    
    type, public :: poisson_1d_sparse_solver_t
        type(mesh_1d_t) :: mesh
        type(triplet_matrix_t) :: K_triplet
        type(csr_matrix_t) :: K_csr
        real(dp), allocatable :: f(:)
        class(linear_solver_t), allocatable :: solver
        character(len=32) :: solver_type = "lapack_dense"
    contains
        procedure :: init
        procedure :: set_solver
        procedure :: solve
        procedure :: deallocate => solver_deallocate
    end type poisson_1d_sparse_solver_t
    
contains

    subroutine init(this, n_nodes, x_min, x_max)
        class(poisson_1d_sparse_solver_t), intent(out) :: this
        integer, intent(in) :: n_nodes
        real(dp), intent(in) :: x_min, x_max
        
        ! Create mesh
        call this%mesh%create_uniform(n_nodes, x_min, x_max)
        
        ! Allocate load vector
        allocate(this%f(n_nodes))
        this%f = 0.0_dp
        
        ! Initialize triplet matrix with estimated non-zeros
        ! For 1D FEM with P1 elements: max 3 non-zeros per row
        call this%K_triplet%init(n=n_nodes, max_nnz=3*n_nodes)
        
        ! Create default solver
        call create_solver(this%solver, this%solver_type)
        call this%solver%init()
        
    end subroutine init
    
    subroutine set_solver(this, solver_type)
        class(poisson_1d_sparse_solver_t), intent(inout) :: this
        character(len=*), intent(in) :: solver_type
        
        ! Destroy old solver
        if (allocated(this%solver)) then
            call this%solver%destroy()
            deallocate(this%solver)
        end if
        
        ! Create new solver
        this%solver_type = solver_type
        call create_solver(this%solver, solver_type)
        call this%solver%init()
        
    end subroutine set_solver
    
    subroutine solve(this, source, u)
        class(poisson_1d_sparse_solver_t), intent(inout) :: this
        procedure(source_function) :: source
        real(dp), allocatable, intent(out) :: u(:)
        
        real(dp) :: K_elem(2,2), f_elem(2)
        real(dp) :: h, x_elem
        integer :: e, i, n1, n2
        integer :: n, info
        
        n = this%mesh%n_nodes
        
        ! Reset matrices
        this%K_triplet%nnz = 0
        this%f = 0.0_dp
        
        ! Assemble in triplet format
        do e = 1, this%mesh%n_elements
            n1 = this%mesh%connectivity(1, e)
            n2 = this%mesh%connectivity(2, e)
            h = this%mesh%element_size(e)
            x_elem = this%mesh%nodes(n1)
            
            ! Get element matrices
            call element_stiffness_matrix(h, K_elem)
            call element_load_vector_physical(h, x_elem, source, f_elem)
            
            ! Add to triplet matrix (accumulates automatically)
            call this%K_triplet%add(n1, n1, K_elem(1, 1))
            call this%K_triplet%add(n1, n2, K_elem(1, 2))
            call this%K_triplet%add(n2, n1, K_elem(2, 1))
            call this%K_triplet%add(n2, n2, K_elem(2, 2))
            
            ! Assemble load vector
            this%f(n1) = this%f(n1) + f_elem(1)
            this%f(n2) = this%f(n2) + f_elem(2)
        end do
        
        ! Apply Dirichlet BC: u(0) = u(1) = 0
        ! Method: zero out rows and set diagonal to 1
        
        ! First node
        call zero_row(this%K_triplet, 1)
        call this%K_triplet%add(1, 1, 1.0_dp)
        this%f(1) = 0.0_dp
        
        ! Last node
        call zero_row(this%K_triplet, n)
        call this%K_triplet%add(n, n, 1.0_dp)
        this%f(n) = 0.0_dp
        
        ! Convert to CSR
        call this%K_triplet%to_csr(this%K_csr)
        
        ! Solve
        allocate(u(n))
        call this%solver%solve(this%K_csr, this%f, u, info)
        
        if (info /= 0) then
            error stop "Solver failed"
        end if
        
    end subroutine solve
    
    subroutine element_load_vector_physical(h, x_start, source, f_elem)
        real(dp), intent(in) :: h, x_start
        procedure(source_function) :: source
        real(dp), intent(out) :: f_elem(2)
        
        ! 2-point Gauss quadrature
        real(dp), parameter :: xi1 = 0.5_dp - 1.0_dp/(2.0_dp*sqrt(3.0_dp))
        real(dp), parameter :: xi2 = 0.5_dp + 1.0_dp/(2.0_dp*sqrt(3.0_dp))
        real(dp), parameter :: w = 0.5_dp
        
        real(dp) :: x1, x2
        
        ! Map to physical coordinates
        x1 = x_start + xi1 * h
        x2 = x_start + xi2 * h
        
        ! Integrate
        f_elem(1) = h * w * (source(x1) * (1.0_dp - xi1) + &
                            source(x2) * (1.0_dp - xi2))
        f_elem(2) = h * w * (source(x1) * xi1 + source(x2) * xi2)
        
    end subroutine element_load_vector_physical
    
    subroutine zero_row(triplet, row)
        type(triplet_matrix_t), intent(inout) :: triplet
        integer, intent(in) :: row
        
        integer :: k
        
        ! Zero out all entries in the specified row
        do k = 1, triplet%nnz
            if (triplet%rows(k) == row) then
                triplet%values(k) = 0.0_dp
            end if
        end do
        
    end subroutine zero_row
    
    subroutine solver_deallocate(this)
        class(poisson_1d_sparse_solver_t), intent(inout) :: this
        
        call this%mesh%deallocate()
        call this%K_triplet%destroy()
        call this%K_csr%destroy()
        if (allocated(this%f)) deallocate(this%f)
        if (allocated(this%solver)) then
            call this%solver%destroy()
            deallocate(this%solver)
        end if
        
    end subroutine solver_deallocate

end module fortfem_poisson_1d_sparse
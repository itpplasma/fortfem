module fortfem_sparse_matrix
    use fortfem_kinds
    implicit none
    private
    
    ! Compressed Sparse Row (CSR) format
    type, public :: csr_matrix_t
        integer :: n = 0           ! Matrix dimension
        integer :: nnz = 0         ! Number of non-zeros
        integer, allocatable :: row_ptr(:)   ! Row pointers (size n+1)
        integer, allocatable :: col_idx(:)   ! Column indices (size nnz)
        real(dp), allocatable :: values(:)   ! Non-zero values (size nnz)
    contains
        procedure :: init => csr_init
        procedure :: destroy => csr_destroy
        procedure :: get => csr_get
        procedure :: matvec => csr_matvec
    end type csr_matrix_t
    
    ! Triplet format for easy assembly
    type, public :: triplet_matrix_t
        integer :: n = 0           ! Matrix dimension
        integer :: nnz = 0         ! Current number of entries
        integer :: max_nnz = 0     ! Maximum capacity
        integer, allocatable :: rows(:)      ! Row indices
        integer, allocatable :: cols(:)      ! Column indices
        real(dp), allocatable :: values(:)   ! Values
    contains
        procedure :: init => triplet_init
        procedure :: destroy => triplet_destroy
        procedure :: add => triplet_add
        procedure :: to_csr => triplet_to_csr
    end type triplet_matrix_t
    
contains

    ! CSR matrix methods
    
    subroutine csr_init(this, n, nnz)
        class(csr_matrix_t), intent(out) :: this
        integer, intent(in) :: n, nnz
        
        this%n = n
        this%nnz = nnz
        allocate(this%row_ptr(n+1))
        allocate(this%col_idx(nnz))
        allocate(this%values(nnz))
        
        ! Initialize to zero
        this%row_ptr = 0
        this%col_idx = 0
        this%values = 0.0_dp
        
    end subroutine csr_init
    
    subroutine csr_destroy(this)
        class(csr_matrix_t), intent(inout) :: this
        
        if (allocated(this%row_ptr)) deallocate(this%row_ptr)
        if (allocated(this%col_idx)) deallocate(this%col_idx)
        if (allocated(this%values)) deallocate(this%values)
        this%n = 0
        this%nnz = 0
        
    end subroutine csr_destroy
    
    function csr_get(this, i, j) result(val)
        class(csr_matrix_t), intent(in) :: this
        integer, intent(in) :: i, j
        real(dp) :: val
        
        integer :: k
        
        val = 0.0_dp
        
        ! Search in row i
        do k = this%row_ptr(i), this%row_ptr(i+1)-1
            if (this%col_idx(k) == j) then
                val = this%values(k)
                return
            end if
        end do
        
    end function csr_get
    
    subroutine csr_matvec(this, x, y)
        class(csr_matrix_t), intent(in) :: this
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        
        integer :: i, k
        
        y = 0.0_dp
        
        do i = 1, this%n
            do k = this%row_ptr(i), this%row_ptr(i+1)-1
                y(i) = y(i) + this%values(k) * x(this%col_idx(k))
            end do
        end do
        
    end subroutine csr_matvec
    
    ! Triplet matrix methods
    
    subroutine triplet_init(this, n, max_nnz)
        class(triplet_matrix_t), intent(out) :: this
        integer, intent(in) :: n, max_nnz
        
        this%n = n
        this%nnz = 0
        this%max_nnz = max_nnz
        
        allocate(this%rows(max_nnz))
        allocate(this%cols(max_nnz))
        allocate(this%values(max_nnz))
        
    end subroutine triplet_init
    
    subroutine triplet_destroy(this)
        class(triplet_matrix_t), intent(inout) :: this
        
        if (allocated(this%rows)) deallocate(this%rows)
        if (allocated(this%cols)) deallocate(this%cols)
        if (allocated(this%values)) deallocate(this%values)
        this%n = 0
        this%nnz = 0
        this%max_nnz = 0
        
    end subroutine triplet_destroy
    
    subroutine triplet_add(this, i, j, val)
        class(triplet_matrix_t), intent(inout) :: this
        integer, intent(in) :: i, j
        real(dp), intent(in) :: val
        
        integer :: k
        
        ! Check if entry already exists
        do k = 1, this%nnz
            if (this%rows(k) == i .and. this%cols(k) == j) then
                ! Accumulate
                this%values(k) = this%values(k) + val
                return
            end if
        end do
        
        ! Add new entry
        if (this%nnz >= this%max_nnz) then
            error stop "Triplet matrix overflow"
        end if
        
        this%nnz = this%nnz + 1
        this%rows(this%nnz) = i
        this%cols(this%nnz) = j
        this%values(this%nnz) = val
        
    end subroutine triplet_add
    
    subroutine triplet_to_csr(this, csr)
        class(triplet_matrix_t), intent(in) :: this
        type(csr_matrix_t), intent(out) :: csr
        
        integer :: i, k, current_pos
        integer, allocatable :: row_count(:)
        integer, allocatable :: idx(:)
        
        ! Initialize CSR matrix
        call csr%init(this%n, this%nnz)
        
        ! Count entries per row
        allocate(row_count(this%n))
        row_count = 0
        
        do k = 1, this%nnz
            row_count(this%rows(k)) = row_count(this%rows(k)) + 1
        end do
        
        ! Build row pointers
        csr%row_ptr(1) = 1
        do i = 1, this%n
            csr%row_ptr(i+1) = csr%row_ptr(i) + row_count(i)
        end do
        
        ! Copy entries (simple insertion, not sorted)
        allocate(idx(this%n))
        idx = csr%row_ptr(1:this%n)
        
        do k = 1, this%nnz
            i = this%rows(k)
            current_pos = idx(i)
            csr%col_idx(current_pos) = this%cols(k)
            csr%values(current_pos) = this%values(k)
            idx(i) = idx(i) + 1
        end do
        
        deallocate(row_count, idx)
        
    end subroutine triplet_to_csr

end module fortfem_sparse_matrix
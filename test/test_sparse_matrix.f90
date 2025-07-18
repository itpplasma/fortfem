program test_sparse_matrix
    use fortfem_kinds
    use fortfem_sparse_matrix
    implicit none
    
    integer :: n_tests_passed = 0
    integer :: n_tests_failed = 0
    
    ! Test 1: Create and manipulate CSR matrix
    call test_csr_creation()
    
    ! Test 2: Add entries to sparse matrix
    call test_add_entries()
    
    ! Test 3: Convert from triplet to CSR format
    call test_triplet_to_csr()
    
    ! Test 4: Matrix-vector multiplication
    call test_matvec()
    
    ! Summary
    print *, "Tests passed: ", n_tests_passed
    print *, "Tests failed: ", n_tests_failed
    
    if (n_tests_failed > 0) then
        error stop "Some tests failed"
    end if
    
contains

    subroutine test_csr_creation()
        type(csr_matrix_t) :: A
        real(dp), parameter :: tol = 1.0e-14_dp
        
        ! Create 3x3 sparse matrix with pattern:
        ! [2 -1  0]
        ! [-1 2 -1]
        ! [0 -1  2]
        call A%init(n=3, nnz=7)
        
        ! Set CSR structure manually
        A%row_ptr = [1, 3, 6, 8]
        A%col_idx = [1, 2, 1, 2, 3, 2, 3]
        A%values = [2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp]
        
        ! Check structure
        if (A%n == 3 .and. A%nnz == 7) then
            print *, "PASS: CSR matrix dimensions"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: CSR matrix dimensions"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check diagonal element
        if (abs(A%get(2,2) - 2.0_dp) < tol) then
            print *, "PASS: CSR get element (2,2)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: CSR get element (2,2)"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call A%destroy()
        
    end subroutine test_csr_creation
    
    subroutine test_add_entries()
        type(triplet_matrix_t) :: T
        type(csr_matrix_t) :: A
        real(dp), parameter :: tol = 1.0e-14_dp
        
        ! Initialize triplet format
        call T%init(n=3, max_nnz=10)
        
        ! Add entries
        call T%add(1, 1, 2.0_dp)
        call T%add(1, 2, -1.0_dp)
        call T%add(2, 1, -1.0_dp)
        call T%add(2, 2, 2.0_dp)
        call T%add(2, 3, -1.0_dp)
        call T%add(3, 2, -1.0_dp)
        call T%add(3, 3, 2.0_dp)
        
        ! Add to existing entry
        call T%add(2, 2, 1.0_dp)  ! Should make (2,2) = 3.0
        
        ! Convert to CSR
        call T%to_csr(A)
        
        if (abs(A%get(2,2) - 3.0_dp) < tol) then
            print *, "PASS: Triplet add accumulation"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Triplet add accumulation"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call T%destroy()
        call A%destroy()
        
    end subroutine test_add_entries
    
    subroutine test_triplet_to_csr()
        type(triplet_matrix_t) :: T
        type(csr_matrix_t) :: A
        integer :: i
        
        ! Create identity matrix in triplet format
        call T%init(n=5, max_nnz=5)
        do i = 1, 5
            call T%add(i, i, 1.0_dp)
        end do
        
        ! Convert to CSR
        call T%to_csr(A)
        
        ! Check conversion
        if (A%nnz == 5) then
            print *, "PASS: Triplet to CSR nnz"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: Triplet to CSR nnz"
            n_tests_failed = n_tests_failed + 1
        end if
        
        ! Check row pointers for identity matrix
        if (all(A%row_ptr == [1, 2, 3, 4, 5, 6])) then
            print *, "PASS: CSR row pointers"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: CSR row pointers"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call T%destroy()
        call A%destroy()
        
    end subroutine test_triplet_to_csr
    
    subroutine test_matvec()
        type(csr_matrix_t) :: A
        real(dp) :: x(3), y(3), y_expected(3)
        real(dp), parameter :: tol = 1.0e-14_dp
        
        ! Create tridiagonal matrix
        call A%init(n=3, nnz=7)
        A%row_ptr = [1, 3, 6, 8]
        A%col_idx = [1, 2, 1, 2, 3, 2, 3]
        A%values = [2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp]
        
        ! Test vector
        x = [1.0_dp, 2.0_dp, 3.0_dp]
        
        ! Expected result: A*x
        y_expected = [0.0_dp, 0.0_dp, 4.0_dp]
        
        ! Compute matrix-vector product
        call A%matvec(x, y)
        
        if (maxval(abs(y - y_expected)) < tol) then
            print *, "PASS: CSR matrix-vector product"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "FAIL: CSR matrix-vector product"
            n_tests_failed = n_tests_failed + 1
        end if
        
        call A%destroy()
        
    end subroutine test_matvec

end program test_sparse_matrix
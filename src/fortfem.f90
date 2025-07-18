module fortfem
    ! Main module that re-exports all public interfaces
    use fortfem_kinds
    use fortfem_mesh_1d
    use fortfem_mesh_2d
    use fortfem_basis_1d
    use fortfem_basis_p1_2d
    use fortfem_assembly_1d
    use fortfem_assembly_2d
    use fortfem_poisson_1d
    use fortfem_poisson_2d
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    use fortfem_umfpack_interface
    use fortfem_poisson_1d_sparse
    
    implicit none
    
    ! Re-export everything
    public
    
end module fortfem

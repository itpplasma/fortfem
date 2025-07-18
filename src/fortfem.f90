module fortfem
    ! Main module that re-exports all public interfaces
    use fortfem_kinds
    use fortfem_mesh_1d
    use fortfem_basis_1d
    use fortfem_assembly_1d
    use fortfem_poisson_1d
    
    implicit none
    
    ! Re-export everything
    public
    
end module fortfem

module fortfem_assembly_2d
    ! Interface module for 2D assembly
    use assembly_2d_module, only: assembly_2d_t, assemble_laplacian, assemble_mass_matrix, assemble_mass_rhs
    implicit none
    public :: assembly_2d_t, assemble_laplacian, assemble_mass_matrix, assemble_mass_rhs
end module fortfem_assembly_2d
# FortFEM Test-Driven Development Plan

## Implementation Philosophy

Each feature follows strict TDD cycle:
1. Write failing test first
2. Implement minimal code to pass
3. Refactor while keeping tests green
4. Never write code without a failing test

## Phase 1: Foundation (Week 1-2)

### 1.1 Basic Types and Constants
- [ ] Test: Define real precision constants
- [ ] Test: Basic error handling module
- [ ] Test: Mathematical constants (pi, etc.)

### 1.2 Mesh Data Structure
- [ ] Test: Create mesh_t type with vertex storage
- [ ] Test: Add triangle connectivity
- [ ] Test: Compute number of edges from connectivity
- [ ] Test: Build edge-to-cell connectivity
- [ ] Test: Identify boundary edges
- [ ] Test: Read simple mesh format (start with hardcoded)

### 1.3 Reference Element
- [ ] Test: Define reference triangle coordinates
- [ ] Test: Define reference quadrilateral coordinates
- [ ] Test: Map from reference to physical element
- [ ] Test: Compute Jacobian of transformation

### 1.4 Quadrature Rules
- [ ] Test: 1-point quadrature on triangle
- [ ] Test: 3-point quadrature on triangle
- [ ] Test: Higher order quadrature rules
- [ ] Test: Quadrature on quadrilaterals
- [ ] Test: Verify exactness for polynomials

## Phase 2: P1 Elements (Week 3-4)

### 2.1 Lagrange P1 Basis Functions
- [ ] Test: Evaluate P1 basis on reference triangle
- [ ] Test: Verify partition of unity
- [ ] Test: Evaluate gradients on reference element
- [ ] Test: Transform gradients to physical element

### 2.2 Local Assembly
- [ ] Test: Mass matrix for single element
- [ ] Test: Stiffness matrix for single element
- [ ] Test: Load vector with constant source
- [ ] Test: Load vector with function source

### 2.3 Global Assembly
- [ ] Test: Assemble 2-element mesh mass matrix
- [ ] Test: Assemble 2-element mesh stiffness matrix
- [ ] Test: DOF numbering and connectivity
- [ ] Test: Assemble load vector

### 2.4 Boundary Conditions
- [ ] Test: Identify Dirichlet DOFs
- [ ] Test: Apply Dirichlet BC by elimination
- [ ] Test: Apply Dirichlet BC by penalty
- [ ] Test: Natural (Neumann) BC assembly

## Phase 3: First Complete Problem (Week 5)

### 3.1 Poisson Solver
- [ ] Test: Solve 1D Poisson with known solution
- [ ] Test: Solve 2D Poisson on unit square
- [ ] Test: Verify convergence order
- [ ] Test: Compare with analytical solution

### 3.2 Visualization
- [ ] Test: Export solution to array format
- [ ] Test: Create fortplotlib contour plot
- [ ] Test: Create mesh plot with solution
- [ ] Test: Animation for time-dependent problems

### 3.3 LAPACK Integration
- [ ] Test: Solve small linear system
- [ ] Test: Handle singular systems gracefully
- [ ] Test: Performance vs naive solver

## Phase 4: Higher Order Elements (Week 6-7)

### 4.1 P2 Lagrange Elements
- [ ] Test: P2 basis on reference element
- [ ] Test: Edge DOFs for P2
- [ ] Test: Assembly with P2 elements
- [ ] Test: Convergence rate verification

### 4.2 Quadrilateral Elements
- [ ] Test: Q1 bilinear basis functions
- [ ] Test: Tensor product quadrature
- [ ] Test: Mixed triangle/quad meshes

## Phase 5: Vector Elements (Week 8-9)

### 5.1 Raviart-Thomas RT0
- [ ] Test: RT0 basis on reference triangle
- [ ] Test: Divergence of RT0 basis
- [ ] Test: Piola transformation
- [ ] Test: Normal continuity across edges

### 5.2 Mixed Formulation
- [ ] Test: Saddle point system assembly
- [ ] Test: Darcy flow problem
- [ ] Test: Verify inf-sup condition

## Phase 6: Advanced Features (Week 10+)

### 6.1 Nedelec Elements
- [ ] Test: Edge element basis
- [ ] Test: Curl conformity
- [ ] Test: Maxwell equation example

### 6.2 DG Elements
- [ ] Test: Discontinuous P1 basis
- [ ] Test: Jump terms assembly
- [ ] Test: Upwind flux

### 6.3 Isoparametric Elements
- [ ] Test: Curved element geometry
- [ ] Test: High-order transformation
- [ ] Test: Integration accuracy

## Example Problems with Visualization

Each phase includes example problems that demonstrate features:

1. `heat_steady.f90` - Steady heat equation with fortplotlib contours
2. `heat_transient.f90` - Time-dependent heat with animation
3. `elasticity.f90` - Linear elasticity with deformation plot
4. `stokes.f90` - Stokes flow with velocity field
5. `darcy.f90` - Mixed formulation with flux vectors
6. `helmholtz.f90` - Wave propagation
7. `poisson_hp.f90` - hp-adaptivity demonstration

## Testing Infrastructure

### Unit Tests (`test/`)
- One test file per module
- Use FPM's built-in test runner
- Aim for >90% code coverage

### Integration Tests (`test/integration/`)
- Complete problem solutions
- Convergence tests
- Performance benchmarks

### Continuous Integration
- GitHub Actions on every push
- Test matrix: gfortran, ifort, nvfortran
- Codecov integration for coverage reports
- Automatic example generation

## Documentation (`doc/`)

### User Guide
- [ ] Installation instructions
- [ ] Tutorial: First FEM problem
- [ ] API reference (auto-generated)
- [ ] Mathematical background

### Developer Guide
- [ ] Architecture overview
- [ ] Adding new element types
- [ ] Testing guidelines
- [ ] Performance tips

## Success Metrics

- All tests pass before moving to next phase
- Code coverage > 90%
- Examples produce correct visualizations
- Performance within 2x of hand-optimized code
- Documentation for every public API
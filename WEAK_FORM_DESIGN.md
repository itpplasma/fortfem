# FortFEM Weak Form Framework Design

## Vision: Even More Elegant Than FreeFEM++

We want to create a Fortran framework that's more consistent and elegant than FreeFEM++ by leveraging:
- Fortran's powerful operator overloading
- Type-bound procedures for clean syntax
- Modern Fortran features for expressiveness

## Proposed Syntax

### Current FreeFEM++ Syntax
```cpp
// FreeFEM++ approach
varf a(u, v) = int2d(Th)( dx(u)*dx(v) + dy(u)*dy(v) );
varf L(v) = int2d(Th)( f*v );
matrix A = a(Vh, Vh);
real[int] b = L(0, Vh);
```

### Proposed FortFEM Syntax
```fortran
! FortFEM approach - more mathematical and consistent
use fortfem

! Mesh and function space
type(mesh_2d_t) :: mesh
type(function_space_t) :: V
type(trial_function_t) :: u
type(test_function_t) :: v

! Create mesh and space
mesh = square(10, 10)
V = P1(mesh)
call u%init(V)
call v%init(V)

! Define source function
f = constant(1.0_dp)

! Bilinear form: a(u,v) = ∫ ∇u·∇v dx
a = grad(u) .dot. grad(v)

! Linear form: L(v) = ∫ f·v dx  
L = f * v

! Alternative: combine forms naturally
problem = a - L  ! Weak form: a(u,v) = L(v)

! Solve with boundary conditions
u = solve(problem, dirichlet(u == 0.0_dp, boundary=all))
```

## Advanced Examples

### Elasticity
```fortran
! Vector-valued functions
type(vector_function_space_t) :: V
type(vector_trial_function_t) :: u
type(vector_test_function_t) :: v

V = P1_vector(mesh, dim=2)
call u%init(V)
call v%init(V)

! Material parameters
mu = 1.0_dp
lambda = 1.0_dp

! Bilinear form for elasticity
a = mu * (grad(u) .dot. grad(v)) + &
    (mu + lambda) * (div(u) * div(v))

! Body force
f = vector([0.0_dp, -1.0_dp])  ! Gravity
L = f .dot. v

! Solve
u = solve(a - L, dirichlet(u == 0.0_dp, boundary=bottom))
```

### Stokes Flow
```fortran
! Mixed function spaces
type(vector_function_space_t) :: V_vel
type(scalar_function_space_t) :: V_pres
type(mixed_function_space_t) :: V

V_vel = P2_vector(mesh, dim=2)
V_pres = P1(mesh)
V = mixed(V_vel, V_pres)

type(mixed_trial_function_t) :: up
type(mixed_test_function_t) :: vq

call up%init(V)
call vq%init(V)

! Extract components
u = up%component(1)  ! velocity
p = up%component(2)  ! pressure
v = vq%component(1)  ! velocity test
q = vq%component(2)  ! pressure test

! Stokes bilinear form
a = grad(u) .dot. grad(v) - p * div(v) - q * div(u)

! No RHS for this example
L = 0.0_dp

! Solve with boundary conditions
up = solve(a - L, [&
    dirichlet(u == [1.0_dp, 0.0_dp], boundary=top), &
    dirichlet(u == 0.0_dp, boundary=[left, right, bottom]) &
])
```

## Implementation Strategy

### 1. Core Types
- `trial_function_t` - represents unknown functions
- `test_function_t` - represents test functions  
- `bilinear_form_t` - represents bilinear forms a(u,v)
- `linear_form_t` - represents linear forms L(v)
- `weak_form_t` - represents complete weak formulation

### 2. Operator Overloading
```fortran
! Arithmetic operators
interface operator(+)
    module procedure add_bilinear_forms
    module procedure add_linear_forms
end interface

interface operator(*)
    module procedure scalar_times_function
    module procedure function_times_scalar
end interface

! Differential operators
interface grad
    module procedure grad_trial_function
    module procedure grad_test_function
end interface

! Inner products
interface operator(.dot.)
    module procedure dot_product_vectors
    module procedure dot_product_tensors
end interface
```

### 3. Integration Domains
```fortran
! Automatic integration over appropriate domains
a = grad(u) .dot. grad(v)  ! Automatically integrates over mesh
a = grad(u) .dot. grad(v) .over. interior(mesh)  ! Explicit domain
a = u * v .over. boundary(mesh, id=1)  ! Boundary integrals
```

### 4. Boundary Conditions
```fortran
! Elegant boundary condition syntax
bc = [&
    dirichlet(u == 0.0_dp, boundary=all), &
    neumann(grad(u) .dot. n == g, boundary=top), &
    robin(u + alpha * grad(u) .dot. n == beta, boundary=right) &
]
```

## Key Advantages Over FreeFEM++

1. **Type Safety**: Fortran's strong typing prevents many errors
2. **Consistency**: Same syntax for scalars, vectors, and mixed problems
3. **Expressiveness**: Natural mathematical notation
4. **Performance**: Compiled Fortran performance
5. **Memory Management**: Automatic with modern Fortran
6. **Extensibility**: Easy to add new element types and forms

## Implementation Phases

### Phase 1: Core Infrastructure
- Function space types
- Trial/test function types
- Basic operator overloading

### Phase 2: Differential Operators
- Gradient, divergence, curl
- Integration domains
- Bilinear/linear form assembly

### Phase 3: Boundary Conditions
- Dirichlet, Neumann, Robin
- Elegant syntax for BC specification

### Phase 4: Advanced Features
- Mixed function spaces
- Vector-valued problems
- Time-dependent problems

### Phase 5: Optimization
- Efficient assembly
- Sparse matrix optimization
- Parallel capabilities

This design would make FortFEM more elegant and mathematically natural than FreeFEM++, while maintaining the performance advantages of compiled Fortran.
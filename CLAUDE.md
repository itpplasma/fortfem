# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FortFEM is a modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS. It provides natural mathematical notation for defining weak forms and supports various element types on triangular and quadrilateral meshes.

## Build and Development Commands

### Build Commands
- **Build the project**: `fpm build`
- **Build with specific compiler**: `fpm build --compiler gfortran`
- **Build with coverage**: `fpm build --flag "--coverage"`
- **Clean build**: `fpm clean && fpm build`

### Run Commands
- **Run the main application**: `fpm run`
- **Run specific example**: `fpm run --example plot_basis`
- **List available examples**: `fpm run --example --list`

### Test Commands
- **Run all tests**: `fpm test`
- **Run specific test**: `fpm test test-name`
- **Run tests with coverage**: `fpm test --flag "--coverage"`

### Development Workflow
1. Follow strict TDD: Write test first, then implementation
2. Source files go in `src/`
3. Tests go in `test/`
4. Examples go in `example/`
5. Documentation goes in `doc/`
6. Check TODO.md for current implementation tasks

## Architecture and Code Structure

### Directory Layout
```
fortfem/
├── src/          # Core library modules
├── test/         # Unit and integration tests
├── example/      # Example programs with visualization
├── app/          # Main applications
└── doc/          # Documentation
```

### Module Organization (Planned)
- **src/fortfem.f90**: Main module exporting public API
- **src/mesh/**: Mesh data structures and I/O
- **src/elements/**: Finite element definitions
- **src/quadrature/**: Integration rules
- **src/forms/**: Weak form abstractions
- **src/assembly/**: Matrix and vector assembly
- **src/solvers/**: Linear solver interfaces

### Dependencies
- **fortplotlib**: For visualization in examples
- **LAPACK**: For linear algebra (initially)
- **SuiteSparse**: For sparse solvers (future)

### Fortran Conventions
- Derived types follow `typename_t` naming convention
- Line limit: 88 characters
- Indentation: 4 spaces
- Arrays extend with: `arr = [arr, new_element]` syntax
- Use associate blocks for unused dummy arguments
- Real precision: `dp = kind(0.0d0)`

### Testing Strategy
- Strict Test-Driven Development (TDD)
- Target >90% code coverage
- Tests must pass before implementation
- Use codecov for coverage tracking
- CI/CD runs on every push

### Example Development
- Examples use fortplotlib for visualization
- Each example demonstrates specific features
- Examples should be self-contained
- Include mathematical description in comments

### Continuous Integration
- GitHub Actions workflow in `.github/workflows/ci.yml`
- Tests on Ubuntu and macOS with gfortran
- Automatic codecov upload
- Example compilation verification
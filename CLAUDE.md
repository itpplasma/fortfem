# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fortfem is a Fortran-based project, likely intended for Finite Element Method implementations. It uses the Fortran Package Manager (FPM) as its build system.

## Build and Development Commands

### Build Commands
- **Build the project**: `fpm build`
- **Build with specific compiler**: `fpm build --compiler gfortran`
- **Clean build**: `fpm clean && fpm build`

### Run Commands
- **Run the main application**: `fpm run`
- **Run with specific executable**: `fpm run --target app-name`

### Test Commands
- **Run all tests**: `fpm test`
- **Run specific test**: `fpm test test-name`

### Development Workflow
1. For new features: Add source files to `src/`
2. For executables: Add programs to `app/`
3. For tests: Add test programs to `test/`
4. FPM automatically discovers and builds new files

## Architecture and Code Structure

### Module Organization
- **src/fortfem.f90**: Main module containing core functionality
  - Uses explicit typing (`implicit none`)
  - Follows private-by-default pattern with explicit public interfaces
  - Currently implements basic module structure ready for expansion

### Build Configuration (fpm.toml)
- Auto-builds executables, tests, and examples
- Enforces modern Fortran practices:
  - No implicit typing
  - No implicit externals
  - Free-form source

### Fortran Conventions
- Derived types should follow `typename_t` naming convention
- Line limit: 88 characters
- Indentation: 4 spaces
- Arrays extend with: `arr = [arr, new_element]` syntax
- Use associate blocks for unused dummy arguments to avoid warnings

### Testing Strategy
- Follow Test-Driven Development (TDD)
- Tests go in `test/` directory as separate programs
- Each test program should focus on a specific module or feature
- Use FPM's automatic test discovery
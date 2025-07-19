# FreeFEM Curl-Curl Benchmark

This directory contains a FreeFEM++ reference implementation for the 2D curl-curl equation, used to validate the FortFEM Nédélec element implementation.

## Problem Description

We solve the curl-curl equation:
```
curl(curl(E)) + E = f    in Ω = [0,1]²
E × n = 0                on ∂Ω
```

**Analytical Solution:**
- `E = [x*y, x²]`
- `curl(E) = ∂E_y/∂x - ∂E_x/∂y = 2x - x = x`
- `curl(curl(E)) = 0`
- `f = curl(curl(E)) + E = [x*y, x²]`

## Finite Element Method

- **Element Type:** RT0Ortho = Nédélec edge elements (not RT0!)
- **DOFs:** Edge-based degrees of freedom  
- **Properties:** H(curl) conforming, tangential continuity
- **Note:** RT0 is for H(div), RT0Ortho is for H(curl)

## Usage

### Prerequisites
Install FreeFEM++:
```bash
# Ubuntu/Debian
sudo apt install freefem++

# Arch Linux  
sudo pacman -S freefem++

# Or download from https://freefem.org
```

### Running the Benchmark

```bash
# Check FreeFEM installation
make check

# Run single benchmark case
make run

# Run full convergence study
make convergence

# Compare with FortFEM results (after running FortFEM)
make compare

# Generate comparison plots
make plot

# Show help
make help
```

### Output Files

- `freefem_results.dat` - Single case results (h, L2_error, Hcurl_error, DOFs)
- `freefem_convergence.dat` - Convergence study data
- `freefem_solution.dat` - Analytical solution values
- `freefem_fortfem_comparison.png` - Comparison plots

## Expected Results

For RT0 elements with the given analytical solution:

- **L2 Error:** Should converge as O(h) for smooth solutions
- **H(curl) Error:** Should converge as O(h) for the curl component
- **DOF Scaling:** Should scale as O(1/h²) for 2D problems

## Comparison with FortFEM

This benchmark provides a reference to validate:

1. **Correctness:** Do both implementations give similar error levels?
2. **Convergence Rates:** Do both achieve theoretical convergence rates?
3. **DOF Efficiency:** Do both use similar numbers of degrees of freedom?

## Technical Notes

### RT0 Elements in FreeFEM

FreeFEM's RT0 elements are equivalent to lowest-order Nédélec edge elements:
- Edge-based DOFs with tangential continuity
- Exact for polynomial fields up to order 0
- Natural enforcement of E × n = 0 boundary conditions

### Curl-Curl Formulation

The weak formulation is:
```
∫_Ω curl(E)·curl(v) dx + ∫_Ω E·v dx = ∫_Ω f·v dx
```

This corresponds to the bilinear form `a(E,v)` and linear form `l(v)` in the FreeFEM script.

## Troubleshooting

1. **FreeFEM not found:** Install FreeFEM++ using package manager or from source
2. **RT0 elements unavailable:** Load `Element_Mixte` library (included in script)
3. **Plotting issues:** Requires Python3 with matplotlib for comparison plots
4. **Large errors:** Check mesh refinement and boundary condition enforcement

## References

- FreeFEM++ documentation: https://freefem.org
- Nédélec elements: J.C. Nédélec, "Mixed finite elements in R³", 1980
- Raviart-Thomas elements: P.A. Raviart, J.M. Thomas, "A mixed finite element method", 1977
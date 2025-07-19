# Makefile for FreeFEM curl-curl benchmark
# Provides reference solution for FortFEM Nédélec element verification

FREEFEM = FreeFem++ -nw
SCRIPT = curl_curl_2d.edp

# Default target: run single case
.PHONY: all run convergence clean compare plot

all: run

# Run single benchmark case
run: $(SCRIPT)
	@echo "Running FreeFEM curl-curl benchmark..."
	$(FREEFEM) $(SCRIPT)
	@echo "FreeFEM benchmark completed."
	@echo "Results:"
	@cat freefem_results.dat

# Run convergence study
convergence: convergence_study.edp
	@echo "Running FreeFEM convergence study..."
	$(FREEFEM) convergence_study.edp
	@echo "Convergence study completed."
	@if [ -f freefem_convergence.dat ]; then \
		echo "Convergence results:"; \
		echo "h		L2_error	H(curl)_error	DOFs"; \
		tail -n +2 freefem_convergence.dat; \
	fi

# Create convergence study script
convergence_study.edp: $(SCRIPT)
	@echo "Creating convergence study script..."
	@sed 's/if (false)/if (true)/' $(SCRIPT) > convergence_study.edp

# Compare with FortFEM results (if available)
compare:
	@echo "Comparing FreeFEM and FortFEM results..."
	@if [ -f ../../../nedelec_convergence.dat ]; then \
		echo "FreeFEM vs FortFEM comparison:"; \
		echo "Source		h	L2_error	H(curl)_error	DOFs"; \
		echo "FreeFEM:"; \
		tail -n +2 freefem_results.dat | head -1 | awk '{print "FreeFEM	" $$0}'; \
		echo "FortFEM:"; \
		tail -n +2 ../../../nedelec_convergence.dat | head -1 | awk '{print "FortFEM	" $$0}'; \
	else \
		echo "FortFEM results not found. Run FortFEM convergence study first."; \
	fi

# Generate comparison plots (requires Python/matplotlib)
plot: freefem_convergence.dat
	@if command -v python3 >/dev/null 2>&1; then \
		echo "Generating comparison plots..."; \
		python3 plot_comparison.py; \
	else \
		echo "Python3 not available for plotting."; \
	fi

# Verify FreeFEM installation
check:
	@echo "Checking FreeFEM installation..."
	@if command -v $(FREEFEM) >/dev/null 2>&1; then \
		echo "FreeFEM found: $$($(FREEFEM) --version | head -1)"; \
		echo "RT0Ortho elements available: $$($(FREEFEM) -ne 'load \"Element_Mixte\"; cout << \"OK\" << endl;' 2>/dev/null || echo \"NOT FOUND\")"; \
	else \
		echo "FreeFEM not found. Please install FreeFEM++"; \
		echo "Ubuntu/Debian: sudo apt install freefem++"; \
		echo "Arch Linux: sudo pacman -S freefem++"; \
		echo "Or download from: https://freefem.org"; \
	fi

# Test constant field case
test_constant: test_constant_field.edp
	@echo "Testing RT0Ortho with constant field..."
	$(FREEFEM) test_constant_field.edp

# Verify weak form implementation
verify: verify_weak_form.edp
	@echo "Verifying curl-curl weak form..."
	$(FREEFEM) verify_weak_form.edp

# Direct comparison with FortFEM test cases
fortfem_compare: fortfem_comparison.edp
	@echo "Running FortFEM comparison tests..."
	$(FREEFEM) fortfem_comparison.edp

# Clean generated files
clean:
	rm -f *.dat *.eps *.png convergence_study.edp
	@echo "Cleaned benchmark files."

# Help
help:
	@echo "FreeFEM curl-curl benchmark Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  run        - Run single benchmark case"
	@echo "  convergence- Run convergence study"
	@echo "  compare    - Compare with FortFEM results"
	@echo "  plot       - Generate comparison plots"
	@echo "  check      - Verify FreeFEM installation"
	@echo "  clean      - Remove generated files"
	@echo "  help       - Show this help"
	@echo ""
	@echo "Usage:"
	@echo "  make run           # Single benchmark"
	@echo "  make convergence   # Full convergence study"
	@echo "  make compare       # Compare with FortFEM"

# Info about the benchmark
info:
	@echo "FreeFEM Curl-Curl Benchmark Information"
	@echo "======================================="
	@echo ""
	@echo "Problem: curl(curl(E)) + E = f in Ω"
	@echo "         E × n = 0 on ∂Ω"
	@echo ""
	@echo "Exact solution: E = [x*y, x^2]"
	@echo "curl(E) = x"
	@echo "curl(curl(E)) = 0"
	@echo "RHS: f = [x*y, x^2]"
	@echo ""
	@echo "Finite element: RT0 (Raviart-Thomas)"
	@echo "Mesh: Unit square [0,1]²"
	@echo ""
	@echo "Output files:"
	@echo "  freefem_results.dat     - Single case results"
	@echo "  freefem_convergence.dat - Convergence study"
	@echo "  freefem_solution.dat    - Analytical solution"
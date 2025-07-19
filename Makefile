# FortFEM Makefile
# Main targets for building, testing, and benchmarking

.PHONY: all build test benchmark clean help

# Default target
all: build

# Build the project
build:
	@echo "Building FortFEM..."
	fpm build

# Run all tests
test: build
	@echo "Running FortFEM tests..."
	fpm test

# Run benchmark comparison against FreeFEM
benchmark: build benchmark-fortfem benchmark-freefem compare-results

# Run FortFEM curl-curl benchmark
benchmark-fortfem:
	@echo "Running FortFEM curl-curl benchmark..."
	@fpm test test_curl_curl_system_solver > fortfem_benchmark.log 2>&1 || true
	@echo "FortFEM benchmark completed. Results saved to fortfem_benchmark.log"

# Run FreeFEM benchmark
benchmark-freefem:
	@echo "Running FreeFEM benchmark..."
	@cd benchmark/freefem/curl-curl && make run > ../../../freefem_benchmark.log 2>&1
	@echo "FreeFEM benchmark completed. Results saved to freefem_benchmark.log"

# Compare benchmark results
compare-results:
	@echo ""
	@echo "=== BENCHMARK COMPARISON ==="
	@echo ""
	@echo "FreeFEM Results:"
	@echo "=================="
	@grep -A 10 "=== FreeFEM Convergence Study ===" freefem_benchmark.log || echo "FreeFEM convergence data not found"
	@echo ""
	@echo "FortFEM Results:"
	@echo "=================="
	@grep -A 10 "h        L2_error        H(curl)_error   DOFs" fortfem_benchmark.log || echo "FortFEM convergence data not found"
	@echo ""
	@echo "=== END COMPARISON ==="

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	fpm clean
	@rm -f *.log
	@rm -f benchmark/freefem/curl-curl/*.dat

# Show help
help:
	@echo "FortFEM Makefile"
	@echo "=================="
	@echo ""
	@echo "Available targets:"
	@echo "  build      - Build the FortFEM project"
	@echo "  test       - Run all tests"
	@echo "  benchmark  - Run curl-curl benchmark comparison (FortFEM vs FreeFEM)"
	@echo "  clean      - Clean build artifacts and logs"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Benchmark comparison:"
	@echo "  This will run both FortFEM and FreeFEM curl-curl solvers"
	@echo "  with identical parameters and compare the convergence results."
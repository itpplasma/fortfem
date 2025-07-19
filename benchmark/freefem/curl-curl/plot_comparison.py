#!/usr/bin/env python3
"""
Comparison plotting script for FreeFEM vs FortFEM curl-curl benchmarks
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def load_data(filename):
    """Load convergence data from file"""
    try:
        data = np.loadtxt(filename)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except:
        return None

def plot_convergence_comparison():
    """Create comparison plots between FreeFEM and FortFEM"""
    
    # Load data
    freefem_data = load_data('freefem_convergence.dat')
    fortfem_data = load_data('../../../nedelec_convergence.dat')
    
    if freefem_data is None and fortfem_data is None:
        print("No convergence data found!")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # L2 error comparison
    if freefem_data is not None:
        h_ff = freefem_data[:, 0]
        L2_ff = freefem_data[:, 1]
        ax1.loglog(h_ff, L2_ff, 'o-', label='FreeFEM (RT0)', linewidth=2, markersize=6)
    
    if fortfem_data is not None:
        h_ft = fortfem_data[:, 0]
        L2_ft = fortfem_data[:, 1]
        ax1.loglog(h_ft, L2_ft, 's-', label='FortFEM (Nédélec)', linewidth=2, markersize=6)
    
    # Theoretical slopes
    if freefem_data is not None:
        h_theory = h_ff
        L2_theory = L2_ff[0] * (h_theory / h_ff[0])
        ax1.loglog(h_theory, L2_theory, '--', color='gray', label='O(h)', alpha=0.7)
    
    ax1.set_xlabel('Mesh size h')
    ax1.set_ylabel('L2 Error')
    ax1.set_title('L2 Error Convergence Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # H(curl) error comparison
    if freefem_data is not None:
        Hcurl_ff = freefem_data[:, 2]
        ax2.loglog(h_ff, np.maximum(Hcurl_ff, 1e-16), 'o-', label='FreeFEM (RT0)', linewidth=2, markersize=6)
    
    if fortfem_data is not None:
        Hcurl_ft = fortfem_data[:, 2]
        ax2.loglog(h_ft, np.maximum(Hcurl_ft, 1e-16), 's-', label='FortFEM (Nédélec)', linewidth=2, markersize=6)
    
    ax2.axhline(1e-14, color='red', linestyle=':', label='Machine precision', alpha=0.7)
    ax2.set_xlabel('Mesh size h')
    ax2.set_ylabel('H(curl) Error')
    ax2.set_title('H(curl) Error Convergence Comparison')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('freefem_fortfem_comparison.png', dpi=300, bbox_inches='tight')
    print("Comparison plot saved as freefem_fortfem_comparison.png")
    
    # Print numerical comparison
    print("\nNumerical Comparison:")
    print("=" * 60)
    print(f"{'Method':<12} {'h':<8} {'L2 Error':<12} {'H(curl) Error':<15} {'DOFs':<6}")
    print("-" * 60)
    
    if freefem_data is not None:
        for i in range(len(freefem_data)):
            h, L2, Hcurl, dofs = freefem_data[i]
            print(f"{'FreeFEM':<12} {h:<8.4f} {L2:<12.6e} {Hcurl:<15.6e} {int(dofs):<6}")
    
    if fortfem_data is not None:
        for i in range(len(fortfem_data)):
            h, L2, Hcurl, dofs = fortfem_data[i]
            print(f"{'FortFEM':<12} {h:<8.4f} {L2:<12.6e} {Hcurl:<15.6e} {int(dofs):<6}")

def plot_single_results():
    """Plot results from single benchmark runs"""
    
    freefem_single = load_data('freefem_results.dat')
    
    if freefem_single is not None:
        print("FreeFEM Single Case Results:")
        h, L2, Hcurl, dofs = freefem_single[0]
        print(f"  Mesh size: {h:.4f}")
        print(f"  L2 error: {L2:.6e}")
        print(f"  H(curl) error: {Hcurl:.6e}")
        print(f"  DOFs: {int(dofs)}")

if __name__ == "__main__":
    print("FreeFEM vs FortFEM Curl-Curl Benchmark Comparison")
    print("=" * 50)
    
    # Check for convergence data
    if os.path.exists('freefem_convergence.dat') or os.path.exists('../../../nedelec_convergence.dat'):
        plot_convergence_comparison()
    else:
        print("No convergence data found. Checking single case results...")
        plot_single_results()
    
    plt.show() if len(sys.argv) > 1 and sys.argv[1] == '--show' else None
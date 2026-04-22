#include "Metrics.hpp"
#include <iostream>
#include <iomanip>
#include <chrono> // For high-precision time measurements

// Convenient shortcut for time-related types
using namespace std::chrono;

void Metrics::run_experiment_1_single_system() {
    std::cout << "\n=== 4.1 Comparison of solving time for a single system ===\n";
    std::cout << "|   n   | Gauss (no pivot) | Gauss (pivot) | LU (Total) | LU (Decomp) | LU (Solve) |\n";
    std::cout << "|-------|------------------|---------------|------------|-------------|------------|\n";

    std::vector<int> sizes = { 100, 200, 500, 1000 };

    for (int n : sizes) {
        Matrix a = Matrix::generate_random(n, -1.0, 1.0, 42);
        Vector b = Matrix::generate_random_vector(n, -1.0, 1.0, 43);

        // 1. Gauss without pivoting
        auto start = high_resolution_clock::now();
        Solver::solve_gauss(a, b);
        auto end = high_resolution_clock::now();
        double time_gauss = duration<double>(end - start).count();

        // 2. Gauss with pivoting
        start = high_resolution_clock::now();
        Solver::solve_gauss_pivot(a, b);
        end = high_resolution_clock::now();
        double time_gauss_pivot = duration<double>(end - start).count();

        // 3. LU decomposition (Total = decomposition + substitution)
        start = high_resolution_clock::now();
        Solver::lu_result lu = Solver::decompose_lu(a);
        auto mid = high_resolution_clock::now(); // Time after decomposition
        Solver::solve_lu(lu, b);
        end = high_resolution_clock::now();      // Time after solving

        double time_lu_decomp = duration<double>(mid - start).count();
        double time_lu_solve = duration<double>(end - mid).count();
        double time_lu_total = time_lu_decomp + time_lu_solve;

        // Table row output
        std::cout << "| " << std::setw(5) << n << " | "
            << std::fixed << std::setprecision(6)
            << std::setw(16) << time_gauss << " | "
            << std::setw(13) << time_gauss_pivot << " | "
            << std::setw(10) << time_lu_total << " | "
            << std::setw(11) << time_lu_decomp << " | "
            << std::setw(10) << time_lu_solve << " |\n";
    }
}

void Metrics::run_experiment_2_multiple_rhs() {
    std::cout << "\n=== 4.2 Time savings for multiple right-hand sides ===\n";
    int n = 500;
    std::cout << "Matrix size n = " << n << "\n";
    std::cout << "|   k   | Gauss Pivot (k times) | LU (1 decomp + k solves) |\n";
    std::cout << "|-------|-----------------------|--------------------------|\n";

    Matrix a = Matrix::generate_random(n, -1.0, 1.0, 42);
    std::vector<int> k_values = { 1, 10, 100 };

    for (int k : k_values) {
        // Generate k random right-hand side vectors
        std::vector<Vector> b_vectors;
        for (int i = 0; i < k; ++i) {
            b_vectors.push_back(Matrix::generate_random_vector(n, -1.0, 1.0, 43 + i));
        }

        // Measure Gauss (solving from scratch for each vector)
        auto start = high_resolution_clock::now();
        for (int i = 0; i < k; ++i) {
            Solver::solve_gauss_pivot(a, b_vectors[i]);
        }
        auto end = high_resolution_clock::now();
        double time_gauss_k = duration<double>(end - start).count();

        // Measure LU (ONE decomposition, k solutions)
        start = high_resolution_clock::now();
        Solver::lu_result lu = Solver::decompose_lu(a);
        for (int i = 0; i < k; ++i) {
            Solver::solve_lu(lu, b_vectors[i]);
        }
        end = high_resolution_clock::now();
        double time_lu_k = duration<double>(end - start).count();

        std::cout << "| " << std::setw(5) << k << " | "
            << std::fixed << std::setprecision(6)
            << std::setw(21) << time_gauss_k << " | "
            << std::setw(24) << time_lu_k << " |\n";
    }
}

void Metrics::run_experiment_3_hilbert() {
    std::cout << "\n=== 4.3 Accuracy check on Hilbert matrices ===\n";
    std::vector<int> sizes = { 5, 10, 15 };

    for (int n : sizes) {
        std::cout << "\nHilbert matrix size n = " << n << "\n";
        Matrix h = Matrix::generate_hilbert(n);

        // Exact solution x = (1, 1, ..., 1)
        Vector x_exact(n, 1.0);

        // Compute right-hand side b = H * x
        Vector b = h.multiply(x_exact);
        double norm_x = calculate_norm(x_exact);

        // --- Gauss without pivoting ---
        try {
            Vector x_tilde_gauss = Solver::solve_gauss(h, b);
            Vector diff_gauss = subtract(x_tilde_gauss, x_exact);
            double rel_err_gauss = calculate_norm(diff_gauss) / norm_x;

            Vector b_tilde_gauss = h.multiply(x_tilde_gauss);
            double residual_gauss = calculate_norm(subtract(b_tilde_gauss, b));

            std::cout << "[Gauss no pivot] Rel. Error: " << std::scientific << rel_err_gauss
                << ", Residual: " << residual_gauss << "\n";
        }
        catch (const std::exception& e) {
            std::cout << "[Gauss no pivot] Error: " << e.what() << "\n";
        }

        // --- Gauss with pivoting ---
        try {
            Vector x_tilde_pivot = Solver::solve_gauss_pivot(h, b);
            Vector diff_pivot = subtract(x_tilde_pivot, x_exact);
            double rel_err_pivot = calculate_norm(diff_pivot) / norm_x;

            Vector b_tilde_pivot = h.multiply(x_tilde_pivot);
            double residual_pivot = calculate_norm(subtract(b_tilde_pivot, b));

            std::cout << "[Gauss pivot   ] Rel. Error: " << std::scientific << rel_err_pivot
                << ", Residual: " << residual_pivot << "\n";
        }
        catch (const std::exception& e) {
            std::cout << "[Gauss pivot   ] Error: " << e.what() << "\n";
        }
    }
}
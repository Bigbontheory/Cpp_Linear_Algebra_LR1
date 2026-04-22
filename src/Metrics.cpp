#include "Metrics.hpp"
#include <iostream>
#include <iomanip>
#include <chrono> 

using namespace std::chrono;

void Metrics::run_experiment_1_single_system() {
    std::cout << "\n=== 4.1 Comparison of solving time for a single system ===\n";
    // Теперь просто Gauss (который с пивотом по умолчанию) и LU
    std::cout << "|   n   |     Gauss      |  LU (Total)  |  LU (Decomp) |  LU (Solve) |\n";
    std::cout << "|-------|----------------|--------------|--------------|-------------|\n";

    std::vector<int> sizes = { 100, 200, 500, 1000 };

    for (int n : sizes) {
        Matrix a = Matrix::generate_random(n, -1.0, 1.0, 42);
        Vector b = Matrix::generate_random_vector(n, -1.0, 1.0, 43);

  
        auto start = high_resolution_clock::now();
        Solver::solve_gauss(a, b);
        auto end = high_resolution_clock::now();
        double time_gauss = duration<double>(end - start).count();

        start = high_resolution_clock::now();
        Solver::lu_result lu = Solver::decompose_lu(a);
        auto mid = high_resolution_clock::now(); 
        Solver::solve_lu(lu, b);
        end = high_resolution_clock::now();      

        double time_lu_decomp = duration<double>(mid - start).count();
        double time_lu_solve = duration<double>(end - mid).count();
        double time_lu_total = time_lu_decomp + time_lu_solve;

        std::cout << "| " << std::setw(5) << n << " | "
            << std::fixed << std::setprecision(6)
            << std::setw(14) << time_gauss << " | "
            << std::setw(12) << time_lu_total << " | "
            << std::setw(12) << time_lu_decomp << " | "
            << std::setw(11) << time_lu_solve << " |\n";
    }
}

void Metrics::run_experiment_2_multiple_rhs() {
    std::cout << "\n=== 4.2 Time savings for multiple right-hand sides ===\n";
    int n = 500;
    std::cout << "Matrix size n = " << n << "\n";
    std::cout << "|    k    |   Gauss (k times)   | LU (1 decomp + k solves) |\n";
    std::cout << "|---------|---------------------|--------------------------|\n";

    Matrix a = Matrix::generate_random(n, -1.0, 1.0, 42);
    std::vector<int> k_values = { 1, 10, 100 };

    for (int k : k_values) {
        std::vector<Vector> b_vectors;
        for (int i = 0; i < k; ++i) {
            b_vectors.push_back(Matrix::generate_random_vector(n, -1.0, 1.0, 43 + i));
        }

        auto start = high_resolution_clock::now();
        for (int i = 0; i < k; ++i) {
            Solver::solve_gauss(a, b_vectors[i]);
        }
        auto end = high_resolution_clock::now();
        double time_gauss_k = duration<double>(end - start).count();

        start = high_resolution_clock::now();
        Solver::lu_result lu = Solver::decompose_lu(a);
        for (int i = 0; i < k; ++i) {
            Solver::solve_lu(lu, b_vectors[i]);
        }
        end = high_resolution_clock::now();
        double time_lu_k = duration<double>(end - start).count();

        std::cout << "| " << std::setw(7) << k << " | "
            << std::fixed << std::setprecision(6)
            << std::setw(19) << time_gauss_k << " | "
            << std::setw(24) << time_lu_k << " |\n";
    }
}

void Metrics::run_experiment_3_hilbert() {
    std::cout << "\n=== 4.3 Accuracy check on Hilbert matrices ===\n";
    std::vector<int> sizes = { 5, 10, 15 };

    for (int n : sizes) {
        std::cout << "\nHilbert matrix size n = " << n << "\n";
        Matrix h = Matrix::generate_hilbert(n);

        Vector x_exact(n, 1.0);
        Vector b = h.multiply(x_exact);
        double norm_x = calculate_norm(x_exact);

        try {
            // Теперь вызываем просто solve_gauss
            Vector x_tilde = Solver::solve_gauss(h, b);
            Vector diff = subtract(x_tilde, x_exact);
            double rel_err = calculate_norm(diff) / norm_x;

            Vector b_tilde = h.multiply(x_tilde);
            double residual = calculate_norm(subtract(b_tilde, b));

            std::cout << "[Gauss Solver] Rel. Error: " << std::scientific << rel_err
                << ", Residual: " << residual << "\n";
        }
        catch (const std::exception& e) {
            std::cout << "[Gauss Solver] Error: " << e.what() << "\n";
        }
    }
}
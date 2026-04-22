#include "Metrics.hpp"
#include <iostream>

int main() {
    try {
        std::cout << "Lab Work: Direct Methods for Solving Linear Systems\n";
        std::cout << "==================================================\n";

        // Experiment 1: Comparison of execution time for different sizes of n
        Metrics::run_experiment_1_single_system();

        // Experiment 2: Time efficiency with multiple right-hand sides (k)
        Metrics::run_experiment_2_multiple_rhs();

        // Experiment 3: Accuracy testing on Hilbert matrices
        Metrics::run_experiment_3_hilbert();

        std::cout << "\n==================================================\n";
        std::cout << "Experiments completed successfully.\n";

    }
    catch (const std::exception& e) {
        // Catching division by zero or other runtime errors
        std::cerr << "\nCritical error during execution: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
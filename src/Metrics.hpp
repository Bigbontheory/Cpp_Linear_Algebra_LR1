#pragma once

#include "matrix.hpp"
#include "solver.hpp"

class Metrics {
public:
 
    static void run_experiment_1_single_system();

    static void run_experiment_2_multiple_rhs();

    static void run_experiment_3_hilbert();
};
#pragma once

#include "Matrix.hpp"

class Solver {
public:
    static Vector solve_gauss(Matrix a, Vector b);
    static Vector solve_gauss_pivot(Matrix a, Vector b);

    struct lu_result {
        Matrix l;
        Matrix u;
    };

    static lu_result decompose_lu(const Matrix& a);

    static Vector solve_lu(const lu_result& lu, const Vector& b);

private:
    static Vector back_substitution(const Matrix& u, const Vector& b);
    static Vector forward_substitution(const Matrix& l, const Vector& b);
};
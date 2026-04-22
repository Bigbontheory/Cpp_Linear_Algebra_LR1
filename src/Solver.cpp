#include "Solver.hpp"
#include <stdexcept>

Vector Solver::solve_gauss(Matrix a, Vector b) {
    int n = a.get_size();
  
    for (int k = 0; k < n; ++k) {
        double max_val = std::abs(a.at(k, k));
        int max_row = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(a.at(i, k)) > max_val) {
                max_val = std::abs(a.at(i, k));
                max_row = i;
            }
        }

        if (max_row != k) {
            for (int j = k; j < n; ++j) {
                std::swap(a.at(k, j), a.at(max_row, j));
            }
            std::swap(b[k], b[max_row]);
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = a.at(i, k) / a.at(k, k);
            for (int j = k; j < n; ++j) {
                a.at(i, j) -= factor * a.at(k, j);
            }
            b[i] -= factor * b[k];
        }
    }

    return back_substitution(a, b);
}

Solver::lu_result Solver::decompose_lu(const Matrix& a) {
    int n = a.get_size();
    Matrix l(n), u(n);

    for (int i = 0; i < n; i++) {
        l.at(i, i) = 1.0;

        for (int j = 0; j < n; j++) {
            double sum = 0;
            if (i <= j) {
                for (int k = 0; k < i; k++) {
                    sum += l.at(i, k) * u.at(k, j);
                }
                u.at(i, j) = a.at(i, j) - sum;
            }
            else {
                for (int k = 0; k < j; k++) {
                    sum += l.at(i, k) * u.at(k, j);
                }
                l.at(i, j) = (a.at(i, j) - sum) / u.at(j, j);
            }
        }
    }
    return { l, u };
}

Vector Solver::solve_lu(const lu_result& lu, const Vector& b) {
    Vector y = forward_substitution(lu.l, b);
    
    return back_substitution(lu.u, y);
}

Vector Solver::back_substitution(const Matrix& u, const Vector& b) {
    int n = u.get_size();
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += u.at(i, j) * x[j];
        }
        x[i] = (b[i] - sum) / u.at(i, i);
    }
    return x;
}

Vector Solver::forward_substitution(const Matrix& l, const Vector& b) {
    int n = l.get_size();
    Vector y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += l.at(i, j) * y[j];
        }
        y[i] = (b[i] - sum) / l.at(i, i);
    }
    return y;
}
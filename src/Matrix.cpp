#include "Matrix.hpp"

Matrix::Matrix(int size) : n(size), data(size, Vector(size, 0.0)) {}

Matrix::Matrix(const std::vector<std::vector<double>>& otherData)
    : n(static_cast<int>(otherData.size())), data(otherData) {
}

int Matrix::get_size() const {
    return n;
}

 
double& Matrix::at(int i, int j) {
    return data[i][j];
}

double Matrix::at(int i, int j) const {
    return data[i][j];
}


Vector Matrix::multiply(const Vector& x) const {
    Vector result(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += data[i][j] * x[j];
        }
    }
    return result;
}

Matrix Matrix::generate_random(int n, double min, double max, unsigned int seed) {
    Matrix mat(n); 
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(min, max);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat.at(i, j) = dis(gen);
        }
    }
    return mat;
}

Matrix Matrix::generate_hilbert(int n) {
    Matrix mat(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat.at(i, j) = 1.0 / (i + j + 1);
        }
    }
    return mat;
}

Vector Matrix::generate_random_vector(int n, double min, double max, unsigned int seed) {
    Vector b(n);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(min, max);

    for (int i = 0; i < n; ++i) {
        b[i] = dis(gen);
    }
    return b;
}

void Matrix::print(int precision) const {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << std::fixed << std::setprecision(precision)
                << std::setw(precision + 8) << data[i][j] << " ";
        }
        std::cout << "\n";
    }
}

double calculate_norm(const Vector& v) {
    double sum = 0;
    for (double x : v) sum += x * x;
    return std::sqrt(sum);
}

Vector subtract(const Vector& a, const Vector& b) {
    int size = static_cast<int>(a.size());
    Vector result(size);
    for (int i = 0; i < size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}
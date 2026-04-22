
#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>


using Vector = std::vector<double>;

class Matrix {
private:
    int n;
    std::vector<std::vector<double>> data;

public:
    Matrix(int size);
    Matrix(const std::vector<std::vector<double>>& otherData);

    int get_size() const;
    double& at(int i, int j);
    double at(int i, int j) const;

    Vector multiply(const Vector& x) const;

    static Matrix generate_random(int n, double min, double max, unsigned int seed = 42);
    static Matrix generate_hilbert(int n);
    static Vector generate_random_vector(int n, double min, double max, unsigned int seed = 43);

    void print(int precision = 4) const;

};

double calculate_norm(const Vector& v);
Vector subtract(const Vector& a, const Vector& b);
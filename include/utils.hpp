#pragma once

#include "consts.hpp"

#include <fstream>
#include <ostream>
#include <vector>
#include <unordered_set>

namespace NUtils {

inline TMatrix ReadMatrix(int n, const std::string& inputFile) {
    std::ifstream input(inputFile);
    TMatrix result;
    result.resize(n, TVector(n));
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            int x = 0;
            input >> x;
            result[row][col] = x;
        }
    }
    return result;
}

inline TVector ReadSizes(int amount, const std::string& inputFile) {
    std::ifstream input(inputFile);
    TVector result(amount);
    for (int i = 0; i < amount; ++i) {
        input >> result[i];
    }
    return result;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vector) {
    for (const auto& i : vector) {
        out << i << " ";
    }
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::unordered_set<T>& vector) {
    for (const auto& i : vector) {
        out << i << " ";
    }
    return out;
}

double CalculateLayoutCost(const TMatrix& matrix, const TVector& layout) {
    double cost = 0.0;
    int n = matrix.size();

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (matrix[i][j] > 0 && layout[i] != layout[j]) {
                cost += matrix[i][j];
            }
        }
    }

    return cost;
}

} // namespace NUtils

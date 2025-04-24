#include "consts.hpp"
#include "utils.hpp"
#include "sequentional.hpp"
#include "iterative.hpp"

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>

namespace {

void RunLayoutAlgorithms(const TMatrix& matrix, const TVector& sizes) {
    std::cout << "Matrix size: " << matrix.size() << "x" << matrix[0].size() << "\n";
    std::cout << "Component sizes: ";
    for (const auto& size : sizes) {
        std::cout << size << " ";
    }
    std::cout << "\n\n";

    TVector sequentialLayout = RunSequentialLayout(matrix, sizes);
    double sequentialCost = NUtils::CalculateLayoutCost(matrix, sequentialLayout);
    std::cout << "Sequential layout cost: " << sequentialCost << "\n\n";

    // Run iterative layout algorithm
    std::cout << "=== Iterative Layout Algorithm ===\n";
    TVector iterativeLayout = RunIterativeLayout(matrix, sizes);
    double iterativeCost = NUtils::CalculateLayoutCost(matrix, iterativeLayout);
    std::cout << "Iterative layout cost: " << iterativeCost << "\n";

    // Compare results
    std::cout << "\n=== Comparison ===\n";
    std::cout << "Sequential layout cost: " << sequentialCost << "\n";
    std::cout << "Iterative layout cost: " << iterativeCost << "\n";
    std::cout << "Improvement: " << (sequentialCost - iterativeCost) / sequentialCost * 100 << "%\n";
}

} // anonymous namespace

int main() {
    // Initialize random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    TMatrix matrix;
    TVector sizes;

    matrix = NUtils::ReadMatrix(SIZE, MATRIX_FILE);
    sizes = NUtils::ReadSizes(AMOUNT, SIZES_FILE);

    RunLayoutAlgorithms(matrix, sizes);

    return 0;
}

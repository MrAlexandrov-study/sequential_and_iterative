#include "consts.hpp"
#include "utils.hpp"
#include "sequentional.hpp"
#include "iterative.hpp"
#include "placement_sequential.hpp"
#include "placement_iterative.hpp"

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

    std::cout << "=== Iterative Layout Algorithm ===\n";
    TVector iterativeLayout = RunIterativeLayout(matrix, sizes);
    double iterativeCost = NUtils::CalculateLayoutCost(matrix, iterativeLayout);
    std::cout << "Iterative layout cost: " << iterativeCost << "\n";

    std::cout << "\n=== Comparison ===\n";
    std::cout << "Sequential layout cost: " << sequentialCost << "\n";
    std::cout << "Iterative layout cost: " << iterativeCost << "\n";
}

void RunPlacementAlgorithms(const TMatrix& matrix) {
    std::cout << "\n=== Placement Algorithms ===\n";
    std::cout << "Matrix size: " << matrix.size() << "x" << matrix[0].size() << "\n\n";

    // Run sequential placement
    std::cout << "=== Sequential Placement Algorithm ===\n";
    SequentialPlacement seqPlacement(matrix);
    auto initialGrid = seqPlacement.run();

    // Run iterative placement
    std::cout << "\n=== Iterative Placement Algorithm ===\n";
    IterativePlacement iterPlacement(matrix, initialGrid);
    auto finalGrid = iterPlacement.run();

    // Print final result
    std::cout << "\n=== Final Placement ===\n";
    for (size_t r = 0; r < finalGrid.size(); ++r) {
        for (size_t c = 0; c < finalGrid[r].size(); ++c) {
            std::cout << std::setw(3) << finalGrid[r][c] + 1;
        }
        std::cout << std::endl;
    }
}

} // anonymous namespace

int main() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    TMatrix matrix;
    TVector sizes;

    // Read matrix for layout algorithms
    matrix = NUtils::ReadMatrix(SIZE, MATRIX_FILE);
    sizes = NUtils::ReadSizes(AMOUNT, SIZES_FILE);

    // Run layout algorithms
    RunLayoutAlgorithms(matrix, sizes);

    // Create matrix for placement algorithms from Python example
    TMatrix placementMatrix = {
        {0, 1, 0, 0, 0, 1, 1, 0, 3, 1, 1, 1, 1, 3, 0, 0, 0, 0, 0, 1},
        {1, 0, 1, 0, 1, 2, 3, 0, 2, 1, 1, 0, 0, 0, 0, 0, 2, 1, 1, 1},
        {0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0},
        {0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 3, 0, 2, 0, 0, 1, 1, 2},
        {0, 1, 1, 3, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3},
        {1, 2, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 3},
        {1, 3, 0, 1, 0, 2, 0, 0, 1, 3, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0},
        {3, 2, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0},
        {1, 1, 2, 0, 0, 0, 3, 1, 0, 0, 0, 1, 0, 1, 0, 0, 4, 1, 0, 2},
        {1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0},
        {1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 3},
        {3, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 3, 1, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 2, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0},
        {0, 2, 0, 0, 0, 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 1, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {1, 1, 0, 2, 3, 3, 1, 0, 0, 2, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0},
    };

    // Run placement algorithms
    RunPlacementAlgorithms(placementMatrix);

    return 0;
}

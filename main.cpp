#include "consts.hpp"
#include "utils.hpp"
#include "layout_sequentional.hpp"
#include "layout_iterative.hpp"
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
    auto [iterativeLayout, iterativeCost] = RunIterativeLayout(matrix, sizes);
    std::cout << "Iterative layout cost: " << iterativeCost << "\n";

    std::cout << "\n=== Comparison ===\n";
    std::cout << "Sequential layout cost: " << sequentialCost << "\n";
    std::cout << "Iterative layout cost: " << iterativeCost << "\n";
}

void RunPlacementAlgorithms(const TMatrix& matrix) {
    std::cout << "\n=== Placement Algorithms ===\n";
    std::cout << "Matrix size: " << matrix.size() << "x" << matrix[0].size() << "\n\n";

    std::cout << "=== Sequential Placement Algorithm ===\n";
    SequentialPlacement seqPlacement(matrix);
    auto initialGrid = seqPlacement.run();

    std::cout << "\n=== Iterative Placement Algorithm ===\n";
    IterativePlacement iterPlacement(matrix, initialGrid);
    auto finalGrid = iterPlacement.run();

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

    matrix = NUtils::ReadMatrix(SIZE, LAYOUT_MATRIX_FILE);
    sizes = NUtils::ReadSizes(AMOUNT, SIZES_FILE);

    RunLayoutAlgorithms(matrix, sizes);

    // TMatrix placementMatrix = NUtils::ReadMatrix(SIZE, PLACEMENT_MATRIX_FILE);

    // RunPlacementAlgorithms(placementMatrix);

    return 0;
}

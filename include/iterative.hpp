#pragma once

#include "consts.hpp"
#include "utils.hpp"

#include "sequentional.hpp"

#include <cmath>
#include <iostream>

double CalculateMoveCost(const TMatrix& matrix, const TVector& layout, int vertex, int newComponent) {
    double costDiff = 0;
    int currentComponent = layout[vertex];
    int n = matrix.size();

    for (int i = 0; i < n; ++i) {
        if (i == vertex) continue;

        if (matrix[vertex][i] > 0) {
            if (layout[i] == currentComponent) {
                costDiff += matrix[vertex][i];
            } else if (layout[i] == newComponent) {
                costDiff -= matrix[vertex][i];
            }
        }
    }

    return costDiff;
}

TVector RunIterativeLayout(const TMatrix& matrix, const TVector& sizes, int maxIterations = 1000, bool printOutput = true) {
    if (printOutput) {
        std::cout << "Iterative Layout Algorithm\n";
        std::cout << "Sizes: ";
        for (const auto& size : sizes) {
            std::cout << size << " ";
        }
        std::cout << "\n";
    }

    TVector bestLayout = RunSequentialLayout(matrix, sizes, false);
    double bestCost = NUtils::CalculateLayoutCost(matrix, bestLayout);

    if (printOutput) {
        std::cout << "Initial layout cost: " << bestCost << "\n";
    }

    TVector currentLayout = bestLayout;

    std::unordered_map<int, int> componentSizes;
    for (int i = 0; i < sizes.size(); ++i) {
        componentSizes[i + 1] = sizes[i];
    }

    // Count vertices in each component
    std::unordered_map<int, int> componentCounts;
    for (const auto& component : currentLayout) {
        componentCounts[component]++;
    }

    // Iterative improvement
    int iterations = 0;
    int noImprovementCount = 0;
    const int maxNoImprovement = 10; // Allow some iterations without improvement

    while (noImprovementCount < maxNoImprovement && iterations < maxIterations) {
        iterations++;
        bool improvedThisIteration = false;

        // Try swapping pairs of vertices between different components
        for (int v1 = 0; v1 < matrix.size(); ++v1) {
            int comp1 = currentLayout[v1];

            for (int v2 = 0; v2 < matrix.size(); ++v2) {
                int comp2 = currentLayout[v2];

                // Only consider vertices in different components
                if (comp1 == comp2) continue;

                // Try swapping v1 and v2
                currentLayout[v1] = comp2;
                currentLayout[v2] = comp1;

                double newCost = NUtils::CalculateLayoutCost(matrix, currentLayout);

                if (newCost < bestCost) {
                    bestCost = newCost;
                    bestLayout = currentLayout;
                    improvedThisIteration = true;

                    if (printOutput) {
                        std::cout << "Iteration " << iterations
                                  << ": Swapped vertex " << v1 << " (component " << comp1
                                  << ") with vertex " << v2 << " (component " << comp2
                                  << "), new cost: " << bestCost << "\n";
                    }
                } else {
                    // Revert the swap
                    currentLayout[v1] = comp1;
                    currentLayout[v2] = comp2;
                }
            }
        }

        // Try moving individual vertices
        for (int vertex = 0; vertex < matrix.size(); ++vertex) {
            int currentComponent = currentLayout[vertex];

            // Try each other component
            for (int newComponent = 1; newComponent <= sizes.size(); ++newComponent) {
                if (newComponent == currentComponent) {
                    continue;
                }

                // Skip if this would make the current component too small
                if (componentCounts[currentComponent] <= componentSizes[currentComponent]) {
                    continue;
                }

                // Skip if the target component is already full
                if (componentCounts[newComponent] >= componentSizes[newComponent]) {
                    continue;
                }

                // Calculate cost change without actually moving
                double costDiff = CalculateMoveCost(matrix, currentLayout, vertex, newComponent);

                if (costDiff < 0) {  // Only make the move if it improves the solution
                    // Make the move
                    currentLayout[vertex] = newComponent;
                    componentCounts[currentComponent]--;
                    componentCounts[newComponent]++;

                    double newCost = bestCost + costDiff;
                    bestCost = newCost;
                    bestLayout = currentLayout;
                    improvedThisIteration = true;

                    if (printOutput) {
                        std::cout << "Iteration " << iterations
                                  << ": Moved vertex " << vertex
                                  << " from component " << currentComponent
                                  << " to component " << newComponent
                                  << ", new cost: " << bestCost << "\n";
                    }

                    break;  // Move to next vertex
                }
            }
        }

        // Try random moves to escape local optima
        if (!improvedThisIteration && noImprovementCount > maxNoImprovement / 2) {
            // Make a random perturbation to escape local optima
            for (int attempt = 0; attempt < 5; ++attempt) {
                int v1 = rand() % matrix.size();
                int v2 = rand() % matrix.size();

                int comp1 = currentLayout[v1];
                int comp2 = currentLayout[v2];

                if (comp1 == comp2) continue;

                // Check if swap would violate size constraints
                if (componentCounts[comp1] <= componentSizes[comp1] ||
                    componentCounts[comp2] <= componentSizes[comp2]) {
                    continue;
                }

                // Make the swap
                currentLayout[v1] = comp2;
                currentLayout[v2] = comp1;

                double newCost = NUtils::CalculateLayoutCost(matrix, currentLayout);

                if (newCost < bestCost) {
                    bestCost = newCost;
                    bestLayout = currentLayout;
                    improvedThisIteration = true;

                    if (printOutput) {
                        std::cout << "Random move in iteration " << iterations
                                  << " improved solution, new cost: " << bestCost << "\n";
                    }

                    break;
                } else {
                    // Keep the move anyway with some probability to escape local optima
                    double p = exp(-(newCost - bestCost) / (iterations * 0.1));
                    if ((double)rand() / RAND_MAX < p) {
                        if (printOutput) {
                            std::cout << "Random move accepted in iteration " << iterations
                                      << " despite worse cost: " << newCost << "\n";
                        }
                        break;
                    } else {
                        // Revert the swap
                        currentLayout[v1] = comp1;
                        currentLayout[v2] = comp2;
                    }
                }
            }
        }

        if (improvedThisIteration) {
            noImprovementCount = 0;
        } else {
            noImprovementCount++;
        }
    }

    if (printOutput) {
        std::cout << "Final layout cost after " << iterations << " iterations: " << bestCost << "\n";
        std::cout << "\nFinal component assignments:\n";
        for (int i = 0; i < matrix.size(); ++i) {
            std::cout << "Vertex " << i << ": Component " << bestLayout[i] << "\n";
        }
    }

    return bestLayout;
}

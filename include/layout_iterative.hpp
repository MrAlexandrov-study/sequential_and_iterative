#pragma once

#include "consts.hpp"
#include "utils.hpp"

#include "layout_sequentional.hpp"

#include <cmath>
#include <iostream>

double CalculateMoveCost(const TMatrix& matrix, const TVector& layout, int vertex, int newComponent) {
    double costDiff = 0;
    int currentComponent = layout[vertex];
    int n = matrix.size();


    for (int j = vertex + 1; j < n; ++j) {
        if (matrix[vertex][j] > 0) {
            if (layout[j] == currentComponent) {

                costDiff += matrix[vertex][j];
            } else if (layout[j] == newComponent) {

                costDiff -= matrix[vertex][j];
            }
        }
    }


    for (int i = 0; i < vertex; ++i) {
        if (matrix[i][vertex] > 0) {
            if (layout[i] == currentComponent) {

                costDiff += matrix[i][vertex];
            } else if (layout[i] == newComponent) {

                costDiff -= matrix[i][vertex];
            }
        }
    }

    return costDiff;
}

void InitializeLayout(
    const TMatrix& matrix,
    const TVector& sizes,
    TVector& currentLayout,
    TVector& bestLayout,
    double& bestCost,
    std::unordered_map<int, int>& componentSizes,
    std::unordered_map<int, int>& componentCounts,
    bool printOutput
) {

    bestLayout = RunSequentialLayout(matrix, sizes, false);
    bestCost = NUtils::CalculateLayoutCost(matrix, bestLayout);

    if (printOutput) {
        std::cout << "Initial layout cost: " << bestCost << "\n";
    }

    currentLayout = bestLayout;


    for (int i = 0; i < sizes.size(); ++i) {
        componentSizes[i + 1] = sizes[i];
    }


    componentCounts.clear();
    for (const auto& component : currentLayout) {
        componentCounts[component]++;
    }


    for (int i = 1; i <= sizes.size(); ++i) {
        if (componentCounts[i] != componentSizes[i]) {
            std::cerr << "Error: Initial layout has incorrect component sizes. "
                      << "Component " << i << " has " << componentCounts[i]
                      << " vertices, but should have " << componentSizes[i] << std::endl;
        }
    }
}

bool TrySwapVertices(
    const TMatrix& matrix,
    TVector& currentLayout,
    TVector& bestLayout,
    double& bestCost,
    std::unordered_map<int, int>& componentCounts,
    std::unordered_map<int, int>& componentSizes,
    int iterations,
    bool printOutput
) {
    bool improved = false;

    for (int v1 = 0; v1 < matrix.size(); ++v1) {
        int comp1 = currentLayout[v1];

        for (int v2 = 0; v2 < matrix.size(); ++v2) {
            int comp2 = currentLayout[v2];

            if (comp1 == comp2) continue;

            currentLayout[v1] = comp2;
            currentLayout[v2] = comp1;

            double newCost = NUtils::CalculateLayoutCost(matrix, currentLayout);

            if (newCost < bestCost) {
                bestCost = newCost;
                bestLayout = currentLayout;
                improved = true;

                if (printOutput) {
                    std::cout << "Iteration " << iterations
                              << ": Swapped vertex " << v1 << " (component " << comp1
                              << ") with vertex " << v2 << " (component " << comp2
                              << "), new cost: " << bestCost << "\n";
                }
            } else {
                currentLayout[v1] = comp1;
                currentLayout[v2] = comp2;
            }
        }
    }

    return improved;
}

bool TryMoveVertices(
    const TMatrix& matrix,
    TVector& currentLayout,
    TVector& bestLayout,
    double& bestCost,
    std::unordered_map<int, int>& componentCounts,
    std::unordered_map<int, int>& componentSizes,
    int iterations,
    bool printOutput
) {
    bool improved = false;

    for (int vertex = 0; vertex < matrix.size(); ++vertex) {
        int currentComponent = currentLayout[vertex];

        for (int newComponent = 1; newComponent <= componentSizes.size(); ++newComponent) {
            if (newComponent == currentComponent) {
                continue;
            }

            if (componentCounts[currentComponent] <= componentSizes[currentComponent]) {
                continue;
            }

            if (componentCounts[newComponent] >= componentSizes[newComponent]) {
                continue;
            }

            double costDiff = CalculateMoveCost(matrix, currentLayout, vertex, newComponent);

            if (costDiff < 0) {
                currentLayout[vertex] = newComponent;
                componentCounts[currentComponent]--;
                componentCounts[newComponent]++;

                double newCost = bestCost + costDiff;
                bestCost = newCost;
                bestLayout = currentLayout;
                improved = true;

                if (printOutput) {
                    std::cout << "Iteration " << iterations
                              << ": Moved vertex " << vertex
                              << " from component " << currentComponent
                              << " to component " << newComponent
                              << ", new cost: " << bestCost << "\n";
                }

                break;
            }
        }
    }

    return improved;
}

bool TryRandomMoves(
    const TMatrix& matrix,
    TVector& currentLayout,
    TVector& bestLayout,
    double& bestCost,
    std::unordered_map<int, int>& componentCounts,
    std::unordered_map<int, int>& componentSizes,
    int iterations,
    bool printOutput
) {
    bool improved = false;

    for (int attempt = 0; attempt < 5; ++attempt) {
        int v1 = rand() % matrix.size();
        int v2 = rand() % matrix.size();

        int comp1 = currentLayout[v1];
        int comp2 = currentLayout[v2];

        if (comp1 == comp2) continue;

        currentLayout[v1] = comp2;
        currentLayout[v2] = comp1;

        double newCost = NUtils::CalculateLayoutCost(matrix, currentLayout);

        if (newCost < bestCost) {
            bestCost = newCost;
            bestLayout = currentLayout;
            improved = true;

            if (printOutput) {
                std::cout << "Random move in iteration " << iterations
                          << " improved solution, new cost: " << bestCost << "\n";
            }

            break;
        } else {
            double p = exp(-(newCost - bestCost) / (iterations * 0.1));
            if ((double)rand() / RAND_MAX < p) {
                if (printOutput) {
                    std::cout << "Random move accepted in iteration " << iterations
                              << " despite worse cost: " << newCost << "\n";
                }
                break;
            } else {
                currentLayout[v1] = comp1;
                currentLayout[v2] = comp2;
            }
        }
    }

    return improved;
}

void PrintFinalResults(
    const TMatrix& matrix,
    const TVector& bestLayout,
    double bestCost,
    int iterations
) {
    std::cout << "Final layout cost after " << iterations << " iterations: " << bestCost << "\n";

    std::unordered_map<int, std::vector<int>> componentToVertices;
    for (int i = 0; i < matrix.size(); ++i) {
        componentToVertices[bestLayout[i]].push_back(i);
    }

    for (int component = 1; component <= componentToVertices.size(); ++component) {
        if (componentToVertices.find(component) != componentToVertices.end()) {
            std::cout << "Component " << component << " vertices: ";

            std::sort(componentToVertices[component].begin(), componentToVertices[component].end());

            for (const auto& vertex : componentToVertices[component]) {
                std::cout << vertex << " ";
            }
            std::cout << "\n";
        }
    }
}

std::pair<TVector, double> RunIterativeLayout(const TMatrix& matrix, const TVector& sizes, int maxIterations = 1000, bool printOutput = true) {
    if (printOutput) {
        std::cout << "Iterative Layout Algorithm\n";
        std::cout << "Sizes: ";
        for (const auto& size : sizes) {
            std::cout << size << " ";
        }
        std::cout << "\n";
    }

    TVector bestLayout;
    TVector currentLayout;
    double bestCost = 0;
    std::unordered_map<int, int> componentSizes;
    std::unordered_map<int, int> componentCounts;


    InitializeLayout(
        matrix,
        sizes,
        currentLayout,
        bestLayout,
        bestCost,
        componentSizes,
        componentCounts,
        printOutput
    );


    TVector originalLayout = bestLayout;
    double originalCost = bestCost;

    int iterations = 0;
    int noImprovementCount = 0;
    const int maxNoImprovement = 10;

    while (noImprovementCount < maxNoImprovement && iterations < maxIterations) {
        iterations++;
        bool improvedThisIteration = false;


        improvedThisIteration |= TrySwapVertices(
            matrix,
            currentLayout,
            bestLayout,
            bestCost,
            componentCounts,
            componentSizes,
            iterations,
            printOutput
        );


        improvedThisIteration |= TryMoveVertices(
            matrix,
            currentLayout,
            bestLayout,
            bestCost,
            componentCounts,
            componentSizes,
            iterations,
            printOutput
        );


        if (!improvedThisIteration && noImprovementCount > maxNoImprovement / 2) {
            improvedThisIteration |= TryRandomMoves(
                matrix,
                currentLayout,
                bestLayout,
                bestCost,
                componentCounts,
                componentSizes,
                iterations,
                printOutput
            );
        }

        if (improvedThisIteration) {

            componentCounts.clear();
            for (const auto& component : bestLayout) {
                componentCounts[component]++;
            }


            currentLayout = bestLayout;

            noImprovementCount = 0;
        } else {
            noImprovementCount++;
        }
    }


    std::unordered_map<int, int> finalComponentCounts;
    for (const auto& component : bestLayout) {
        finalComponentCounts[component]++;
    }


    std::vector<int> actualSizes;
    for (int i = 1; i <= sizes.size(); ++i) {
        actualSizes.push_back(finalComponentCounts[i]);
    }


    std::vector<int> expectedSizes(sizes.begin(), sizes.end());


    std::sort(actualSizes.begin(), actualSizes.end());
    std::sort(expectedSizes.begin(), expectedSizes.end());


    bool sizesMismatch = (actualSizes != expectedSizes);

    if (sizesMismatch) {
        if (printOutput) {
            std::cerr << "Warning: Final layout has incorrect component sizes.\n";
            std::cerr << "Expected sizes (sorted): ";
            for (const auto& size : expectedSizes) {
                std::cerr << size << " ";
            }
            std::cerr << "\nActual sizes (sorted): ";
            for (const auto& size : actualSizes) {
                std::cerr << size << " ";
            }
            std::cerr << std::endl;
        }
    }


    if (sizesMismatch) {
        if (printOutput) {
            std::cout << "Reverting to original layout due to component size constraints.\n";
        }
        bestLayout = originalLayout;
        bestCost = originalCost;
    }

    if (printOutput) {
        PrintFinalResults(matrix, bestLayout, bestCost, iterations);
    }

    return {bestLayout, bestCost};
}

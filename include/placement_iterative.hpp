#pragma once

#include "consts.hpp"
#include "utils.hpp"
#include "placement_sequential.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <unordered_set>
#include <tuple>

// Iterative placement algorithm based on the Python implementation
class IterativePlacement {
private:
    TMatrix matrix;
    std::vector<std::vector<int>> grid;
    int rows;
    int cols;
    int n;
    TVector rho;

public:
    IterativePlacement(const TMatrix& matrix, const std::vector<std::vector<int>>& initialGrid)
        : matrix(matrix), grid(initialGrid), rows(initialGrid.size()), cols(initialGrid[0].size()), n(matrix.size()) {
        rho = calculateRho();
    }

    // Calculate row sums (rho) for the matrix
    TVector calculateRho() const {
        TVector rho(n, 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                rho[i] += matrix[i][j];
            }
        }
        return rho;
    }

    // Get position (row, col) of a vertex in the grid
    std::pair<int, int> getPosition(int vertex) const {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (grid[r][c] == vertex) {
                    return {r, c};
                }
            }
        }
        return {-1, -1}; // Not found
    }

    // Compute the total cost Q of the placement
    double computeQ() const {
        double Q = 0.0;

        // Create a map of vertex to position
        std::vector<std::pair<int, int>> positions(n);
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int vertex = grid[r][c];
                positions[vertex] = {r, c};
            }
        }

        // Calculate total cost
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int w = matrix[i][j];
                if (w > 0) {
                    auto [r1, c1] = positions[i];
                    auto [r2, c2] = positions[j];
                    int dist = std::abs(r1 - r2) + std::abs(c1 - c2);
                    Q += w * dist;
                }
            }
        }

        return Q / 2.0; // Divide by 2 because each edge is counted twice
    }

    // Compute Li for a vertex
    double computeLi(int i) const {
        auto [r_i, c_i] = getPosition(i);
        double total = 0.0;

        for (int j = 0; j < n; ++j) {
            int w = matrix[i][j];
            if (w > 0) {
                auto [r_j, c_j] = getPosition(j);
                int dist = std::abs(r_i - r_j) + std::abs(c_i - c_j);
                total += w * dist;
            }
        }

        return total / rho[i];
    }

    // Run the iterative placement algorithm
    std::vector<std::vector<int>> run(bool printOutput = true) {
        std::vector<std::vector<int>> currentGrid = grid;
        std::unordered_set<std::string> seen;

        // Convert grid to string for seen set
        auto gridToString = [](const std::vector<std::vector<int>>& g) {
            std::string result;
            for (const auto& row : g) {
                for (int val : row) {
                    result += std::to_string(val) + ",";
                }
            }
            return result;
        };

        seen.insert(gridToString(currentGrid));
        int iteration = 0;

        double Q_current = computeQ();
        if (printOutput) {
            std::cout << "\nStart Q = " << Q_current << std::endl;
        }

        while (true) {
            // Calculate L values for all vertices
            std::vector<double> L(n);
            for (int i = 0; i < n; ++i) {
                L[i] = computeLi(i);
            }

            // Find vertex with worst L value
            int m = std::distance(L.begin(), std::max_element(L.begin(), L.end()));
            double worst_L = L[m];

            // Find candidates for swap
            std::vector<std::tuple<int, double, double, double, std::vector<std::vector<int>>>> candidates;

            for (int k = 0; k < n; ++k) {
                if (k == m) continue;

                // Create a copy of the current grid
                auto g = currentGrid;

                // Get positions of vertices m and k
                auto [r_m, c_m] = getPosition(m);
                auto [r_k, c_k] = getPosition(k);

                // Swap vertices m and k
                std::swap(g[r_m][c_m], g[r_k][c_k]);

                // Check if this configuration has been seen before
                std::string gridStr = gridToString(g);
                if (seen.find(gridStr) != seen.end()) {
                    continue;
                }

                // Create a temporary placement object to compute new values
                IterativePlacement temp(matrix, g);

                // Calculate new L value for vertex m
                double new_L = temp.computeLi(m);
                double delta_L = new_L - worst_L;

                // Calculate new total cost
                double new_Q = temp.computeQ();

                // Add to candidates if it improves both L and Q
                if (delta_L < 0 && new_Q < Q_current) {
                    candidates.push_back(std::make_tuple(k, delta_L, new_L, new_Q, g));
                }
            }

            // If no improving candidates, we're done
            if (candidates.empty()) {
                if (printOutput) {
                    std::cout << "\nEnded" << std::endl;
                }
                break;
            }

            // Sort candidates by delta_L and then by vertex index
            std::sort(candidates.begin(), candidates.end(),
                [](const auto& a, const auto& b) {
                    if (std::get<1>(a) != std::get<1>(b)) {
                        return std::get<1>(a) < std::get<1>(b);
                    }
                    return std::get<0>(a) < std::get<0>(b);
                });

            // Select the best candidate
            auto [k_best, best_delta, L_new, Q_new, grid_new] = candidates[0];

            iteration++;
            if (printOutput) {
                std::cout << "\nIteration " << iteration << ":" << std::endl;
                std::cout << "  Worst edge v" << m + 1 << ", L_old=" << std::fixed << std::setprecision(3) << worst_L
                          << ", Q_old=" << std::fixed << std::setprecision(1) << Q_current << std::endl;
                std::cout << "  cadidats (v_k): delta L, L_new, Q_new:" << std::endl;

                for (const auto& [k, dL, Ln, Qn, _] : candidates) {
                    std::cout << "    v" << m + 1 << "<-> v" << k + 1 << ": delta L=" << std::fixed << std::setprecision(3) << dL
                              << ", L_new=" << std::fixed << std::setprecision(3) << Ln
                              << ", Q_new=" << std::fixed << std::setprecision(1) << Qn << std::endl;
                }

                std::cout << "  Choose change v" << m + 1 << "↔v" << k_best + 1
                          << ", delta L=" << std::fixed << std::setprecision(3) << best_delta
                          << ", L_new=" << std::fixed << std::setprecision(3) << L_new
                          << ", Q_new=" << std::fixed << std::setprecision(1) << Q_new << std::endl;

                std::cout << "  New variant " << rows << "×" << cols << ":" << std::endl;
                for (int r = 0; r < rows; ++r) {
                    for (int c = 0; c < cols; ++c) {
                        std::cout << std::setw(3) << grid_new[r][c] + 1;
                    }
                    std::cout << std::endl;
                }
            }

            // Update current grid and Q
            currentGrid = grid_new;
            seen.insert(gridToString(currentGrid));
            Q_current = Q_new;
        }

        return currentGrid;
    }
};

#pragma once

#include "consts.hpp"
#include "utils.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>

class SequentialPlacement {
private:
    TMatrix matrix;
    int rows;
    int cols;
    int n;

public:
    SequentialPlacement(const TMatrix& matrix, int rows = 4, int cols = 5)
        : matrix(matrix), rows(rows), cols(cols), n(matrix.size()) {}

    TVector calculateRho() const {
        TVector rho(n, 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                rho[i] += matrix[i][j];
            }
        }
        return rho;
    }

    void printStep(int step, const std::vector<int>& selected, const TVector& rho, const std::vector<double>& K) const {
        std::cout << "\nStep " << step << ", layouted [";
        for (size_t i = 0; i < selected.size(); ++i) {
            std::cout << "v" << selected[i] + 1;
            if (i < selected.size() - 1) std::cout << ", ";
        }
        std::cout << "]:" << std::endl;

        std::cout << std::setw(6) << " ";
        for (int j = 0; j < n; ++j) {
            std::cout << std::setw(6) << "v" << j + 1;
        }
        std::cout << std::setw(8) << "rho" << std::setw(8) << "K" << std::endl;

        for (int i = 0; i < n; ++i) {
            std::cout << "v" << i + 1 << std::setw(3) << " ";
            for (int j = 0; j < n; ++j) {
                std::cout << std::setw(6) << matrix[i][j];
            }
            std::cout << std::setw(8) << rho[i];

            if (std::find(selected.begin(), selected.end(), i) != selected.end()) {
                std::cout << std::setw(8) << "*";
            } else {
                std::cout << std::setw(8) << std::fixed << std::setprecision(1) << K[i];
            }
            std::cout << std::endl;
        }
    }

    std::vector<std::vector<int>> run(bool printOutput = true) {
        TVector rho = calculateRho();
        std::vector<int> selected = {0};
        std::vector<int> order = {0};

        for (int step = 1; step < n; ++step) {
            std::vector<double> K(n, std::numeric_limits<double>::quiet_NaN());

            for (int i = 0; i < n; ++i) {
                if (std::find(selected.begin(), selected.end(), i) != selected.end()) {
                    continue;
                }

                double P_sel = 0.0;
                for (int sel : selected) {
                    P_sel += matrix[i][sel];
                }

                K[i] = 2 * P_sel - rho[i];
            }

            if (printOutput) {
                printStep(step, selected, rho, K);
            }

            // Find maximum K value
            double maxK = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < n; ++i) {
                if (std::find(selected.begin(), selected.end(), i) == selected.end() &&
                    !std::isnan(K[i]) && K[i] > maxK) {
                    maxK = K[i];
                }
            }

            std::vector<int> candidates;
            for (int i = 0; i < n; ++i) {
                if (std::find(selected.begin(), selected.end(), i) == selected.end() &&
                    !std::isnan(K[i]) && std::abs(K[i] - maxK) < 1e-6) {
                    candidates.push_back(i);
                }
            }

            int next_v = *std::min_element(candidates.begin(), candidates.end());
            selected.push_back(next_v);
            order.push_back(next_v);
        }

        std::vector<std::vector<int>> grid(rows, std::vector<int>(cols));
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int idx = r * cols + c;
                grid[r][c] = order[idx];
            }
        }

        if (printOutput) {
            std::cout << "\nStart layout " << rows << "×" << cols << " (edge numbers):" << std::endl;
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    std::cout << std::setw(3) << grid[r][c] + 1;
                }
                std::cout << std::endl;
            }
        }

        return grid;
    }
};

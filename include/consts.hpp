#pragma once

#include <vector>

constexpr int SIZE = 20;
constexpr int AMOUNT = 4;

constexpr auto MATRIX_FILE = "matrix.txt";
constexpr auto SIZES_FILE = "sizes.txt";

using Type = int;
using TVector = std::vector<Type>;
using TMatrix = std::vector<TVector>;

#include "consts.hpp"
#include "utils.hpp"

// #include "project_template.hpp"
#include <algorithm>

#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include <cassert>
#include <fstream>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>

namespace {

TVector CountRho(const TMatrix& matrix) {
    int n = static_cast<int>(matrix.size());
    TVector result(n, 0);
    for (int row = 0; row < n; ++row) {
        result[row] = boost::accumulate(matrix[row], 0);
    }
    return result;
}

int CountDelta(
    const TVector& rho
    , const TMatrix& matrix
    , const std::unordered_set<Type>& newComponent
    , int vertex
) {
    Type sum = 0;
    for (const auto& to : newComponent) {
        sum += matrix[vertex][to];
    }
    return rho[vertex] - 2 * sum;
}

auto GetNeibours(const TMatrix& matrix, int vertex) {
    std::unordered_set<Type> result;
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        if (matrix[vertex][i] != 0) {
            result.insert(i);
        }
    }
    return result;
}

auto CountOuterConnections(
    const TMatrix& matrix
    , const std::unordered_set<int>& newComponent
    , int vertex
) {
    int result = 0;
    for (const auto& to : newComponent) {
        result += matrix[vertex][to];
    }
    return result;
}

auto RemoveVertexes(
    const TMatrix& matrix
    , const std::unordered_set<int>& newComponent
    , std::unordered_map<int, int>& localToGlobalNumber
) {
    int newSize = matrix.size() - newComponent.size();
    TMatrix newMatrix;
    newMatrix.resize(newSize, TVector(newSize));

    std::unordered_map<int, int> oldToNewIndex;
    int newIdx = 0;
    for (int i = 0, endi = matrix.size(); i < endi; ++i) {
        if (newComponent.contains(i)) {
            continue;
        }
        oldToNewIndex[i] = newIdx++;
    }

    std::unordered_map<int, int> updatedLocalToGlobal;
    for (const auto& [local, global] : localToGlobalNumber) {
        if (!newComponent.contains(local)) {
            updatedLocalToGlobal[oldToNewIndex[local]] = global;
        }
    }
    localToGlobalNumber = updatedLocalToGlobal;

    for (int i = 0, endi = matrix.size(); i < endi; ++i) {
        if (newComponent.contains(i)) {
            continue;
        }

        int newRow = oldToNewIndex[i];
        for (int j = 0, endj = matrix.size(); j < endj; ++j) {
            if (newComponent.contains(j)) {
                continue;
            }

            int newCol = oldToNewIndex[j];
            newMatrix[newRow][newCol] = matrix[i][j];
        }
    }

    return newMatrix;
}

TVector ReadSizes(int amount = 4, const std::string& inputFile = "sizes.txt") {
    std::ifstream input(inputFile);
    int size = 0;
    TVector result(amount);
    for (int i = 0; i < amount; ++i) {
        input >> result[i];
    }
    return result;
}

struct GraphState {
    TMatrix matrix;
    TVector color;
    std::unordered_map<int, int> localToGlobalNumber;
    std::unordered_set<int> assignedVertices;
    std::unordered_set<int> availableVertices;
    TVector sizes;
};

GraphState InitializeGraphState() {
    GraphState state;
    state.matrix = NUtils::ReadMatrix(SIZE, "matrix.txt");
    state.color.resize(SIZE, 0);
    state.sizes = ReadSizes();

    for (int i = 0; i < SIZE; ++i) {
        state.localToGlobalNumber[i] = i;
    }

    for (int i = 0; i < SIZE; ++i) {
        state.availableVertices.insert(i);
    }

    return state;
}

int FindMinimumElement(const TVector& rho) {
    auto minElement = 0;
    auto minValue = std::numeric_limits<int>::max();

    for (int i = 0; i < rho.size(); ++i) {
        if (rho[i] < minValue) {
            minValue = rho[i];
            minElement = i;
        }
    }

    return minElement;
}

std::unordered_set<int> FindConnectedElements(
    const TMatrix& matrix,
    const std::unordered_set<int>& component,
    int currentSize
) {
    std::unordered_set<int> connectedElements;

    for (int i = 0; i < currentSize; ++i) {
        if (component.contains(i)) continue;

        bool isConnected = false;
        for (const auto& vertex : component) {
            if (matrix[i][vertex] > 0) {
                isConnected = true;
                break;
            }
        }

        if (isConnected) {
            connectedElements.insert(i);
        }
    }

    return connectedElements;
}

std::optional<int> FindElementToAdd(
    const TMatrix& matrix,
    const std::unordered_set<int>& component,
    int currentSize
) {
    TVector outerConnections(currentSize, 0);
    std::optional<int> elementToAdd;

    for (int i = 0; i < currentSize; ++i) {
        if (component.contains(i)) continue;

        outerConnections[i] = CountOuterConnections(matrix, component, i);

        if (elementToAdd.has_value()) {
            if (outerConnections[i] < outerConnections[elementToAdd.value()]) {
                elementToAdd = i;
            }
        } else {
            elementToAdd = i;
        }
    }

    return elementToAdd;
}

void GrowComponentToTargetSize(
    std::unordered_set<int>& component,
    const TMatrix& matrix,
    int target,
    int currentSize,
    int& minElement
) {
    auto elementToAdd = FindElementToAdd(matrix, component, currentSize);

    if (!elementToAdd.has_value()) {
        return;
    }

    component.insert(elementToAdd.value());
    minElement = elementToAdd.value();

    if (component.size() == target) {
        return;
    }

    auto connectedElements = FindConnectedElements(matrix, component, currentSize);

    for (const auto& element : connectedElements) {
        if (component.size() >= target) break;
        component.insert(element);
    }
}

void ShrinkComponentToTargetSize(
    std::unordered_set<int>& component,
    std::unordered_set<int>& neighbors,
    const TMatrix& matrix,
    const TVector& rho,
    int target,
    int minElement
) {
    while (component.size() > target) {
        std::unordered_map<int, int> deltas;
        deltas.reserve(neighbors.size());

        for (const auto& vertex : neighbors) {
            deltas[vertex] = CountDelta(rho, matrix, component, vertex);
        }

        auto elementToDelete = deltas.cbegin();
        for (auto element = deltas.begin(), end = deltas.end(); element != end; ++element) {
            if (element->second > elementToDelete->second) {
                elementToDelete = element;
            }
        }

        neighbors.erase(elementToDelete->first);

        component = neighbors;
        component.insert(minElement);
    }
}

void UpdateAndPrintComponentInfo(
    const std::unordered_set<int>& component,
    const std::unordered_map<int, int>& localToGlobalNumber,
    int componentIndex,
    TVector& color,
    std::unordered_set<int>& assignedVertices,
    std::unordered_set<int>& availableVertices
) {
    std::cout << "Component " << componentIndex + 1 << " vertices: ";
    std::unordered_set<int> globalIndices;

    for (const auto& element : component) {
        auto realIndex = localToGlobalNumber.at(element);
        globalIndices.insert(realIndex);
        assignedVertices.insert(realIndex);
        color[realIndex] = componentIndex + 1;
    }

    std::vector<int> sortedIndices(globalIndices.begin(), globalIndices.end());
    std::ranges::sort(sortedIndices);
    for (const auto& index : sortedIndices) {
        std::cout << index << " ";
    }
    std::cout << "\n";

    for (const auto& element : component) {
        auto realIndex = localToGlobalNumber.at(element);
        availableVertices.erase(realIndex);
    }
}

void ProcessComponent(
    GraphState& state,
    int componentIndex
) {
    auto rho = CountRho(state.matrix);

    auto minElement = FindMinimumElement(rho);

    auto neighbors = GetNeibours(state.matrix, minElement);

    auto component = neighbors;
    component.insert(minElement);

    const int target = state.sizes[componentIndex];
    const int currentSize = state.matrix.size();

    while (component.size() != target) {
        if (component.size() < target) {
            GrowComponentToTargetSize(component, state.matrix, target, currentSize, minElement);
        } else {
            ShrinkComponentToTargetSize(component, neighbors, state.matrix, rho, target, minElement);
        }
    }

    UpdateAndPrintComponentInfo(
        component,
        state.localToGlobalNumber,
        componentIndex,
        state.color,
        state.assignedVertices,
        state.availableVertices
    );

    state.matrix = RemoveVertexes(state.matrix, component, state.localToGlobalNumber);
}

void PrintFinalAssignments(const TVector& color) {
    std::cout << "\nFinal component assignments:" << "\n";
    for (int i = 0; i < SIZE; ++i) {
        std::cout << "Vertex " << i << ": Component " << color[i] << "\n";
    }
}

} // anonymous namespace

int main() {
    using namespace NUtils;
    auto state = InitializeGraphState();

    std::cout << "sizes: " << state.sizes << "\n";

    for (int currentComponent = 0; currentComponent < state.sizes.size(); ++currentComponent) {
        ProcessComponent(state, currentComponent);
    }

    PrintFinalAssignments(state.color);
    return 0;
}

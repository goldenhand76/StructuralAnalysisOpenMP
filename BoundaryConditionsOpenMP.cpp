#include "BoundaryConditionsOpenMP.h"
#include <omp.h>

// Function to apply boundary conditions (fixed nodes and loads)
void BoundaryConditions::applyBoundaryConditions(std::vector<std::vector<double>>& K_global, std::vector<double>& F, const std::vector<int>& fixedNodes, const std::vector<Load>& loads) {
    // Apply boundary conditions for fixed nodes in parallel
#pragma omp parallel for
    for (int i = 0; i < fixedNodes.size(); ++i) {
        int node = fixedNodes[i];

        // Set entire row corresponding to the fixed node to 0
        for (int j = 0; j < K_global[node].size(); ++j) {
            K_global[node][j] = 0.0;
        }
        // Set the diagonal entry to 1 to maintain matrix invertibility
        K_global[node][node] = 1.0;

        // Set the corresponding force to 0 for fixed nodes
        F[node] = 0.0;
    }

    // Apply loads to the force vector in parallel
#pragma omp parallel for
    for (int i = 0; i < loads.size(); ++i) {
        const auto& load = loads[i];
        F[load.node] = load.value;
    }
}

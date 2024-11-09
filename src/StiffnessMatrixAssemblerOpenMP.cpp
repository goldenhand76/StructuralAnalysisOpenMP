#include "StiffnessMatrixAssemblerOpenMP.h"
#include <iostream>
#include <omp.h>

// Helper function to compute the stiffness matrix for an individual element
void StiffnessMatrixAssembler::computeElementStiffness(const Triangle& element, double E, double nu, double Ke[3][3]) {
    // Placeholder: Just for demonstration, not real stiffness calculation
    Ke[0][0] = E / (1 - nu * nu);
    Ke[0][1] = E * nu / (1 - nu * nu);
    Ke[0][2] = 0;

    Ke[1][0] = E * nu / (1 - nu * nu);
    Ke[1][1] = E / (1 - nu * nu);
    Ke[1][2] = 0;

    Ke[2][0] = 0;
    Ke[2][1] = 0;
    Ke[2][2] = E / (2 * (1 + nu)); // Shear modulus for the 2D case
}

// Function to assemble the global stiffness matrix with OpenMP parallelization
void StiffnessMatrixAssembler::assembleGlobalStiffness(const Graph& mesh, double E, double nu, std::vector<std::vector<double>>& K_global) {
    // Get the number of nodes in the mesh
    size_t totalNodes = mesh.points.size();

    // Initialize global stiffness matrix (size: totalNodes x totalNodes)
    K_global.resize(totalNodes, std::vector<double>(totalNodes, 0.0));

    // Loop over each triangle element in the mesh in parallel
#pragma omp parallel for
    for (int elemIdx = 0; elemIdx < mesh.triangles.size(); ++elemIdx) {
        const auto& element = mesh.triangles[elemIdx];

        // Local stiffness matrix for the current triangle element (3x3 matrix for 2D triangles)
        double Ke[3][3];
        computeElementStiffness(element, E, nu, Ke);

        // Indices of the triangle nodes in the global system
        int globalIndices[3] = {
            mesh.getNodeIndex(element.a),
            mesh.getNodeIndex(element.b),
            mesh.getNodeIndex(element.c)
        };

        // Assemble the local stiffness matrix into the global stiffness matrix
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int global_i = globalIndices[i]; // Get the global index of node 'i'
                int global_j = globalIndices[j]; // Get the global index of node 'j'

                // Critical section to prevent race conditions when writing to K_global
#pragma omp atomic
                K_global[global_i][global_j] += Ke[i][j];
            }
        }
    }
}

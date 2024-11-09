#include "DelaunayTriangulationOpenMP.h"
#include "StiffnessMatrixAssemblerOpenMP.h"
#include "BoundaryConditionsOpenMP.h"
#include "SolverOpenMP.h"
#include "PostProcessorOpenMP.h"
#include "chrono"

int main() {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    Graph graph;
    std::vector<Point> points;

    int count = 0;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            std::cout << "Adding Point :" << count << std::endl;
            //points.push_back(Point(i, j));
            graph.addPoint(Point(i, j));

            count++;
        }
    }

    int numDomains = 4;  // Example: divide the graph into 4 domains for parallelism
    //graph.addPointsWithDomainDecomposition(points, numDomains);

    //graph.addPoint(Point(4.5, 0.5));

    // Print triangles
    //graph.printTriangles();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto phase1_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> phase1_duration = phase1_end - start;
    std::cout << "Time taken for Phase 1 (Printing triangles): " << phase1_duration.count() << " seconds." << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StiffnessMatrixAssembler assembler;
    BoundaryConditions bc;
    Solver solver;
    PostProcessor postProcessor;


    // Define material properties
    double E = 210e9; // Young's modulus in Pascals
    double nu = 0.3;  // Poisson's ratio

    // Create a global stiffness matrix (initialize as an empty 2D vector)
    std::vector<std::vector<double>> K_global;
    std::vector<double> F(graph.points.size(), 0.0); // Initialize force vector with zeros

    // Assemble the global stiffness matrix
    assembler.assembleGlobalStiffness(graph, E, nu, K_global);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto phase2_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> phase2_duration = phase2_end - phase1_end;
    std::cout << "Time taken for Phase 2 (Assembling global stiffness matrix): " << phase2_duration.count() << " seconds." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Print the global stiffness matrix
    //std::cout << "\nGlobal Stiffness Matrix:" << std::endl;
    //for (const auto& row : K_global) {
    //    for (const auto& value : row) {
    //        std::cout << value << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    // Define fixed nodes (for example: fix node 0 and node 1)
    std::vector<int> fixedNodes = { 0, 1 };

    // Define loads (for example: apply a load of 1000 N at node 2)
    std::vector<Load> loads = {
        { 2, 1000.0 } // Load of 1000 N applied at node 2
    };

    // Apply boundary conditions
    bc.applyBoundaryConditions(K_global, F, fixedNodes, loads);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto phase3_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> phase3_duration = phase3_end - phase2_end;
    std::cout << "Time taken for Phase 3 (Applying boundary conditions): " << phase3_duration.count() << " seconds." << std::endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Print the modified global stiffness matrix and force vector
    //std::cout << "\nGlobal Stiffness Matrix after Boundary Conditions:" << std::endl;
    //for (const auto& row : K_global) {
    //    for (const auto& value : row) {
    //        std::cout << value << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    std::cout << "\nForce Vector:" << std::endl;
    for (const auto& value : F) {
        std::cout << value << "\t";
    }
    std::cout << std::endl;

    // Solve the system for displacements
    std::vector<double> displacements = solver.solveSystem(K_global, F);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto phase4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> phase4_duration = phase4_end - phase3_end;
    std::cout << "Time taken for Phase 4 (Solving system for displacements): " << phase4_duration.count() << " seconds." << std::endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Print the results
    //std::cout << "\nDisplacement Vector:" << std::endl;
    //for (const auto& value : displacements) {
    //    std::cout << value << "\t";
    //}
    //std::cout << std::endl;

    // Post-process to compute the stress field
    std::vector<std::vector<double>> stress_field = postProcessor.computeStressField(displacements, graph.triangles, graph.points, E, nu);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto phase5_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> phase5_duration = phase5_end - phase4_end;
    std::cout << "Time taken for Phase 5 (Post-processing to compute stress field): " << phase5_duration.count() << " seconds." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Print the stress results
    //std::cout << "\nStress Field:" << std::endl;
    //for (const auto& stress : stress_field) {
    //    std::cout << "Stress: ";
    //    for (const auto& s : stress) {
    //        std::cout << s << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    return 0;
}


// 100 Points : 7 Sec
// 400 Points : 35 Minutes
// CPU Usage  : 80-90%
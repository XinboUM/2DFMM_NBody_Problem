// This is the driver code that implements 2D MLFMM (Greengard, Rokhlin).
// The tree is implemented using 2D vector. The reason it is not implemented using a tree structure is 
// because I don't know how to traverse a certain level of a tree easily.
// Since the box number is represented using unsigned int, maximum # of total boxes allowed is ~4 billion, 
// i.e., this code can hold up to 65536 boxes in x and y directions, max number of level is 16.

#include <algorithm>
#include <iostream>
#include <iterator>
#include <ostream>
#include <time.h>

#include <math.h> // log
#include <vector>

#include "BoxTreeNode.hpp"
#include "FMMUtilities.hpp"
#include <cstdlib> // exit, strtod

int main(int argc, char **argv)
{

    //----------
    // Parse command line inputs: # of bodies and desired precision
    //----------
    if (argc < 3)
    {
        std::cerr << "You must pass the number of bodies and desired precision as command line parameters" << std::endl;
        std::exit(1);
    }
    else if (argc > 3)
    {
        std::cerr << "Too many command line parameters" << std::endl;
        std::exit(1);
    }

    const unsigned int N = std::stoul(argv[1]);       // total number of points

    const unsigned int nLevel = ceil(log(N)/log(4))+1; // total number of levels. Levels are indexed as 0 to nLevel-1.
    const unsigned int maxLevel = nLevel-1; // 0-based index of max level

    if (nLevel==1) // exit if FMM is not applicable. For the following we assume nLevel >= 2.
    {
        std::cout << "No FMM implementable. Exit." << std::endl;
        return 0;
    }
    const double precision = std::strtod(argv[2], NULL);
    const int nTerm = ceil(-log(precision)/log(2));
    
    std::cout << "==========\n";
    std::cout << "User input: N-body problem and FMM tree statistics:\n";
    std::cout << "----------\n";
    std::cout << "\t# of bodies = "<< N << "\n";
    std::cout << "\t# of box tree levels = "<< nLevel << "\n";
    std::cout << "\tprecision = "<< precision << "\n";
    std::cout << "\t# of terms in source and target series expansion = " << nTerm << std::endl;
    std::cout << "==========\n";

    std::cout << "\n";
    const BoxGeometry bodiesRegion; // default bodies region

    //----------
    // FOR DEBUGGING PURPOSE: N Fixed Bodies Creation
    // The locations are normalized (between 0 to 1), the associated function is 0 to 1000.
    //----------
    // std::vector<Body> bodies;
    // generateAndPrintFixedBodies(bodies, N, bodiesRegion);

    //----------
    // N Random Bodies Creation
    // The locations are normalized (between 0 to 1), the associated function is 0 to 1000
    //----------
    std::vector<Body> bodies;
    srand(clock());
    generateRandomBodies(bodies, N, bodiesRegion);
    
    //----------
    // FOR DEBUGGING PURPOSE: Output the statistics of the tree
    //----------
    std::cout << "\n";
    std::cout << "==========\n";
    std::cout << "FMM tree statistics:\n";
    std::cout << "----------\n";

    std::vector<BoxTreeFixedLevelStatistics*> boxTreeStat(nLevel);
    for (unsigned char level = 0; level < nLevel; level++)
    {
        boxTreeStat[level] = new BoxTreeFixedLevelStatistics(bodiesRegion, level);
        std::cout << "level l=" << int(level) << "\n";
        boxTreeStat[level]->printBoxTreeFixedLevelStatistics();
    }

    std::cout << "==========\n";
    
    //----------
    // Build the nLevel-level Tree
    //----------
    std::cout << "\n";

    std::cout << "==========\n";
    std::cout << "FMM tree generation...\n";
    std::cout << "==========\n";

    std::vector<std::vector<BoxTreeNode*>> boxTree(nLevel);
    constructBoxTree(boxTree, bodiesRegion, boxTreeStat, nTerm, maxLevel);

    // std::cout << "\n";
    // printBoxTree(boxTree);
    
    //----------
    // Add bodies to the leaf level boxes
    //----------
    assignBodyToLeafBox(bodies, boxTree, bodiesRegion, boxTreeStat, maxLevel);
    std::cout << "\n";
    printBoxTree(boxTree);
    
    //----------
    // Upward Pass
    //----------
    // compute multipole expansion coefficients in the leaf level boxes
    for (unsigned int iLeafBox = 0; iLeafBox < boxTree[maxLevel].size(); iLeafBox++)
        boxTree[maxLevel][iLeafBox]->computeMultipoleExpansionCoefficients_leafBox();
    
    // compute multipole expansion coefficients in higher level boxes from its children
    for (short level = maxLevel-1; level >= 0; level--) // level has type short to contain its possible max value being max value of unsigned char
        for (unsigned int iBox = 0; iBox < boxTree[level].size(); iBox++)
            boxTree[level][iBox]->computeMultipoleExpansionCoefficientsFromChildren(); 
        
    //----------
    // Downward Pass
    //----------
    // From top level down, for each box, compute power series coefficients in the potential expansion due to source boxes in its interaction list
    for (unsigned int level = 2; level <= maxLevel; level++)
        for (unsigned int iBox = 0; iBox < boxTree[level].size(); iBox++)
            boxTree[level][iBox]->computePowerSeriesCoefficientsDueToInteractionList();

    // children box update its power series coefficients from its parent box
    for (unsigned int level = 3; level <= maxLevel; level++)
        for (unsigned int iBox = 0; iBox < boxTree[level].size(); iBox++)
            boxTree[level][iBox]->updatePowerSeriesCoefficientsFromParent();
        
    // Evaluate the power series expansion of the potential for each body in each leaf box
    for (unsigned int iBox = 0; iBox < boxTree[maxLevel].size(); iBox++)
        boxTree[maxLevel][iBox]->computePotentialFMM();

    //==========
    // Delete the tree and statistics tree
    //==========
    // delete the statistics tree
    for (BoxTreeFixedLevelStatistics* boxTreeFixedLevel : boxTreeStat)
        delete boxTreeFixedLevel; 
    boxTreeStat.clear();

    // delete the box tree from top down
    delete boxTree[0][0];

    // // Method 2. delete each box in the box tree. Using this method requires deleting the destructor in BoxTreeNode class.
    // // This is not as natural as deleting the box tree from top down because BoxTreeNode class should have a destructor, since
    // // generateChildren calls the new operator. 
    // for (unsigned int level = 0; level < boxTree.size(); level++)
    // {
    //     for (BoxTreeNode* box : boxTree[level])
    //         delete box;
    //     boxTree[level].clear();
    // }
    // boxTree.clear();
    //==========
    // Validate the potential results
    //==========
    std::cout << "==========\n";

    // Directly compute the potentials and print the results to console
    unsigned int nPotentialMax_NA = 10; // maximum amount of direct computed potential
    computeNBodyPotential_NA(bodies, nPotentialMax_NA);

    // print potential for each body
    for (unsigned int n = 0; n < N; n++)
    {
        std::cout << "Body ID = " << n << ": z=" << bodies[n].z << ", f=" << bodies[n].f << "\n";
        std::cout << "\tpotential: naive computed vs FMM = "
                  << bodies[n].potential_NA << " vs " << bodies[n].potential_NA << "\n";
    }

    std::cout << "==========\n";

    return 0;
}































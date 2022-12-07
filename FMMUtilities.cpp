#include "./FMMUtilities.hpp"
#include "BoxTreeNode.hpp"
#include <algorithm>
#include <complex>
#include "iostream"

void generateAndPrintFixedBodies(std::vector<Body>& bodies, unsigned int N, const BoxGeometry& bodiesRegion)
{
    bodies.resize(N);

    std::cout << "==========\n";
    std::cout << "The given bodies region is x range = (" 
              << bodiesRegion.minCoord[0] << ", " << bodiesRegion.maxCoord[0]
              << "); y range = ("
              << bodiesRegion.minCoord[1] << ", " << bodiesRegion.maxCoord[1] << ").\n";
    std::cout << "In this region, fixed bodies are generated/initialized as:\n";
    std::cout << "----------\n";

    double bodiesRegion_x = bodiesRegion.maxCoord[0] - bodiesRegion.minCoord[0];
    double bodiesRegion_y = bodiesRegion.maxCoord[1] - bodiesRegion.minCoord[1];

    for (unsigned int iBody = 0; iBody < N; iBody++)
    {
        bodies[iBody].r[0] = bodiesRegion.minCoord[0] + double(iBody+1) / double(N+1) * bodiesRegion_x;
        bodies[iBody].r[1] = bodiesRegion.minCoord[1] + double(iBody+1) / double(N+1) * bodiesRegion_y;
        bodies[iBody].z = bodies[iBody].r[0] + imagUnit* bodies[iBody].r[1];

        bodies[iBody].f = (iBody+1)*10;

        bodies[iBody].potential_FMM = complexOrigin;
        bodies[iBody].potential_NA = complexOrigin;

        std::cout << "Body ID = " << iBody << ": z=" << bodies[iBody].z << ", f=" << bodies[iBody].f 
                  << ", initilized FMM potential=" << bodies[iBody].potential_FMM 
                  << ", initilized NA potential=" << bodies[iBody].potential_NA << std::endl;
    }
    std::cout << "==========\n";
}

void generateRandomBodies(std::vector<Body>& bodies, unsigned int N, const BoxGeometry& bodiesRegion)
{
    bodies.resize(N);

    std::cout << "==========\n";
    std::cout << "Random bodies generated/initialized:\n";
    std::cout << "----------\n";

    double bodiesRegion_x = bodiesRegion.maxCoord[0] - bodiesRegion.minCoord[0];
    double bodiesRegion_y = bodiesRegion.maxCoord[1] - bodiesRegion.minCoord[1];

    for (unsigned int n = 0; n < N; n++)
    {
        bodies[n].r[0] = bodiesRegion.minCoord[0] + (double)rand() / (double)RAND_MAX * bodiesRegion_x;
        bodies[n].r[1] = bodiesRegion.minCoord[1] + (double)rand() / (double)RAND_MAX * bodiesRegion_y;
        bodies[n].z = bodies[n].r[0] + imagUnit* bodies[n].r[1];

        bodies[n].f = (double)rand() / (double)RAND_MAX * 1000; // f being 0 not handled (even though not statistically possible)

        bodies[n].potential_FMM = complexOrigin;
        bodies[n].potential_NA = complexOrigin;
        
        std::cout << "Body ID = " << n << ": z=" << bodies[n].z << ", f=" << bodies[n].f 
                  << ", initilized FMM potential=" << bodies[n].potential_FMM
                  << ", initilized NA potential=" << bodies[n].potential_NA << std::endl;
    }
    
}

// find nearest neighbours and interaction list for each box in tree
void findNearestNeighboursAndInteractionListForEachBoxInTree(std::vector<std::vector<BoxTreeNode*>>& boxTree, 
                                   const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStat)
{
    // for boxes on level 1
    for (unsigned char iBox_l1 = 0; iBox_l1 < 4; iBox_l1++)
    {
        for (unsigned char jBox_l1 = 0; jBox_l1 < 4; jBox_l1++)
        {
            if (jBox_l1 != iBox_l1)
                boxTree[1][iBox_l1]->addBoxToNearestNeighbourList(boxTree[1][jBox_l1]);
        }
    }
    // for boxes on higher levels
    // fixed level 
    for (unsigned char level = 2; level < boxTree.size(); level++)
    {
        double largestNearestNeighbourCenterDistance = 
                    (boxTreeStat[level]->edgeLength[0]*boxTreeStat[level]->edgeLength[0]+
                     boxTreeStat[level]->edgeLength[1]*boxTreeStat[level]->edgeLength[1]);
        
        unsigned char parentLevel = level-1;
        // for each box on this level, find boxes in its nearest neighbour list and interaction list
        for (unsigned int iBox = 0; iBox < boxTree[level].size(); iBox++)
        {   
            const BoxTreeNode* parent = boxTree[level][iBox]->getParent();
            // parent's children excluding self box are all in the nearest neighbour list
            std::array<const BoxTreeNode*, 4> children_parent = parent->getChildren();
            for (const BoxTreeNode* child_parent : children_parent)
                if (child_parent != boxTree[level][iBox])
                    boxTree[level][iBox]->addBoxToNearestNeighbourList(child_parent);

            // nearest neighbours of parent
            std::vector<const BoxTreeNode*> nearestNeighbours_parent = parent->getNearestNeighbours();
            
            // loop over children of parent's nearest neighbours, and add boxes to the nearest neighbour or interaction list of this box
            for (const BoxTreeNode* nearestNeighbour_parent : nearestNeighbours_parent)
            {
                std::array<const BoxTreeNode*, 4> children_parentNearestNeighbour = nearestNeighbour_parent->getChildren();

                // compute the squared center distance between the box and the candidate box to decide if the 
                // candidate is target box's nearest neighbour or in its interaction list
                for (const BoxTreeNode* candidate : children_parentNearestNeighbour)
                {
                    double centerDistSq = findDistanceSquared(candidate->getBoxCenter(),boxTree[level][iBox]->getBoxCenter());
                    
                    if (centerDistSq <= largestNearestNeighbourCenterDistance)
                        boxTree[level][iBox]->addBoxToNearestNeighbourList(candidate);
                    else
                        boxTree[level][iBox]->addBoxToInteractionList(candidate);
                }
            }
        }
    }
}

void printBoxTree(const std::vector<std::vector<BoxTreeNode*>>& boxTree)
{
    std::cout << "==========\n";
    std::cout << "Statistics of each box in the FMM tree:\n";
    std::cout << "----------\n";
    std::cout << "The number of rows/levels in the box tree is " << boxTree.size() << ".\n";
    std::cout << "----------\n";

    // output statistics for a fixed level
    for (unsigned char iRow = 0; iRow < boxTree.size(); iRow++)
    {
        std::cout << "==========\n";
        std::cout << "Level " << int(iRow) << ":\n";
        std::cout << "# of columns/boxs = " << boxTree[iRow].size() << "\n";

        // output statistics for each box at a given level
        for (unsigned int iCol = 0; iCol < boxTree[iRow].size(); iCol++)
        {
            std::cout << "----------\n";
            std::cout << "Box ID = " << iCol << "\n";
            boxTree[iRow][iCol]->printBox(boxTree.size()-1);
        }

    }
    std::cout << "==========\n";
}

unsigned int findLeafBoxIDForGivenBody(const Body& body, const BoxGeometry& bodiesRegion, 
                                    const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats, unsigned char maxLevel)
{   
    unsigned short xID = (body.r[0]-bodiesRegion.minCoord[0])/(bodiesRegion.maxCoord[0]-bodiesRegion.minCoord[0])
                        *boxTreeStats[maxLevel]->nBox_1D;
    unsigned short yID = (body.r[1]-bodiesRegion.minCoord[1])/(bodiesRegion.maxCoord[1]-bodiesRegion.minCoord[1])
                        *boxTreeStats[maxLevel]->nBox_1D;

    return (xID<<maxLevel)+yID;
}

unsigned int findLeafBoxIDForGivenPoint(const std::complex<double>& point, const BoxGeometry& bodiesRegion, 
                                    const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats, unsigned char maxLevel)
{   
    unsigned short xID = (point.real()-bodiesRegion.minCoord[0])/(bodiesRegion.maxCoord[0]-bodiesRegion.minCoord[0])
                        *boxTreeStats[maxLevel]->nBox_1D;
    unsigned short yID = (point.imag()-bodiesRegion.minCoord[1])/(bodiesRegion.maxCoord[1]-bodiesRegion.minCoord[1])
                        *boxTreeStats[maxLevel]->nBox_1D;

    return (xID<<maxLevel)+yID;
}

void assignBodyToLeafBox(std::vector<Body>& bodies, std::vector<std::vector<BoxTreeNode*>>& boxTree, 
                         const BoxGeometry& bodiesRegion, const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats,
                         unsigned char maxLevel)
{
    for (unsigned int iBody = 0; iBody < bodies.size(); iBody++)
    {
        unsigned int iBox = findLeafBoxIDForGivenBody(bodies[iBody], bodiesRegion, boxTreeStats, maxLevel);
        // std::cout << "The box ID for body " << iBody << " is "<< iBox << "\n";
        boxTree[maxLevel][iBox]->addBody(&(bodies[iBody]));
    }

};

double findDistanceSquared(const std::complex<double>& z1, const std::complex<double>& z2)
{
    return (z1.real()-z2.real())*(z1.real()-z2.real()) + (z1.imag()-z2.imag())*(z1.imag()-z2.imag());
}

unsigned int binomialCoeff(unsigned int n, unsigned int k)
{
    unsigned int C[k + 1];
    memset(C, 0, sizeof(C));
 
    C[0] = 1; // nC0 is 1
 
    for (unsigned int i = 1; i <= n; i++) {
        // Compute next row of pascal triangle using the previous row
        for (int j = (i<k)?i:k; j > 0; j--)
            C[j] = C[j] + C[j - 1];
    }
    return C[k];
}

// compute the potential directly as the naive approach
void computeNBodyPotential_NA(std::vector<Body>& bodies, unsigned int nPotentialMax_NA)
{
    for (unsigned int iTar = 0; iTar < bodies.size(); iTar++)
    {
        for (unsigned int iSrc = 0; iSrc < bodies.size(); iSrc++)
        {
            if (iSrc == iTar)
                continue;
            else
                bodies[iTar].potential_NA += std::log(bodies[iTar].z-bodies[iSrc].z) * bodies[iSrc].f;
        }
        
        if (iTar >= nPotentialMax_NA)
            break;
    }
}

template< class T >
void reorder(std::vector<T>& v, std::vector<unsigned int> const &order )  
{   
    for ( int s = 1, d; s < order.size(); ++ s ) 
    {
        for ( d = order[s]; d < s; d = order[d] ) ;
        if ( d == s ) while ( d = order[d], d != s ) std::swap( v[s], v[d] );
    }
}

// O(N)
void rearrangeLeafBoxes(std::vector<BoxTreeNode*>& leafBoxes, const BoxGeometry& bodiesRegion, 
                        const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats, unsigned char maxLevel)
{
    // find the sequential ID for each leaf level box
    std::vector<unsigned int> boxID_sq(leafBoxes.size());
    for (unsigned int iLeafBox = 0; iLeafBox < leafBoxes.size(); iLeafBox++)
        boxID_sq[iLeafBox] = findLeafBoxIDForGivenPoint(leafBoxes[iLeafBox]->getBoxCenter(), bodiesRegion, boxTreeStats, maxLevel);
    
    // reorder the leaf boxes according to the sequential ID
    reorder(leafBoxes, boxID_sq);
};


void constructBoxTree(std::vector<std::vector<BoxTreeNode*>>& boxTree, const BoxGeometry& bodiesRegion, 
                          const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStat, unsigned int nTerm, unsigned char maxLevel)
{
    // generate the root
    boxTree[0].resize(1);
    boxTree[0][0] = new BoxTreeNode(bodiesRegion, nTerm);
    
    // for each box in the current level, construct the children of that box and fill the children to the next level
    for (int parentLevel = 0; parentLevel < maxLevel; parentLevel++)
    {
        int childrenLevel = parentLevel+1;
        for (int iParent = 0; iParent < boxTree[parentLevel].size(); iParent++)
        {
            // generate children for iParent box on parentLevel
            std::array<BoxTreeNode*, 4> children;
            boxTree[parentLevel][iParent]->generateChildren(nTerm, children);
        
            // assign children for iParent to the tree
            boxTree[childrenLevel].insert(boxTree[childrenLevel].begin()+boxTree[childrenLevel].size(), 
                                          std::begin(children), std::end(children));
        }
    }
    rearrangeLeafBoxes(boxTree[maxLevel], bodiesRegion, boxTreeStat, maxLevel);

    // Now the boxes in tree are constructed, we then find the nearest neighbours and interaction list for each box.
    findNearestNeighboursAndInteractionListForEachBoxInTree(boxTree, boxTreeStat);
}
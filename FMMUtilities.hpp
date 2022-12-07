#ifndef _FMMUtilities_HPP_
#define _FMMUtilities_HPP_

// #include <complex>
// #include <vector>
// #include <iostream>
#include "BoxTreeNode.hpp"
#include <vector>

void generateAndPrintFixedBodies(std::vector<Body>& bodies, unsigned int N, const BoxGeometry& bodiesRegion);
void generateRandomBodies(std::vector<Body>& bodies, unsigned int N, const BoxGeometry& bodiesRegion);

void findNearestNeighboursAndInteractionListForEachBoxInTree(std::vector<std::vector<BoxTreeNode*>>& boxTree, 
                                        const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStat);

void printBoxTree(const std::vector<std::vector<BoxTreeNode*>>& boxTree);

unsigned int findLeafBoxIDForGivenBody(const Body& body, const BoxGeometry& bodiesRegion, 
                                    const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats, unsigned char maxLevel);
unsigned int findLeafBoxIDForGivenPoint(const std::complex<double>& center, const BoxGeometry& bodiesRegion, 
                                    const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats, unsigned char maxLevel);

void assignBodyToLeafBox(std::vector<Body>& bodies, std::vector<std::vector<BoxTreeNode*>>& boxTree, 
                         const BoxGeometry& bodiesRegion, const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats,
                         unsigned char maxLevel);

double findDistanceSquared(const std::complex<double>& z1, const std::complex<double>& z2);

void computeNBodyPotential_NA(std::vector<Body>& bodies, unsigned int nPotentialMax_NA);

// swap two boxes. Now I am having trouble putting it into other files
template< class T >
void reorder(std::vector<T>& v, std::vector<unsigned int> const &order);

// O(N)
void rearrangeLeafBoxes(std::vector<BoxTreeNode*>& leafBoxes, const BoxGeometry& bodiesRegion, 
                        const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStats, unsigned char maxLevel);

void constructBoxTree(std::vector<std::vector<BoxTreeNode*>>& boxTree, const BoxGeometry& bodiesRegion, 
                          const std::vector<BoxTreeFixedLevelStatistics*>& boxTreeStat, unsigned int nTerm, unsigned char maxLevel);

#endif
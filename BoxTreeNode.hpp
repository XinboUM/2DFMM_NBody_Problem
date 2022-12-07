#ifndef _TREENODE_H_
#define _TREENODE_H_

#include <algorithm>
#include <array>
#include <complex>
#include <iostream>
#include <vector>
#include <numeric>

const std::complex<double> complexOrigin(0.0,0.0), imagUnit(0.0,1.0);

unsigned int binomialCoeff(unsigned int n, unsigned int k);

struct Body
{
    double r[2];// it might seem redundant to keep both real and complex locations, but we keep them for now.
    double f; // f can represent charge or mass

    std::complex<double> z = complexOrigin;
    std::complex<double> potential_NA = complexOrigin;
    std::complex<double> potential_FMM = complexOrigin;
};

struct BoxGeometry
{
    // The hard-coded 2 is because we are solving 2D problem.
    // 0th entry is x-related; 1st entry is y-related.
    std::array<double, 2> minCoord;
    std::array<double, 2> maxCoord;

    // default constructor that sets x range = y range = (0,1)
    BoxGeometry() 
    {
        minCoord[0] = 0; minCoord[1] = 0; 
        maxCoord[0] = 1; maxCoord[1] = 1; 
    }

    // constructor that allows the user to 
    BoxGeometry(double xmin, double xmax, double ymin, double ymax) 
    {
        minCoord[0] = xmin; minCoord[1] = ymin; 
        maxCoord[0] = xmax; maxCoord[1] = ymax; 
    }
};

struct BoxTreeFixedLevelStatistics
{
public:
    std::array<double, 2> edgeLength; // the edge length of each box
    unsigned int nBox_1D; // # of boxes in the x and y directions
    unsigned long nBox; // total # of boxes at this level
    
    BoxTreeFixedLevelStatistics(BoxGeometry box, unsigned char level)
    {
        nBox_1D = pow(2,level);    nBox = pow(4,level);
        edgeLength[0] = (box.maxCoord[0]-box.minCoord[0])/nBox_1D;    
        edgeLength[1] = (box.maxCoord[1]-box.minCoord[0])/nBox_1D; 
    };
    void printBoxTreeFixedLevelStatistics()
    {
        std::cout << "\ttotal # of boxes at this level = " << nBox << "\n";
        std::cout << "\t# of boxes in x and y direction = " << nBox_1D << "\n";
        std::cout << "\teach box has edge length in x direction = " << edgeLength[0] 
                  << ", edge length in y direction = " << edgeLength[1] << "\n";
    };
};

class BoxTreeNode
{
public: 
    BoxTreeNode(const BoxGeometry& boxGeom, unsigned char nTerm);
    BoxTreeNode(double xmin, double xmax, double ymin, double ymax, unsigned char nTerm, const BoxTreeNode* parent);
    ~BoxTreeNode()
    {
        for (const BoxTreeNode* child : children_)
        {
            if (child != nullptr)
                delete child;
            child = nullptr;
            bodies_.clear();
        }
    };
    void addBody(Body* body) {bodies_.push_back(body);};
    void addBoxToNearestNeighbourList(const BoxTreeNode* box) {nearestNeighbourList_.push_back(box);};
    void addBoxToInteractionList(const BoxTreeNode* box) {interactionList_.push_back(box);};
    void generateChildren(unsigned char nTerm, std::array<BoxTreeNode*, 4>& chilren); // children is an output
    
    //----------
    // const member functions
    //----------
    void printBox(unsigned int maxLevel) const;
    // getters
    const BoxGeometry& getBoxRegion() const {return boxGeom_;};
    const std::complex<double>& getBoxCenter() const {return center_;};
    const std::vector<std::complex<double>>& getMultipoleCoeff() const {return multipole_coeff_; };
    const std::vector<std::complex<double>>& getPowerSeriesCoeff() const {return powerSeries_coeff_; };
    const std::vector<Body*>& getBodies() const {return bodies_;};
    const std::array<const BoxTreeNode*, 4>& getChildren() const {return children_;};
    const BoxTreeNode* getParent() const {return parent_;};
    const std::vector<const BoxTreeNode*>& getNearestNeighbours() const {return nearestNeighbourList_;};
    
    // leaf box function
    void computeMultipoleExpansionCoefficients_leafBox(); 
    void computePotentialFMM();

    // non-leaf box function
    void computeMultipoleExpansionCoefficientsFromChildren();
    
    // nTerm is the number of power series coefficients. -1^k can be optimized.
    void computePowerSeriesCoefficientsDueToInteractionList();
    void updatePowerSeriesCoefficientsFromParent();
    
private:
    // geometrical properties
    BoxGeometry boxGeom_;
    std::complex<double> center_;

    // Pointer to other related boxes in the tree. We add const to all of them because a box does not change other box's property.
    const BoxTreeNode* parent_;    // parent ID in parent's level
    std::array<const BoxTreeNode*, 4> children_;    // children ID in children's level
    std::vector<const BoxTreeNode*> nearestNeighbourList_; // can not initialize with size in a class?    
    std::vector<const BoxTreeNode*> interactionList_;  

    // multipole coefficients (source representation) and power series coefficients (target representation)
    std::vector<std::complex<double>> multipole_coeff_;
    std::vector<std::complex<double>> powerSeries_coeff_;
    
    // Pointer to the bodies for leaf level boxes.
    std::vector<Body*> bodies_;
};

#endif
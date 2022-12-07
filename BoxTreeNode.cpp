
#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <vector>
#include "BoxTreeNode.hpp"
#include "FMMUtilities.hpp"

// used for root box initialization
BoxTreeNode::BoxTreeNode(const BoxGeometry& boxGeom, unsigned char nTerm)
{
    // set up geometrical location
    boxGeom_ = boxGeom;
    center_ = std::complex<double>(0.5*(boxGeom_.minCoord[0]+boxGeom_.maxCoord[0]), 0.5*(boxGeom_.minCoord[1]+boxGeom_.maxCoord[1]));

    // Initialize the related boxes pointers
    parent_ = nullptr;
    std::fill(children_.begin(), children_.end(), nullptr);
    nearestNeighbourList_.resize(0);
    interactionList_.resize(0);

    // Initialize coefficients list
    multipole_coeff_.resize(nTerm);
    std::fill(multipole_coeff_.begin(), multipole_coeff_.end(), 0);
    powerSeries_coeff_.resize(nTerm);
    std::fill(powerSeries_coeff_.begin(), powerSeries_coeff_.end(), 0);

    // Initialize bodies list
    bodies_.resize(0);
}

// used for children generation
BoxTreeNode::BoxTreeNode(double xmin, double xmax, double ymin, double ymax, unsigned char nTerm, const BoxTreeNode* parent)
{
    // set up geometrical location
    boxGeom_.minCoord[0] = xmin;    boxGeom_.maxCoord[0] = xmax;
    boxGeom_.minCoord[1] = ymin;    boxGeom_.maxCoord[1] = ymax;
    center_ = std::complex<double>(0.5*(boxGeom_.minCoord[0]+boxGeom_.maxCoord[0]), 0.5*(boxGeom_.minCoord[1]+boxGeom_.maxCoord[1]));

    // Initialize the related boxes pointers
    parent_ = parent;
    std::fill(children_.begin(), children_.end(), nullptr);
    nearestNeighbourList_.resize(0);
    interactionList_.resize(0);

    // Initialize coefficients list
    multipole_coeff_.resize(nTerm);
    powerSeries_coeff_.resize(nTerm);

    // Initialize bodies list
    bodies_.resize(0);
}

void BoxTreeNode::printBox(unsigned int maxLevel) const
{
    // ID
    // std::cout << "\tparent ID = " << parentID_ << "\n";
    // // std::cout << "\tparent xID = " << parent_xID_ << "\n";
    // // std::cout << "\tparent yID = " << parent_yID_ << "\n";
    // std::cout << "\tchildren ID = " << childrenID_[0] << ", " << childrenID_[1] << ", " << childrenID_[2] << ", " << childrenID_[3] << "\n";
    
    // nearest neighbours
    if (nearestNeighbourList_.empty()) 
        std::cout << "\tNo nearest neighbours.\n";
    else    
        std::cout << "\t# of nearest neighbours = " << nearestNeighbourList_.size() << "\n";

    // interaction list
    if (interactionList_.empty()) 
        std::cout << "\tNo interaction List.\n";
    else    
        std::cout << "\t# of interaction List = " << interactionList_.size() << "\n";
    
    // geometrical range
    std::cout << "\tx range = (" << boxGeom_.minCoord[0] << ", " << boxGeom_.maxCoord[0] << ")\n";
    std::cout << "\ty range = (" << boxGeom_.minCoord[1] << ", " << boxGeom_.maxCoord[1] << ")\n";

    // center 
    std::cout << "\tcenter = " << center_.real() << "+" << center_.imag() << "i\n";

    // // bodies contained
    // if (level_ == maxLevel)
    // {
    //     if (bodies_.empty())
    //         std::cout << "No body is contained in this box.\n";
    //     else
    //     {
    //         std::cout << "Contains bodies: \n";
    //         for (const Body* body : bodies_)
    //             std::cout << "\tz=" << body->z << ", f=" << body->f << ", FMM potential=" << body->potential_FMM << "\n";
    //     }
    // }
    std::cout << "----------\n";

}

void BoxTreeNode::generateChildren(unsigned char nTerm, std::array<BoxTreeNode*, 4>& children)
{
    children[0] = new BoxTreeNode(boxGeom_.minCoord[0], center_.real(), boxGeom_.minCoord[1], center_.imag(), nTerm, this);
    children[1] = new BoxTreeNode(boxGeom_.minCoord[0], center_.real(), center_.imag(), boxGeom_.maxCoord[1], nTerm, this);
    children[2] = new BoxTreeNode(center_.real(), boxGeom_.maxCoord[0], boxGeom_.minCoord[1], center_.imag(), nTerm, this);
    children[3] = new BoxTreeNode(center_.real(), boxGeom_.maxCoord[0], center_.imag(), boxGeom_.maxCoord[1], nTerm, this);
    
    std::copy(std::begin(children), std::end(children), std::begin(children_));
}

void BoxTreeNode::computeMultipoleExpansionCoefficients_leafBox()
{
    for (unsigned int iBody = 0; iBody < bodies_.size(); iBody++)
    {
        multipole_coeff_[0] += bodies_[iBody]->f;
        for (unsigned int iTerm = 1; iTerm < multipole_coeff_.size(); iTerm++) // if nTerm == 0 or 1, this for loop wouldn't be executed
            multipole_coeff_[iTerm] += - bodies_[iBody]->f * pow(bodies_[iBody]->z-center_, iTerm) / double(iTerm);
    }
};

void BoxTreeNode::computeMultipoleExpansionCoefficientsFromChildren()
{
    for (const BoxTreeNode* child : children_)
    {
        const std::vector<std::complex<double>>& multipoleCoeff_child = child->getMultipoleCoeff(); // child's multipole expansion coefficients
        std::complex<double> z0 = child->getBoxCenter() - center_;// the coordinate of child box's center when complex plane is originated at my center
        // 0th multipole coefficient
        multipole_coeff_[0] = multipoleCoeff_child[0];
        // 1 to (nTerm-1)th entries of multipole coefficients
        for (unsigned int l = 1; l < multipole_coeff_.size(); l++)
        {
            for (unsigned int k = 1; k <= l; k++)
                multipole_coeff_[l] += multipoleCoeff_child[k] * pow(z0, l-k) * double(binomialCoeff(l-1, k-1));
            
            multipole_coeff_[l] -= multipoleCoeff_child[0] * pow(z0, l) / double(l);
        }
    }
    
};

void BoxTreeNode::computePowerSeriesCoefficientsDueToInteractionList() 
{
    // loop over all boxes in the interaction list of current box
    for (const BoxTreeNode* srcBoxIL : interactionList_)
    {
        const std::vector<std::complex<double>>& mutipoleCoeff_srcBoxIL = srcBoxIL->getMultipoleCoeff();
        std::complex<double> z0 = srcBoxIL->getBoxCenter() - center_;

        // 0th entry of the power series coefficients
        powerSeries_coeff_[0] += mutipoleCoeff_srcBoxIL[0] * std::log(-z0);
        for (unsigned int k = 1; k < multipole_coeff_.size(); k++)
            powerSeries_coeff_[0] += mutipoleCoeff_srcBoxIL[k] * pow(-1.0/z0, k);

        // 1 to (nTerm-1)th entries of the power series coefficients
        for (unsigned int l = 1; l < powerSeries_coeff_.size(); l++)
        {
            powerSeries_coeff_[l] -= mutipoleCoeff_srcBoxIL[0] / double(l) / pow(z0, l);
            
            for (unsigned int k = 1; k < multipole_coeff_.size(); k++)
                powerSeries_coeff_[l] += mutipoleCoeff_srcBoxIL[k] * double(binomialCoeff(l+k-1, k-1)) * pow(-1, k) / pow(z0, k+l);
        }
    }
};

void BoxTreeNode::updatePowerSeriesCoefficientsFromParent()
{
    const std::vector<std::complex<double>>& powerSeriesCoeff_parent = parent_->getPowerSeriesCoeff();
    std::complex<double> z0 = parent_->getBoxCenter() - center_;

    for (unsigned int l = 0; l < powerSeries_coeff_.size(); l++)
        for (unsigned int k = l; k < powerSeries_coeff_.size(); k++)
            powerSeries_coeff_[l] += powerSeriesCoeff_parent[l] * double(binomialCoeff(k, l)) * pow(-z0, k-l);

};

// compute potential for each body in the leaf box
void BoxTreeNode::computePotentialFMM() 
{
    for (unsigned int iBody = 0; iBody < bodies_.size(); iBody++)
    {
        // compute far-field (FF) potential
        std::complex<double> z = bodies_[iBody]->z - center_; // location of the point in the current box
        for (unsigned int l = 0; l < powerSeries_coeff_.size(); l++)
            bodies_[iBody]->potential_FMM += powerSeries_coeff_[l] * pow(z, l);
    
        // compute near-field (NF) potential
        // 1. NF potential due to sources in the same box
        for (unsigned int jBody = 0; jBody < bodies_.size(); jBody++)
        {
            if (jBody == iBody)
                continue;
            else
                bodies_[iBody]->potential_FMM += std::log(bodies_[iBody]->z-bodies_[jBody]->z) * bodies_[jBody]->f;
        }
        
        // 2. NF potential due to sources in the nearest neighbours
        for (const BoxTreeNode* nearestNeighbour : nearestNeighbourList_)
        {
            const std::vector<Body*>& bodyList = nearestNeighbour->getBodies();

            for (const Body* src : bodyList)
                bodies_[iBody]->potential_FMM += std::log(bodies_[iBody]->z-src->z) * src->f;
        }
    }
};
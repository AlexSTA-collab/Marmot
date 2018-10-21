#pragma once
#include "bftTypedefs.h"
#include "bftFiniteElement.h"
#include <map>

template< int nDim, int  nNodes>
class BftGeometryElement {
    /* This is the Geometry Base element, which serves as a bases for all BftUels.
     * It corresponds to the GeometryElement in mpFEM,
     * although this as a static templated version.
     *
     * BftUels (corresponding do DofElements in mpFEM) can inherit from this element,
     * and access shape functions, derivatives and B Operator
     *
     * The element automatically determines its shape by the given  nDimension and number of nodes
     * */

public:

    /*Typedefs*/
    static constexpr int    VoigtSize = ( ( ( nDim*nDim ) +  nDim ) / 2 );

    typedef Matrix<double,  nDim, 1 >                                   XiSized;
    typedef Matrix<double,  nDim*nNodes, 1 >                            CoordinateVector;
    typedef Matrix<double,  nDim,  nDim>                                JacobianSized;

    typedef Matrix<double, 1, nNodes>                                   NSized;
    typedef Matrix<double, nDim, nNodes*nDim>                           NBSized;
    typedef Matrix<double, nDim, nNodes>                                dNdXiSized;
    typedef Matrix<double, VoigtSize, nNodes*nDim >                     BSized;


    /*Properties*/
    Map<const CoordinateVector>                                         coordinates;
    const bft::FiniteElement::ElementShapes                             shape;

    /*Methods*/
    BftGeometryElement() :
        coordinates ( nullptr ),
        shape ( bft::FiniteElement::getElementShapeByMetric ( nDim, nNodes ) )
    {};

    std::string getElementShape()
    {
        static std::map< bft::FiniteElement::ElementShapes, std::string> shapes = {
            {bft::FiniteElement::Truss2,    "truss2"},
            {bft::FiniteElement::Quad4,     "quad4"},
            {bft::FiniteElement::Quad8,     "quad8"},
            {bft::FiniteElement::Hexa8,     "hexa8"},
            {bft::FiniteElement::Hexa20,    "hexa20"}
        };

        return shapes[this->shape];
    }

    void initializeYourself ( const double* coords )
    {
        new ( &coordinates ) Map<const CoordinateVector> ( coords );
    }

    /*Please specialize these functions for each element individially
     *.cpp file.
     *Fully specialized templates are precompiled in bftMechanics (rather than the unspecialized and partially specialized templates)
     * */
    NSized   N ( const XiSized&   xi );
    dNdXiSized dNdXi ( const XiSized&   xi );
    BSized   B ( const dNdXiSized&  dNdX );

    /*These functions are equal for each element and independent of node number and  nDimension*/
    NBSized NB ( const NSized&  N )
    {
        return bft::FiniteElement::NB<nDim, nNodes> ( N );
    }
    JacobianSized Jacobian ( const dNdXiSized& dNdXi )
    {
        return bft::FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
    }
    dNdXiSized dNdX ( const dNdXiSized& dNdXi, const JacobianSized& JacobianInverse )
    {
        return ( dNdXi.transpose() * JacobianInverse ).transpose();
    }

};

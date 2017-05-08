#pragma once 
#include <string>
#include "bftUel.h"

namespace bft{
    /*copy of bftTypedefs.h*/
        typedef void (*pUmatType)( double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&,const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);

        typedef void (*pSimpleUelWithUmatType)( double [], double [], double [], const int &, const double [], const int &, const double [], const double [], const double [], const double [2], const double &, const int& , double &, const int [], const int &, bft::pUmatType, int );

}

namespace userLibrary{

    bft::pUmatType getUmatById(int id);
    bft::pUmatType getUmatByName(const std::string& nameUpperCase);

    BftUel* UelFactory(int id, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    void umatPlaneStressWrapped(const bft::pUmatType, double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&, const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);

    static const int sizeGeostaticDefinition = 5;
}

namespace Abaqus{
    
    enum UelFlags1{
        Geostatic=61,
    };
}



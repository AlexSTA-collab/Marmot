#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotExceptions.h"

#include <algorithm>

MarmotElement::~MarmotElement() {}

void MarmotElement::assignProperty( const ElementProperties& property ) {}

void MarmotElement::assignProperty( const MarmotMaterialSection& property ) {}

void MarmotElement::computeYourself( const double* QTotal,
                                     const double* dQ,
                                     double*       Pint,
                                     double*       K,
                                     const double* time,
                                     double        dT,
                                     double&       pNewDT )
{
  try {
    computeKernels( QTotal, dQ, Pint, K, time[1], dT );
  }
  catch ( const Marmot::StressUpdateFailed& ) {
    pNewDT = std::min( pNewDT, 0.25 );
  }
}

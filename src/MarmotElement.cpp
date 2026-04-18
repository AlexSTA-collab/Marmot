#include "Marmot/MarmotElement.h"
#include "Marmot/Marmot.h"

MarmotElement::~MarmotElement() {}

void MarmotElement::assignProperty( const ElementProperties& property ) {}

void MarmotElement::assignProperty( const MarmotMaterialSection& property ) {}

void MarmotElement::assignMaterial( const std::string& materialName,
                                    const double*      materialProperties,
                                    int                nMaterialProperties )
{
  int code = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( materialName );
  assignProperty( MarmotMaterialSection( code, materialProperties, nMaterialProperties ) );
}

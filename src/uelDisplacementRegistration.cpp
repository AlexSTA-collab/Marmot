#include "MarmotFiniteElement.h"
#include "MarmotFiniteElementUelSpatialWrapper.h"
#include "uelDisplacement.h"
#include "userLibrary.h"

namespace UelDisplacementRegistration {

    template <class T, Marmot::NumIntegration::IntegrationTypes integrationType, typename T::SectionType sectionType>
    userLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
    {
        return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
    }

    using namespace userLibrary;
    using namespace Marmot::NumIntegration;

    const static bool UelCPS4_isRegistered = MarmotElementFactory::
        registerElement( "UELCPS4",
                         userLibrary::ElementCode::UelCPS4,
                         makeFactoryFunction<UelDisplacement<2, 4>,
                                             FullIntegration,
                                             UelDisplacement<2, 4>::PlaneStress>() );

    const static bool UelCPE4_isRegistered = MarmotElementFactory::
        registerElement( "UELCPE4",
                         userLibrary::ElementCode::UelCPE4,
                         makeFactoryFunction<UelDisplacement<2, 4>,
                                             FullIntegration,
                                             UelDisplacement<2, 4>::PlaneStrain>() );

    const static bool UelCPS8R_isRegistered = MarmotElementFactory::
        registerElement( "UELCPS8R",
                         userLibrary::ElementCode::UelCPS8R,
                         makeFactoryFunction<UelDisplacement<2, 8>,
                                             ReducedIntegration,
                                             UelDisplacement<2, 8>::PlaneStress>() );

    const static bool UelCPE8R_isRegistered = MarmotElementFactory::
        registerElement( "UELCPE8R",
                         userLibrary::ElementCode::UelCPE8R,
                         makeFactoryFunction<UelDisplacement<2, 8>,
                                             ReducedIntegration,
                                             UelDisplacement<2, 8>::PlaneStrain>() );

    const static bool UelC3D8_isRegistered = userLibrary::MarmotElementFactory::
        registerElement( "UELC3D8", userLibrary::ElementCode::UelC3D8, []( int elementID ) -> MarmotElement* {
            return new UelDisplacement<3, 8>( elementID,
                                              Marmot::NumIntegration::IntegrationTypes::FullIntegration,
                                              UelDisplacement<3, 8>::SectionType::Solid );
        } );

    const static bool UelC3D20_isRegistered = userLibrary::MarmotElementFactory::
        registerElement( "UELC3D20", userLibrary::ElementCode::UelC3D20, []( int elementID ) -> MarmotElement* {
            return new UelDisplacement<3, 20>( elementID,
                                               Marmot::NumIntegration::IntegrationTypes::FullIntegration,
                                               UelDisplacement<3, 20>::SectionType::Solid );
        } );

    const static bool UelC3D20R_isRegistered = userLibrary::MarmotElementFactory::
        registerElement( "UELC3D20R", userLibrary::ElementCode::UelC3D20R, []( int elementID ) -> MarmotElement* {
            return new UelDisplacement<3, 20>( elementID,
                                               Marmot::NumIntegration::IntegrationTypes::ReducedIntegration,
                                               UelDisplacement<3, 20>::SectionType::Solid );
        } );

    MarmotElement* generateUelT2D2( int elementID )
    {
        auto uelT2D2 = std::unique_ptr<MarmotElement>(
            new UelDisplacement<1, 2>( elementID,
                                       Marmot::NumIntegration::IntegrationTypes::FullIntegration,
                                       UelDisplacement<1, 2>::SectionType::UniaxialStress ) );
        constexpr static int indicesToBeWrapped[] = {0, 1};
        constexpr static int nIndicesToBeWrapped  = 2;
        return new MarmotElementSpatialWrapper( 2,
                                             1,
                                             2,
                                             2,
                                             indicesToBeWrapped,
                                             nIndicesToBeWrapped,
                                             std::move( uelT2D2 ) );
    }
    const static bool
        UelT2D2_isRegistered = userLibrary::MarmotElementFactory::registerElement( "UELT2D2",
                                                                                userLibrary::ElementCode::UelT2D2,
                                                                                generateUelT2D2 );
} // namespace UelDisplacementRegistration

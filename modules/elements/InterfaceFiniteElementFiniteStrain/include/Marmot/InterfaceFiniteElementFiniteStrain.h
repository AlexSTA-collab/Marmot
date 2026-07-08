#pragma once

#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotInterfaceMaterialFiniteStrain.h"

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  class InterfaceFiniteElementFiniteStrain : public InterfaceFiniteElement< nDim, nNodes > {
  public:
    using Base = InterfaceFiniteElement< nDim, nNodes >;

    using typename Base::QuadraturePoint;
    using typename Base::SectionType;

    using typename Base::ForceSized;
    using typename Base::InterfaceDisplSized;
    using typename Base::InterfaceSurfaceGradSized;
    using typename Base::KeSizedMatrix;
    using typename Base::QMatrixSized;
    using typename Base::RhsSized;
    using typename Base::SurfaceStressSized;
    using typename Base::ZMatrixSized;
    using HAverageJumpMatrixSized = Eigen::Matrix< double, Base::nTensor, nDim, Eigen::RowMajor >;
    using HMatrixSized            = Eigen::Matrix< double, nDim, Base::nTensor, Eigen::RowMajor >;

    using Material = MarmotInterfaceMaterialFiniteStrain;

    static constexpr int nTensor = Base::nTensor;

    using Base::assignProperty;

    InterfaceFiniteElementFiniteStrain( int                                         elementID,
                                        FiniteElement::Quadrature::IntegrationTypes integrationType,
                                        SectionType sectionType = SectionType::Interface )
      : Base( elementID, integrationType, sectionType )
    {
      materials.resize( this->qps.size() );
    }

    int getNumberOfRequiredStateVars()
    {
      if ( materials.empty() || !materials[0] ) {
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": no finite-strain interface material assigned." );
      }

      return ( Base::QuadraturePoint::QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
               materials[0]->getNumberOfRequiredStateVars() ) *
             this->qps.size();
    }

    void assignProperty( const MarmotMaterialSection& section )
    {
      assignMaterial( section.materialName, section.materialProperties, section.nMaterialProperties );
    }

    void assignMaterial( const std::string& materialName, const double* materialProperties, int nMaterialProperties )
    {
      for ( auto& material : materials ) {
        material = std::unique_ptr< Material >(
          MarmotLibrary::MarmotInterfaceMaterialFiniteStrainFactory::createMaterial( materialName,
                                                                                     materialProperties,
                                                                                     nMaterialProperties,
                                                                                     this->elLabel ) );
      }
    }

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = this->qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      if ( stateName == "sdv" ) {
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": named state views are not implemented for finite-strain "
                                                   "interface materials." );
    }

    void setInitialConditions( MarmotElement::StateTypes state, const double* values )
    {
      (void)values;

      switch ( state ) {
      case MarmotElement::MarmotMaterialInitialization: {
        for ( size_t q = 0; q < this->qps.size(); ++q ) {
          materials[q]->initializeYourself( this->qps[q].managedStateVars->materialStateVars.data(),
                                            this->qps[q].managedStateVars->materialStateVars.size() );
        }
        break;
      }

      case MarmotElement::MarmotMaterialStateVars: {
        throw std::invalid_argument( "Please use initializeStateVars directly on material" );
      }

      default:
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition for InterfaceFiniteElementFiniteStrain" );
      }
    }

    void computeYourself( const double* QTotal_,
                          const double* dQ_,
                          double*       Pe_,
                          double*       Ke_,
                          const double* time,
                          double        dT,
                          double&       pNewDT )
    {
      (void)QTotal_;

      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );

      constexpr int halfSize = nNodes * nDim / 2;

      for ( size_t q = 0; q < this->qps.size(); ++q ) {
        QuadraturePoint& qp = this->qps[q];

        const auto& Nside = qp.NmatSide;
        const auto& Bside = qp.BmatSide;
        const auto& Njump = qp.NmatJump;
        const auto& Bavg  = qp.BmatAverage;

        const auto dQBottom = dQ.template segment< halfSize >( 0 );
        const auto dQTop    = dQ.template segment< halfSize >( halfSize );

        InterfaceDisplSized dU_GPs;
        dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
        dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

        InterfaceSurfaceGradSized dSurface_strain_GPs;
        dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
        dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

        ForceSized         force          = qp.managedStateVars->force;
        SurfaceStressSized surface_stress = qp.managedStateVars->surfaceStress;

        QMatrixSized            Q_ik;
        HMatrixSized            HJumpAverage;
        HAverageJumpMatrixSized HAverageJump;
        ZMatrixSized            AAverage;

        Q_ik.setZero();
        HJumpAverage.setZero();
        HAverageJump.setZero();
        AAverage.setZero();

        if constexpr ( nDim == 3 ) {
          typename Material::State         materialState{ force.data(),
                                                  surface_stress.data(),
                                                  qp.managedStateVars->materialStateVars.data() };
          typename Material::Tangents      materialTangents{ Q_ik.data(),
                                                        HJumpAverage.data(),
                                                        HAverageJump.data(),
                                                        AAverage.data() };
          typename Material::Deformation   materialDeformation{ dU_GPs.data(),
                                                              dSurface_strain_GPs.data(),
                                                              qp.normal.data() };
          typename Material::TimeIncrement materialTimeIncrement{ time[0], dT };

          try {
            materials[q]->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );
          }
          catch ( const Marmot::StressUpdateFailed& ) {
            pNewDT = 0.5;
            return;
          }
        }
        else if constexpr ( nDim == 2 ) {
          Eigen::Vector3d                                force3d = Eigen::Vector3d::Zero();
          Eigen::Matrix< double, 9, 1 >                  surfaceStress3d;
          Eigen::Matrix< double, 6, 1 >                  dU3d;
          Eigen::Matrix< double, 18, 1 >                 dSurfaceStrain3d;
          Eigen::Vector3d                                normal3d = Eigen::Vector3d::Zero();
          Eigen::Matrix< double, 3, 3, Eigen::RowMajor > Q3d;
          Eigen::Matrix< double, 3, 9, Eigen::RowMajor > HJumpAverage3d;
          Eigen::Matrix< double, 9, 3, Eigen::RowMajor > HAverageJump3d;
          Eigen::Matrix< double, 9, 9, Eigen::RowMajor > AAverage3d;

          surfaceStress3d.setZero();
          dU3d.setZero();
          dSurfaceStrain3d.setZero();
          Q3d.setZero();
          HJumpAverage3d.setZero();
          HAverageJump3d.setZero();
          AAverage3d.setZero();

          for ( int i = 0; i < nDim; ++i ) {
            force3d( i )  = force( i );
            normal3d( i ) = qp.normal( i );
            dU3d( i )     = dU_GPs( i );
            dU3d( 3 + i ) = dU_GPs( nDim + i );

            for ( int j = 0; j < nDim; ++j ) {
              const int index2d = i * nDim + j;
              const int index3d = i * 3 + j;

              surfaceStress3d( index3d )      = surface_stress( index2d );
              dSurfaceStrain3d( index3d )     = dSurface_strain_GPs( index2d );
              dSurfaceStrain3d( 9 + index3d ) = dSurface_strain_GPs( nTensor + index2d );
            }
          }

          typename Material::State         materialState{ force3d.data(),
                                                  surfaceStress3d.data(),
                                                  qp.managedStateVars->materialStateVars.data() };
          typename Material::Tangents      materialTangents{ Q3d.data(),
                                                        HJumpAverage3d.data(),
                                                        HAverageJump3d.data(),
                                                        AAverage3d.data() };
          typename Material::Deformation   materialDeformation{ dU3d.data(), dSurfaceStrain3d.data(), normal3d.data() };
          typename Material::TimeIncrement materialTimeIncrement{ time[0], dT };

          try {
            materials[q]->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );
          }
          catch ( const Marmot::StressUpdateFailed& ) {
            pNewDT = 0.5;
            return;
          }

          for ( int i = 0; i < nDim; ++i ) {
            force( i ) = force3d( i );

            for ( int j = 0; j < nDim; ++j ) {
              const int index2d = i * nDim + j;
              const int index3d = i * 3 + j;

              surface_stress( index2d ) = surfaceStress3d( index3d );
              Q_ik( i, j )              = Q3d( i, j );

              for ( int k = 0; k < nDim; ++k ) {
                const int tensorCol2d = j * nDim + k;
                const int tensorCol3d = j * 3 + k;

                HJumpAverage( i, tensorCol2d ) = HJumpAverage3d( i, tensorCol3d );
                HAverageJump( index2d, k )     = HAverageJump3d( index3d, k );

                for ( int l = 0; l < nDim; ++l ) {
                  const int tensorRow2d  = i * nDim + j;
                  const int tensorRow3d  = i * 3 + j;
                  const int tensorCol2d4 = k * nDim + l;
                  const int tensorCol3d4 = k * 3 + l;

                  AAverage( tensorRow2d, tensorCol2d4 ) = AAverage3d( tensorRow3d, tensorCol3d4 );
                }
              }
            }
          }
        }

        qp.managedStateVars->force         = force;
        qp.managedStateVars->surfaceStress = surface_stress;
        qp.managedStateVars->displacement += dU_GPs;
        qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

        Pe -= Njump.transpose() * force * qp.J0xW;
        Pe -= Bavg.transpose() * surface_stress * qp.J0xW;

        Ke += ( Njump.transpose() * Q_ik * Njump + Njump.transpose() * HJumpAverage * Bavg +
                Bavg.transpose() * HAverageJump * Njump + Bavg.transpose() * AAverage * Bavg ) *
              qp.J0xW;
      }
    }

  private:
    std::vector< std::unique_ptr< Material > > materials;
  };

} // namespace Marmot::Elements

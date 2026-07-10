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
      Eigen::Map< const RhsSized > QTotal( QTotal_ );
      Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );

      constexpr int halfSize = nNodes * nDim / 2;

      for ( size_t q = 0; q < this->qps.size(); ++q ) {
        QuadraturePoint& qp = this->qps[q];

        const auto& Nside = qp.NmatSide;
        const auto& Bside = qp.BmatSide;
        const auto& Njump = qp.NmatJump;
        const auto& Bavg  = qp.BmatAverage;

        const auto QTotalBottom = QTotal.template segment< halfSize >( 0 );
        const auto QTotalTop    = QTotal.template segment< halfSize >( halfSize );

        InterfaceDisplSized totalU_GPs;
        totalU_GPs.template segment< nDim >( 0 )    = Nside * QTotalTop;
        totalU_GPs.template segment< nDim >( nDim ) = Nside * QTotalBottom;

        InterfaceSurfaceGradSized totalSurfaceGradient_GPs;
        totalSurfaceGradient_GPs.template segment< nTensor >( 0 )       = Bside * QTotalTop;
        totalSurfaceGradient_GPs.template segment< nTensor >( nTensor ) = Bside * QTotalBottom;

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
          typename Material::Deformation   materialDeformation{ totalU_GPs.data(),
                                                              totalSurfaceGradient_GPs.data(),
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

          using namespace Fastor;
          using namespace Marmot::FastorIndices;
          using namespace Marmot::FastorStandardTensors;

          Tensor< double, 3, nDim > E( 0.0 );
          E( 0, 0 ) = 1.0;
          E( 1, 1 ) = 1.0;

          const Tensor< double, nDim >       forceTensor( force.data() );
          const Tensor< double, nDim >       normalTensor( qp.normal.data() );
          const Tensor< double, nDim >       topDisplacementTensor( totalU_GPs.data() );
          const Tensor< double, nDim >       bottomDisplacementTensor( totalU_GPs.data() + nDim );
          const Tensor< double, nDim, nDim > surfaceStressTensor( surface_stress.data() );
          const Tensor< double, nDim, nDim > topSurfaceGradientTensor( totalSurfaceGradient_GPs.data() );
          const Tensor< double, nDim, nDim > bottomSurfaceGradientTensor( totalSurfaceGradient_GPs.data() + nTensor );

          const Tensor3d  forceEmbedded3dTensor              = einsum< Ii, i, to_I >( E, forceTensor );
          const Tensor3d  normalEmbedded3dTensor             = einsum< Ii, i, to_I >( E, normalTensor );
          const Tensor3d  topDisplacementEmbedded3dTensor    = einsum< Ii, i, to_I >( E, topDisplacementTensor );
          const Tensor3d  bottomDisplacementEmbedded3dTensor = einsum< Ii, i, to_I >( E, bottomDisplacementTensor );
          const Tensor33d surfaceStressEmbedded3dTensor      = einsum< Ii, Jj, ij, to_IJ >( E, E, surfaceStressTensor );
          const Tensor33d topSurfaceGradientEmbedded3dTensor = einsum< Ii, Jj, ij, to_IJ >( E,
                                                                                            E,
                                                                                            topSurfaceGradientTensor );
          const Tensor33d
            bottomSurfaceGradientEmbedded3dTensor = einsum< Ii, Jj, ij, to_IJ >( E, E, bottomSurfaceGradientTensor );

          force3d         = Marmot::mapEigenToFastor( forceEmbedded3dTensor );
          normal3d        = Marmot::mapEigenToFastor( normalEmbedded3dTensor );
          surfaceStress3d = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( surfaceStressEmbedded3dTensor.data() );
          Eigen::Map< Eigen::Vector3d >( dU3d.data() ) = Marmot::mapEigenToFastor( topDisplacementEmbedded3dTensor );
          Eigen::Map< Eigen::Vector3d >( dU3d.data() +
                                         3 )           = Marmot::mapEigenToFastor( bottomDisplacementEmbedded3dTensor );
          Eigen::Map< Eigen::Matrix< double, 9, 1 > >( dSurfaceStrain3d.data() ) = Eigen::Map<
            const Eigen::Matrix< double, 9, 1 > >( topSurfaceGradientEmbedded3dTensor.data() );
          Eigen::Map< Eigen::Matrix< double, 9, 1 > >(
            dSurfaceStrain3d.data() +
            9 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( bottomSurfaceGradientEmbedded3dTensor.data() );

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

          const TensorMap3d    force3dTensor( force3d.data() );
          const TensorMap33d   surfaceStress3dTensor( surfaceStress3d.data() );
          const TensorMap33d   Q3dTensor( Q3d.data() );
          const TensorMap333d  HJumpAverage3dTensor( HJumpAverage3d.data() );
          const TensorMap333d  HAverageJump3dTensor( HAverageJump3d.data() );
          const TensorMap3333d AAverage3dTensor( AAverage3d.data() );

          const Tensor< double, nDim > forceReducedTensor = einsum< Ii, I, to_i >( E, force3dTensor );
          const Tensor< double, nDim, nDim >
            surfaceStressReducedTensor               = einsum< Ii, Jj, IJ, to_ij >( E, E, surfaceStress3dTensor );
          const Tensor< double, nDim, nDim > QTensor = einsum< Ii, Jj, IJ, to_ij >( E, E, Q3dTensor );
          const Tensor< double, nDim, nDim, nDim >
            HJumpAverageTensor = einsum< Ii, Jj, Kk, IJK, to_ijk >( E, E, E, HJumpAverage3dTensor );
          const Tensor< double, nDim, nDim, nDim >
            HAverageJumpTensor = einsum< Ii, Jj, Kk, IJK, to_ijk >( E, E, E, HAverageJump3dTensor );
          const Tensor< double, nDim, nDim, nDim, nDim >
            AAverageTensor = einsum< Ii, Jj, Kk, Ll, IJKL, to_ijkl >( E, E, E, E, AAverage3dTensor );

          force          = Marmot::mapEigenToFastor( forceReducedTensor );
          surface_stress = Eigen::Map< const SurfaceStressSized >( surfaceStressReducedTensor.data() );
          Q_ik           = Marmot::mapEigenToFastor( QTensor );
          HJumpAverage   = Eigen::Map< const HMatrixSized >( HJumpAverageTensor.data() );
          HAverageJump   = Eigen::Map< const HAverageJumpMatrixSized >( HAverageJumpTensor.data() );
          AAverage       = Eigen::Map< const ZMatrixSized >( AAverageTensor.data() );
        }

        qp.managedStateVars->force         = force;
        qp.managedStateVars->surfaceStress = surface_stress;
        qp.managedStateVars->displacement  = totalU_GPs;
        qp.managedStateVars->surfaceStrain = totalSurfaceGradient_GPs;

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

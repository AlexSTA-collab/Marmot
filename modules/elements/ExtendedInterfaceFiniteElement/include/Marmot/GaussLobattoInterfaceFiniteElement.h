/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

/**
 * @file GaussLobattoInterfaceFiniteElement.h
 * @brief Gauss-Lobatto through-thickness interface finite element (GLIQUAD4/GLILINE2).
 *
 * This file defines the templated class `Marmot::Elements::GaussLobattoInterfaceFiniteElement`
 * and its full in-header implementation. It is a clone of
 * `YInterfaceFiniteElement` with only the material type changed: it routes
 * the SAME top/bottom surface-gradient B-matrices to
 * `MarmotGaussLobattoInterfaceMaterialHypoElastic`, which internally
 * resolves a fixed five-point Gauss-Lobatto through-thickness integration
 * (five independent station materials, a shared-traction local
 * equilibrium problem, static condensation of the per-station normal
 * gradients) and returns a generally FULLY POPULATED 3x3 block tangent
 * over (w, A+, A-) -- including the A+/A- cross-coupling blocks. No
 * element-level kinematics change relative to YIQUAD4/XIQUAD4: the
 * through-thickness resolution is entirely a material-level concern.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGaussLobattoInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class GaussLobattoInterfaceFiniteElement
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Interface finite element with Gauss-Lobatto through-thickness resolution.
   */
  template < int nDim, int nNodes, int NStations = 5 >
  class GaussLobattoInterfaceFiniteElement : public MarmotElement,
                                             public MarmotGeometryInterfaceElement< nDim, nNodes > {

  public:
    enum SectionType {
      Interface,
    };

    static constexpr int nDofPerNodeU    = nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;

    static constexpr int sizeLoadVector = nNodes * nDim;
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int nSideDofs = nInterfaceNodes * nDim;
    static constexpr int nTensor   = nDim * nDim;

    using ParentGeometryElement = MarmotGeometryInterfaceElement< nDim, nNodes >;

    using XiSized              = typename ParentGeometryElement::XiSized;
    using NSized               = typename ParentGeometryElement::NSized;
    using dNdXiSized           = typename ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized = typename ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized          = typename ParentGeometryElement::MetricSized;
    using GradSized            = typename ParentGeometryElement::GradSized;

    using VectorDim = typename ParentGeometryElement::VectorDim;
    using TensorDim = typename ParentGeometryElement::TensorDim;

    using NMatrixSized     = typename ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized = typename ParentGeometryElement::NJumpMatrixSized;

    using BSurfaceSized    = typename ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized = typename ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using QwwMatrixSized = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using QwAMatrixSized = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using QAwMatrixSized = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;
    using QAAMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;

    using Material = MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;

      double detJ;
      double sqrtDetG;
      double J0xW;

      NSized               N;
      dNdXiSized           dNdXi;
      SurfaceJacobianSized J;
      MetricSized          G;
      GradSized            gradN;

      VectorDim normal;
      TensorDim normalProjection;
      TensorDim tangentProjection;

      VectorDim separationVector;
      VectorDim tangentialSeparation;
      double    normalSeparation;

      NMatrixSized  NmatSide;
      BSurfaceSized BmatSide;

      NJumpMatrixSized NmatJump;
      BAvgSurfaceSized BmatAverage;

      /**
       * @brief Named state-variable manager for interface quadrature points.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = nDim },
          { .name = "surfaceStressPlus", .length = nDim * nDim },
          { .name = "surfaceStressMinus", .length = nDim * nDim },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surfaceStrain", .length = 2 * nDim * nDim },
          { .name   = "state block alignment padding",
            .length = ( 4 - ( ( nDim + nDim * nDim + nDim * nDim + 2 * nDim + 2 * nDim * nDim ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< ForceSized >                generalizedForce;
        Eigen::Map< SurfaceStressSized >        surfaceStressPlus;
        Eigen::Map< SurfaceStressSized >        surfaceStressMinus;
        Eigen::Map< InterfaceDisplSized >       displacement;
        Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            generalizedForce( &find( "generalizedForce" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surfaceStrain" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;
      std::unique_ptr< Material >          material;

      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      }

      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          sqrtDetG( 0.0 ),
          J0xW( 0.0 ),
          N( NSized::Zero() ),
          dNdXi( dNdXiSized::Zero() ),
          J( SurfaceJacobianSized::Zero() ),
          G( MetricSized::Zero() ),
          gradN( GradSized::Zero() ),
          normal( VectorDim::Zero() ),
          normalProjection( TensorDim::Zero() ),
          tangentProjection( TensorDim::Zero() ),
          separationVector( VectorDim::Zero() ),
          tangentialSeparation( VectorDim::Zero() ),
          normalSeparation( 0.0 ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    GaussLobattoInterfaceFiniteElement( int                                         elementID,
                                        FiniteElement::Quadrature::IntegrationTypes integrationType,
                                        SectionType sectionType = SectionType::Interface );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape()
    {
      if constexpr ( nDim == 3 && nNodes == 8 ) {
        return "hexa8";
      }
      else if constexpr ( nDim == 2 && nNodes == 4 ) {
        return "bar2";
      }
      else {
        return ParentGeometryElement::getElementShape();
      }
    }

    void assignStateVars( double* stateVars, int nStateVars );

    void assignProperty( const ElementProperties& marmotElementProperty );

    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    void assignMaterial( const std::string& materialName, const double* materialProperties, int nMaterialProperties );

    void assignNodeCoordinates( const double* coordinates );

    void initializeYourself();

    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
                                 double                              dT );

    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    void computeConsistentInertia( double* M );

    void computeLumpedInertia( double* M );

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      if ( stateName == "sdv" ) {
        std::cout << __PRETTY_FUNCTION__ << " on 'sdv' is discouraged and deprecated, please use precise state name";
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

  ///@{

  template < int nDim, int nNodes, int NStations >
  GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::GaussLobattoInterfaceFiniteElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    const auto qpInfos = FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType );

    for ( const auto& qpInfo : qpInfos ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes, int NStations >
  int GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes, int NStations >
  std::vector< std::vector< std::string > > GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::
    getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }
    }

    return nodeFields;
  }

  template < int nDim, int nNodes, int NStations >
  std::vector< int > GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::assignStateVars( double* stateVars,
                                                                                       int     nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::assignProperty(
    const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( section.materialName,
                                                  section.materialProperties,
                                                  section.nMaterialProperties,
                                                  elLabel );
    }
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::assignMaterial( const std::string& materialName,
                                                                                      const double* materialProperties,
                                                                                      int nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( materialName, materialProperties, nMaterialProperties, elLabel );
    }
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const auto geom = this->evaluateAt( qp.xi, 0 );

      qp.N                 = geom.N;
      qp.dNdXi             = geom.dNdXi;
      qp.J                 = geom.J;
      qp.G                 = geom.G;
      qp.sqrtDetG          = geom.sqrtDetG;
      qp.detJ              = geom.sqrtDetG;
      qp.gradN             = geom.gradN;
      qp.normal            = geom.n;
      qp.normalProjection  = geom.normalProjection;
      qp.tangentProjection = geom.tangentProjection;

      qp.NmatSide    = geom.NmatSide;
      qp.BmatSide    = geom.BmatSide;
      qp.NmatJump    = geom.NmatJump;
      qp.BmatAverage = geom.BmatAverage;

      const VectorDim xBottom = qp.NmatSide * this->getSideCoordinates( 0 );
      const VectorDim xTop    = qp.NmatSide * this->getSideCoordinates( 1 );
      qp.separationVector     = xTop - xBottom;

      constexpr double geometryTolerance = 1.0e-12;
      if ( qp.separationVector.norm() > geometryTolerance && qp.separationVector.dot( qp.normal ) < 0.0 ) {
        qp.normal *= -1.0;
      }

      qp.normalProjection     = qp.normal * qp.normal.transpose();
      qp.tangentProjection    = TensorDim::Identity() - qp.normalProjection;
      qp.normalSeparation     = qp.separationVector.dot( qp.normal );
      qp.tangentialSeparation = qp.tangentProjection * qp.separationVector;

      if ( qp.separationVector.norm() > geometryTolerance && qp.normalSeparation <= geometryTolerance ) {
        throw std::invalid_argument( "GaussLobattoInterfaceFiniteElement: paired faces have no positive separation in "
                                     "the interface-normal direction." );
      }

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      if ( qp.material ) {
        if constexpr ( nDim == 3 ) {
          qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
        }
        else if constexpr ( nDim == 2 ) {
          qp.material->setCharacteristicElementLength( qp.sqrtDetG );
        }
      }
    }
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::computeKernels( const double* QTotal_,
                                                                                      const double* dQ_,
                                                                                      double*       Pe_,
                                                                                      double*       Ke_,
                                                                                      double        time,
                                                                                      double        dT )
  {
    (void)QTotal_;

    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    constexpr int halfSize = nNodes * nDim / 2;

    for ( QuadraturePoint& qp : qps ) {
      const auto& Nside = qp.NmatSide;
      const auto& Bside = qp.BmatSide;
      const auto& Njump = qp.NmatJump;

      BAvgSurfaceSized BPlusFull                                   = BAvgSurfaceSized::Zero();
      BAvgSurfaceSized BMinusFull                                  = BAvgSurfaceSized::Zero();
      BMinusFull.template block< nTensor, halfSize >( 0, 0 )       = Bside;
      BPlusFull.template block< nTensor, halfSize >( 0, halfSize ) = Bside;

      const auto dQBottom = dQ.template segment< halfSize >( 0 );
      const auto dQTop    = dQ.template segment< halfSize >( halfSize );

      InterfaceDisplSized dU_GPs;
      dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
      dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

      InterfaceSurfaceGradSized dSurface_strain_GPs;
      dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
      dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

      ForceSized         generalizedForce   = qp.managedStateVars->generalizedForce;
      SurfaceStressSized surfaceStressPlus  = qp.managedStateVars->surfaceStressPlus;
      SurfaceStressSized surfaceStressMinus = qp.managedStateVars->surfaceStressMinus;

      QwwMatrixSized Q_ww;
      QwAMatrixSized Q_wAp, Q_wAm;
      QAwMatrixSized Q_Apw, Q_Amw;
      QAAMatrixSized Q_ApAp, Q_ApAm, Q_AmAp, Q_AmAm;

      Q_ww.setZero();
      Q_wAp.setZero();
      Q_wAm.setZero();
      Q_Apw.setZero();
      Q_Amw.setZero();
      Q_ApAp.setZero();
      Q_ApAm.setZero();
      Q_AmAp.setZero();
      Q_AmAm.setZero();

      if constexpr ( nDim == 3 ) {
        typename Material::State         materialState{ generalizedForce.data(),
                                                surfaceStressPlus.data(),
                                                surfaceStressMinus.data(),
                                                qp.managedStateVars->materialStateVars.data() };
        typename Material::Tangents      materialTangents{ Q_ww.data(),
                                                      Q_wAp.data(),
                                                      Q_wAm.data(),
                                                      Q_Apw.data(),
                                                      Q_ApAp.data(),
                                                      Q_ApAm.data(),
                                                      Q_Amw.data(),
                                                      Q_AmAp.data(),
                                                      Q_AmAm.data() };
        typename Material::Deformation   materialDeformation{ dU_GPs.data(),
                                                            dSurface_strain_GPs.data(),
                                                            qp.normal.data(),
                                                            qp.separationVector.data() };
        typename Material::TimeIncrement materialTimeIncrement{ time, dT };

        qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );
      }
      else if constexpr ( nDim == 2 ) {
        Eigen::Vector3d                                force3d = Eigen::Vector3d::Zero();
        Eigen::Matrix< double, 9, 1 >                  surfaceStressPlus3d;
        Eigen::Matrix< double, 9, 1 >                  surfaceStressMinus3d;
        Eigen::Matrix< double, 6, 1 >                  dU3d;
        Eigen::Matrix< double, 18, 1 >                 dSurfaceStrain3d;
        Eigen::Vector3d                                normal3d = Eigen::Vector3d::Zero();
        Eigen::Matrix< double, 3, 3, Eigen::RowMajor > Qww3d;
        Eigen::Matrix< double, 3, 9, Eigen::RowMajor > QwAp3d, QwAm3d;
        Eigen::Matrix< double, 9, 3, Eigen::RowMajor > QApw3d, QAmw3d;
        Eigen::Matrix< double, 9, 9, Eigen::RowMajor > QApAp3d, QApAm3d, QAmAp3d, QAmAm3d;
        Eigen::Vector3d                                separation3d = Eigen::Vector3d::Zero();

        surfaceStressPlus3d.setZero();
        surfaceStressMinus3d.setZero();
        dU3d.setZero();
        dSurfaceStrain3d.setZero();
        Qww3d.setZero();
        QwAp3d.setZero();
        QwAm3d.setZero();
        QApw3d.setZero();
        QAmw3d.setZero();
        QApAp3d.setZero();
        QApAm3d.setZero();
        QAmAp3d.setZero();
        QAmAm3d.setZero();

        for ( int i = 0; i < nDim; ++i ) {
          force3d( i )      = generalizedForce( i );
          normal3d( i )     = qp.normal( i );
          separation3d( i ) = qp.separationVector( i );
          dU3d( i )         = dU_GPs( i );
          dU3d( 3 + i )     = dU_GPs( nDim + i );

          for ( int j = 0; j < nDim; ++j ) {
            const int index2d = i * nDim + j;
            const int index3d = i * 3 + j;

            surfaceStressPlus3d( index3d )  = surfaceStressPlus( index2d );
            surfaceStressMinus3d( index3d ) = surfaceStressMinus( index2d );
            dSurfaceStrain3d( index3d )     = dSurface_strain_GPs( index2d );
            dSurfaceStrain3d( 9 + index3d ) = dSurface_strain_GPs( nTensor + index2d );
          }
        }

        typename Material::State         materialState{ force3d.data(),
                                                surfaceStressPlus3d.data(),
                                                surfaceStressMinus3d.data(),
                                                qp.managedStateVars->materialStateVars.data() };
        typename Material::Tangents      materialTangents{ Qww3d.data(),
                                                      QwAp3d.data(),
                                                      QwAm3d.data(),
                                                      QApw3d.data(),
                                                      QApAp3d.data(),
                                                      QApAm3d.data(),
                                                      QAmw3d.data(),
                                                      QAmAp3d.data(),
                                                      QAmAm3d.data() };
        typename Material::Deformation   materialDeformation{ dU3d.data(),
                                                            dSurfaceStrain3d.data(),
                                                            normal3d.data(),
                                                            separation3d.data() };
        typename Material::TimeIncrement materialTimeIncrement{ time, dT };

        qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        for ( int i = 0; i < nDim; ++i ) {
          generalizedForce( i ) = force3d( i );
          for ( int k = 0; k < nDim; ++k ) {
            Q_ww( i, k ) = Qww3d( i, k );
          }
        }

        for ( int i = 0; i < nDim; ++i ) {
          for ( int j = 0; j < nDim; ++j ) {
            const int row2d = i * nDim + j;
            const int row3d = i * 3 + j;

            surfaceStressPlus( row2d )  = surfaceStressPlus3d( row3d );
            surfaceStressMinus( row2d ) = surfaceStressMinus3d( row3d );

            for ( int k = 0; k < nDim; ++k ) {
              Q_Apw( row2d, k ) = QApw3d( row3d, k );
              Q_Amw( row2d, k ) = QAmw3d( row3d, k );
            }

            for ( int k = 0; k < nDim; ++k ) {
              for ( int l = 0; l < nDim; ++l ) {
                const int col2d = k * nDim + l;
                const int col3d = k * 3 + l;

                Q_ApAp( row2d, col2d ) = QApAp3d( row3d, col3d );
                Q_ApAm( row2d, col2d ) = QApAm3d( row3d, col3d );
                Q_AmAp( row2d, col2d ) = QAmAp3d( row3d, col3d );
                Q_AmAm( row2d, col2d ) = QAmAm3d( row3d, col3d );
              }
            }
          }
        }

        for ( int i = 0; i < nDim; ++i ) {
          for ( int k = 0; k < nDim; ++k ) {
            for ( int l = 0; l < nDim; ++l ) {
              const int col2d = k * nDim + l;
              const int col3d = k * 3 + l;

              Q_wAp( i, col2d ) = QwAp3d( i, col3d );
              Q_wAm( i, col2d ) = QwAm3d( i, col3d );
            }
          }
        }
      }

      qp.managedStateVars->generalizedForce   = generalizedForce;
      qp.managedStateVars->surfaceStressPlus  = surfaceStressPlus;
      qp.managedStateVars->surfaceStressMinus = surfaceStressMinus;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      Pe -= Njump.transpose() * generalizedForce * qp.J0xW;
      Pe -= BPlusFull.transpose() * surfaceStressPlus * qp.J0xW;
      Pe -= BMinusFull.transpose() * surfaceStressMinus * qp.J0xW;

      // Full nine-block assembly -- unlike XInterfaceFiniteElement, the
      // cross blocks Q_ApAm/Q_AmAp are generally nonzero after
      // equilibrium condensation and must NOT be omitted.
      Ke += ( Njump.transpose() * Q_ww * Njump + Njump.transpose() * Q_wAp * BPlusFull +
              Njump.transpose() * Q_wAm * BMinusFull + BPlusFull.transpose() * Q_Apw * Njump +
              BPlusFull.transpose() * Q_ApAp * BPlusFull + BPlusFull.transpose() * Q_ApAm * BMinusFull +
              BMinusFull.transpose() * Q_Amw * Njump + BMinusFull.transpose() * Q_AmAp * BPlusFull +
              BMinusFull.transpose() * Q_AmAm * BMinusFull ) *
            qp.J0xW;
    }
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::setInitialConditions( StateTypes    state,
                                                                                            const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      }
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument(
        MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition for GaussLobattoInterfaceFiniteElement" );
    }
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             P,
    double*                             K,
    const int                           elementFace,
    const double*                       load,
    const double*                       QTotal,
    double                              time,
    double                              dT )
  {
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__
                                 << ": distributed loads are not implemented for GaussLobattoInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::computeBodyForce( double*       P,
                                                                                        double*       K,
                                                                                        const double* load,
                                                                                        const double* QTotal,
                                                                                        double        time,
                                                                                        double        dT )
  {
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__
                                 << ": body forces are not implemented for GaussLobattoInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::computeConsistentInertia( double* M )
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__
                                           << ": inertia is not implemented for GaussLobattoInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes, int NStations >
  void GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::computeLumpedInertia( double* M )
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__
                                           << ": inertia is not implemented for GaussLobattoInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes, int NStations >
  std::vector< double > GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< VectorDim > coordsMap( coords.data() );

    const auto centerXi = XiSized::Zero();
    const auto Ncenter  = this->N( centerXi );
    const auto Nmat     = this->NMatrix( Ncenter );

    const auto xSide = this->getSideCoordinates( 0 );

    coordsMap = Nmat * xSide;

    return coords;
  }

  template < int nDim, int nNodes, int NStations >
  std::vector< std::vector< double > > GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::
    getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    for ( const auto& qp : qps ) {
      std::vector< double > coords( nDim );

      Eigen::Map< VectorDim > coordsMap( coords.data() );

      const auto xSide = this->getSideCoordinates( 0 );
      coordsMap        = qp.NmatSide * xSide;

      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes, int NStations >
  int GaussLobattoInterfaceFiniteElement< nDim, nNodes, NStations >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  ///@}

} // namespace Marmot::Elements

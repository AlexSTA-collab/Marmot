#pragma once
#include "MarmotFunctions.h"
#include "MarmotTypedefs.h"

#define isNaN( x ) ( x != x )
#define VOIGTFROMDIM( x ) ( ( ( x * x ) + x ) >> 1 )
//#define DIMFROMVOIGT( x ) ( x<<1  )

namespace Marmot {
    namespace Vgt {

        // TODO: Remove, only valid for 3D!
        const int VoigtSize = 6;

        extern const Marmot::Vector6 P;
        extern const Marmot::Vector6 PInv;

        extern const Marmot::Vector6 I;
        extern const Marmot::Vector6 IHyd;
        extern const Matrix6      IDev;

        // Plane Stress handling
        Eigen::Vector3d voigtToPlaneVoigt( const Marmot::Vector6& voigt );
        Marmot::Vector6    planeVoigtToVoigt( const Eigen::Vector3d& voigtPlane );

        template <int voigtSize>
        Eigen::Matrix<double, voigtSize, 1> reduce3DVoigt( const Marmot::Vector6& Voigt3D )
        {
            if constexpr ( voigtSize == 1 )
                return ( Eigen::Matrix<double, 1, 1>() << Voigt3D( 0 ) ).finished();
            else if constexpr ( voigtSize == 3 )
                return voigtToPlaneVoigt( Voigt3D );
            else if constexpr ( voigtSize == 6 )
                return Voigt3D;
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        template <int voigtSize>
        Marmot::Vector6 make3DVoigt( const Eigen::Matrix<double, voigtSize, 1>& Voigt )
        {
            if constexpr ( voigtSize == 1 )
                return ( Marmot::Vector6() << Voigt( 0 ), 0, 0, 0, 0, 0 ).finished();
            else if constexpr ( voigtSize == 3 )
                return planeVoigtToVoigt( Voigt );
            else if constexpr ( voigtSize == 6 )
                return Voigt;
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        /*compute E33 for a given elastic strain, to compute the compensation for
         * planeStress = Cel : (elasticStrain + compensationStrain) */
        Marmot::Vector6 planeStressCompensationStrain( const Marmot::Vector6& elasticStrain, double nu );
        /* Returns the transformation Matrix T which fullfills
         * planeStressIncrement = C : (T : arbitraryStrainIncrement) */
        Matrix6                     planeStressTangentTransformationMatrix( const Matrix6& tangent );
        Eigen::Matrix<double, 6, 3> dStrainDStrainPlaneStress( const Matrix6& tangent );
        Eigen::Matrix<double, 3, 6> dStressPlaneStressDStress();
        Eigen::Matrix<double, 6, 3> dStrainDStrainPlaneStrain();

        // function prototypes for  Marmot::Vector6 handling
        Eigen::Matrix3d voigtToStrain( const Marmot::Vector6& strainVector );
        Eigen::Matrix3d voigtToStress( const Marmot::Vector6& stressVector );
        Marmot::Vector6    strainToVoigt( const Eigen::Matrix3d& strainTensor );
        Marmot::Vector6    stressToVoigt( const Eigen::Matrix3d& stressTensor );
        Eigen::Vector3d haighWestergaard( const Marmot::Vector6& stress );
        Eigen::Vector3d haighWestergaardStrain( const Marmot::Vector6& strain );

        template <int nDim>
        Eigen::Matrix<double, nDim, nDim> StressMatrixFromVoigt(
            const Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1>& Voigt )
        {
            if constexpr ( nDim == 1 )
                return ( Eigen::Matrix<double, nDim, nDim>() << Voigt( 0 ) ).finished();
            else if constexpr ( nDim == 2 )
                return ( Eigen::Matrix<double, nDim, nDim>() << Voigt( 0 ), Voigt( 2 ), Voigt( 2 ), Voigt( 1 ) )
                    .finished();
            else if constexpr ( nDim == 3 )
                return voigtToStress( Voigt );
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        template <int nDim>
        Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1> VoigtFromStrainMatrix(
            const Eigen::Matrix<double, nDim, nDim>& strain )
        {
            if constexpr ( nDim == 1 )
                return ( Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1>() << strain( 0, 0 ) ).finished();
            else if constexpr ( nDim == 2 )
                return ( Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1>() << strain( 0, 0 ),
                         strain( 1, 1 ),
                         2 * strain( 0, 1 ) )
                    .finished();
            else if constexpr ( nDim == 3 )
                return strainToVoigt( strain );
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        // principal strains calculated by solving eigenvalue problem ( !NOT sorted! )
        Eigen::Vector3d principalStrains( const Marmot::Vector6& strain );
        // principal strains calculated from haigh westergaard strains ( sorted --> e1 > e2 > e3 )
        Eigen::Vector3d principalStrainsHW( const Marmot::Vector6& strain );
        // principal stresses calculated by solving eigenvalue problem ( !NOT sorted! )
        Eigen::Vector3d principalStresses( const Marmot::Vector6& stress );
        // principal stressDirections calculated by solving eigenvalue problem ( !NOT sorted! )
        Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6& stress );
        // rotate a 2nd order tensor T in voigt notation by : T' = Q * T * Q^T
        Marmot::Vector6 rotateVoigtStress( const Eigen::Matrix3d Q, const Marmot::Vector6& stress );

        // equivalent von Mises stress
        double vonMisesEquivalentStress( const Marmot::Vector6& stress );
        // equivalent von Mises strain
        double vonMisesEquivalentStrain( const Marmot::Vector6& strain );
        // Euclidian norm of strain
        double normStrain( const Marmot::Vector6& strain );
        // Euclidian norm of stress
        double normStress( const Marmot::Vector6& stress );
        // Trace of compressive strains
        double Evolneg( const Marmot::Vector6& strain );

        // Invariants - keep atention: different for stress/strain tensor
        double I1( const Marmot::Vector6& stress );
        double I2( const Marmot::Vector6& stress );
        double I2strain( const Marmot::Vector6& strain );
        double I3( const Marmot::Vector6& stress );
        double I3strain( const Marmot::Vector6& strain );

        // Invariants of the deviatoric part of the stress/strain tensor
        double J2( const Marmot::Vector6& stress );
        double J2strain( const Marmot::Vector6& strain );
        double J3( const Marmot::Vector6& stress );
        double J3strain( const Marmot::Vector6& strain );

        // derivatives of Haigh Westergaard stresses with respect to cauchy stress in eng. notation
        Marmot::Vector6 dSigmaMdSigma();
        Marmot::Vector6 dRhodSigma( double rho, const Marmot::Vector6& stress );
        Marmot::Vector6 dThetadSigma( double theta, const Marmot::Vector6& stress );
        // derivatives of Haigh Westergaard stresses with respect to deviatoric invariants
        double dTheta_dJ2( const Marmot::Vector6& stress );
        double dTheta_dJ3( const Marmot::Vector6& stress );
        // derivatives of Haigh Westergaard strains with respect to deviatoric invariants
        double dThetaE_dJ2E( const Marmot::Vector6& strain );
        double dThetaE_dJ3E( const Marmot::Vector6& strain );
        // derivatives of deviatoric invariants with respect to eng. stresses
        Marmot::Vector6 dJ2_dStress( const Marmot::Vector6& stress );
        Marmot::Vector6 dJ3_dStress( const Marmot::Vector6& stress );
        // derivatives of deviatoric invariants with respect to eng. strains
        Marmot::Vector6 dJ2E_dE( const Marmot::Vector6& strain );
        Marmot::Vector6 dJ3E_dE( const Marmot::Vector6& strain );
        // derivatives of Haigh Westergaard strains with respect to eng. strains
        Marmot::Vector6 dThetaE_dE( const Marmot::Vector6& strain );

        // derivatives of principalStess with respect to stress
        Marmot::Matrix36 dStressPrincipals_dStress( const Marmot::Vector6& stress );

        // derivatives of plastic strains with respect to strains
        Eigen::Vector3d dDeltaEpvneg_dDeltaEpPrincipals( const Marmot::Vector6& strain );
        Matrix6         dEp_dE( const Matrix6& CelInv, const Matrix6& Cep );
        RowVector6d     dDeltaEpv_dE( const Matrix6& CelInv, const Matrix6& Cep );
        Marmot::Matrix36   dDeltaEpPrincipals_dDeltaEp( const Marmot::Vector6& dEp );
        RowVector6d     dDeltaEpvneg_dE( const Marmot::Vector6& dEp, const Matrix6& CelInv, const Matrix6& Cep );

    } // namespace Vgt

    namespace mechanics {
        template <typename T>
        int sgn( T val )
        {
            return ( T( 0 ) < val ) - ( val < T( 0 ) );
        }

        Matrix6         Cel( double E, double nu );
        Matrix6         CelInverse( double E, double nu );
        Eigen::Matrix3d getPlaneStressTangent( const Matrix6& C );

        Eigen::Matrix3d getPlaneStrainTangent( const Matrix6& C );
        double          getUniaxialStressTangent( const Eigen::Ref<const Matrix6>& C );
        double          E( const double K, const double G );
        double          nu( const double K, const double G );

        namespace PlaneStrain {
            Tensor322d dStressdDeformationGradient( const Tensor633d& dStressdDeformationGradient3D );
        }
        namespace PlaneStress {
            Tensor322d dStressdDeformationGradient( const Tensor633d& dStressdDeformationGradient3D );
        }
    } // namespace mechanics
} // namespace Marmot

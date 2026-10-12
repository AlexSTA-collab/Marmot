#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"

#include "Fastor/Fastor.h"

#include <stdexcept>

namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions {

  using namespace Marmot::FastorIndices;
  using namespace Fastor;

  NormalFibre normalFibre( const Tensor6d&  U,
                           const Tensor18d& surfaceGradient,
                           const Tensor3d&  normal,
                           const Tensor3d&  connector )
  {
    NormalFibre fibre;

    const Tensor3d jumpU                      = U( seq( 0, 3 ) ) - U( seq( 3, last ) );
    const Tensor9d averageSurfaceGradientFlat = 0.5 *
                                                ( surfaceGradient( seq( 0, 9 ) ) + surfaceGradient( seq( 9, last ) ) );
    fibre.averageSurfaceGradient = Tensor33d( reshape< 3, 3 >( averageSurfaceGradientFlat ) );

    fibre.ell = inner( connector, normal );

    if ( fibre.ell <= 0.0 )
      throw std::invalid_argument(
        "MarmotInterfaceMaterialFiniteStrain: the paired faces of the interface element must be separated along the "
        "interface normal, and the normal must point from the first side to the second. A coincident-face mesh "
        "carries no normal fibre to differentiate along." );

    fibre.dTau = connector - fibre.ell * normal;

    fibre.normalGradient = ( 1.0 / fibre.ell ) *
                           ( jumpU - einsum< ij, j, to_i >( fibre.averageSurfaceGradient, fibre.dTau ) );

    fibre.F = Tensor33d( 0.0 );
    fibre.F.eye();
    fibre.F += fibre.averageSurfaceGradient + einsum< i, j, to_ij >( fibre.normalGradient, normal );

    return fibre;
  }

  void firstPiolaKirchhoff( const Tensor33d&   tau,
                            const Tensor3333d& dTau_dF,
                            const Tensor33d&   F,
                            Tensor33d&         P,
                            Tensor3333d&       A )
  {
    const Tensor33d FInvT = transpose( inverse( F ) );
    P                     = einsum< ij, jk, to_ik >( tau, FInvT );
    A = einsum< imkl, mj, to_ijkl >( dTau_dF, FInvT ) - einsum< im, kj, ml, to_ijkl >( tau, FInvT, FInvT );
  }

  ReducedResponse reduceToInterface( const Tensor33d&   P,
                                     const Tensor3333d& A,
                                     const Tensor3d&    normal,
                                     const NormalFibre& fibre,
                                     double             h )
  {
    const Tensor3d& dTau = fibre.dTau;
    const double    ell  = fibre.ell;

    Tensor33d T( 0.0 );
    T.eye();
    T -= einsum< i, j, to_ij >( normal, normal );

    const double hOverEll = h / ell;

    const Tensor3d    traction = einsum< ij, j, to_i >( P, normal );
    const Tensor33d   ANN      = einsum< ijkl, j, l, to_ik >( A, normal, normal );
    const Tensor333d  ANT      = einsum< ijkl, j, ml, to_ikm >( A, normal, T );
    const Tensor333d  ATN      = einsum< ijkl, mj, l, to_imk >( A, T, normal );
    const Tensor3333d ATT      = einsum< ijkl, mj, nl, to_imkn >( A, T, T );

    // ANN with the tangential offset appended, in both index orders the tangent blocks need
    const Tensor333d ANNdTauJumpAverage = einsum< ik, m, to_ikm >( ANN, dTau );
    const Tensor333d ANNdTauAverageJump = einsum< ik, m, to_imk >( ANN, dTau );

    ReducedResponse response;
    response.force         = hOverEll * traction;
    response.surfaceStress = h * einsum< ij, kj, to_ik >( P, T ) - hOverEll * einsum< i, k, to_ik >( traction, dTau );

    response.Q_ik   = ( hOverEll / ell ) * ANN;
    response.H_ikR  = hOverEll * ( ANT - ( 1.0 / ell ) * ANNdTauJumpAverage );
    response.H_iRk  = hOverEll * ( ATN - ( 1.0 / ell ) * ANNdTauAverageJump );
    response.A_iRkP = h * ATT -
                      hOverEll * ( einsum< imk, n, to_imkn >( ATN, dTau ) + einsum< m, ikn, to_imkn >( dTau, ANT ) ) +
                      ( hOverEll / ell ) * einsum< imk, n, to_imkn >( ANNdTauAverageJump, dTau );
    return response;
  }

} // namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions

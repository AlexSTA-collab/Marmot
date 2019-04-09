#pragma once
#include "bftTypedefs.h"

/*
 * convinience functions for umats in Abaqus
 * */

namespace bft {
    [[deprecated]]
    void backToAbaqus( const Matrix6&               jacobian,
                       Eigen::Map<Eigen::MatrixXd>& ABQJacobian,
                       const bft::Vector6&          stress,
                       Eigen::Map<Eigen::VectorXd>& ABQStress,
                       int                          nTensor );

    [[deprecated]]
    void backToAbaqusPlaneStress( const Matrix6&               jacobian,
                                  Eigen::Map<Eigen::MatrixXd>& ABQJacobian,
                                  const bft::Vector6&          stress,
                                  Eigen::Map<Eigen::VectorXd>& ABQStress );

    [[deprecated]]
    void backToAbaqusNonLocal( const Matrix6&              dStressdStrain,
                               Eigen::Ref<Eigen::MatrixXd> ABQdStressDStrain,
                               const bft::Vector6&         stress,
                               Eigen::Ref<Eigen::VectorXd> ABQStress,
                               double                      intParameterLocal,
                               double&                     ABQParameterLocal,
                               const bft::Vector6&         dStressDIntParamNonLocal,
                               Eigen::Ref<Eigen::VectorXd> ABQDStressDIntParamNonLocal,
                               const bft::Vector6&         dIntParamLocalDStrain,
                               Eigen::Ref<Eigen::VectorXd> ABQDIntParameterLocalDStrain,
                               double                      nonLocalRadius,
                               double&                     ABQNonLocalRadius,
                               int                         nTensor );

    void discardTheIncrement( double& pNewDT, double value, const std::string& message );
} // namespace bft

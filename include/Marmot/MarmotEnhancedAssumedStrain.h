#pragma once
#include "Marmot/MarmotTypedefs.h"

namespace Marmot {
    namespace FiniteElementTechnology::EAS {

        enum EASType {
            DeBorstEAS2,
            DeBorstEAS2_P2,
            EAS3,
            DeBorstEAS6b,
            DeBorstEAS9,
            SimoRifaiEAS5,
            SimoRifaiEAS4,
        };

        Eigen::MatrixXd F( const Eigen::MatrixXd& J );

        Eigen::MatrixXd EASInterpolation( EASType type, const Eigen::VectorXd& xi );

    } // namespace FiniteElementTechnology::EAS
} // namespace Marmot

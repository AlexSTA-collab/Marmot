#!/usr/bin/env bash
# Build the plug-ins of the IJES revision against an installed Marmot of this branch.
# usage: build.sh <Marmot install prefix> <conda env prefix with Eigen and autodiff> [out dir = ./lib]
# Load them before EdelweissFE starts: LD_PRELOAD=<prefix>/lib/libMarmot.so.1:<out>/libMarmotBBar.so:...
set -e
PFX=$1; ENV=$2; OUT=${3:-$PWD/lib}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
CXX=${CXX:-g++}
FLAGS="-std=c++20 -O3 -fPIC -shared -fopenmp -I$PFX/include -I$ENV/include/eigen3 -I$ENV/include"
build() {  # <dir> <source> <library>
  (cd "$HERE/$1" && $CXX $FLAGS -I. "$2" -L"$PFX/lib" -lMarmot -o "$OUT/$3") && echo "built $3"
}
build bbar          DisplacementBBarRegistration.cpp        libMarmotBBar.so           # C3D8B
build mixed         DisplacementMixedBubbleRegistration.cpp libMarmotMixed.so          # C3D8P
(cd "$HERE/druckerprager" && $CXX $FLAGS -I. DruckerPrager.cpp Registration.cpp -L"$PFX/lib" -lMarmot \
  -o "$OUT/libMarmotDRUCKERPRAGER.so") && echo "built libMarmotDRUCKERPRAGER.so"                      # DRUCKERPRAGER
build psie-avg      Registration.cpp                        libMarmotPSIEAVG.so        # YIQUAD4_STABP_MINI_AVG
build psie-squeeze  Registration.cpp                        libMarmotPSIESQ.so         # YIQUAD4_STABP_MINI_SQ
build twophase      Registration.cpp                        libMarmotTWOPHASE.so       # YIQUAD4_TWOPHASE

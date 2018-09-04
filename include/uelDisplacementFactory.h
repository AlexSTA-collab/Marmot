#pragma once
#include "userLibrary.h"
#include "bftUel.h"
#include "bftTypedefs.h"
#include <string>

namespace UelDisplacementFactory{

    BftUel* generateUelT2D2(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties);
    
    BftUel* generateUelCPS4(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPS8(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPE4(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode umat, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPE8(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPS4R(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPS8R(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPE4R(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelCPE8R(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelC3D8(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelC3D8R(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelC3D20(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

    BftUel* generateUelC3D20R(const double* coordinates, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode hypoElasticMaterialCode, const double* materialProperties, int nMaterialProperties);

}

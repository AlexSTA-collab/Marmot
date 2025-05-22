# Install script for directory: /home/alexsta1993/alexandros/Marmot

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/alexsta1993/miniforge3/envs/alexandros")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMarmot.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMarmot.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMarmot.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/alexsta1993/alexandros/Marmot/build/lib/libMarmot.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMarmot.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMarmot.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMarmot.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/Marmot.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotElement.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotElementProperty.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotJournal.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterial.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialRegistrationHelper.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotTesting.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotUtils.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMath.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotConstants.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotTensor.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotTypedefs.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialGradientEnhancedHypoElastic.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialMechanical.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialHypoElastic.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialHypoElasticAD.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialHyperElastic.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialGradientEnhancedMechanical.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotVoigt.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialHypoElasticInterface.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/MarmotMaterialMechanicalInterface.h;/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot/LinearElasticInterface.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/alexsta1993/miniforge3/envs/alexandros/include/Marmot" TYPE FILE FILES
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/Marmot.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotElement.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotElementProperty.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotJournal.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotMaterial.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotMaterialRegistrationHelper.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotTesting.h"
    "/home/alexsta1993/alexandros/Marmot/include/Marmot/MarmotUtils.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/include/Marmot/MarmotMath.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/include/Marmot/MarmotConstants.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/include/Marmot/MarmotTensor.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/include/Marmot/MarmotTypedefs.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialGradientEnhancedHypoElastic.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialMechanical.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialHypoElastic.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialHypoElasticAD.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialHyperElastic.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialGradientEnhancedMechanical.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotVoigt.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialHypoElasticInterface.h"
    "/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/include/Marmot/MarmotMaterialMechanicalInterface.h"
    "/home/alexsta1993/alexandros/Marmot/modules/materials/LinearElasticInterface/include/Marmot/LinearElasticInterface.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/alexsta1993/alexandros/Marmot/build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/alexsta1993/alexandros/Marmot/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()

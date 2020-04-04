## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Compile and install enabled bundled features:
#

MESSAGE(STATUS "Setting up bundled features")


IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${BOOST_FOLDER}/include/boost
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    )

  ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/serialization/src)
  ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/system/src)

  IF(DEAL_II_WITH_ZLIB)
    ADD_SUBDIRECTORY(${BOOST_FOLDER}/libs/iostreams/src)
  ELSE()
    MESSAGE(STATUS "BOOST::Iostreams will not be available because zlib or its header files could not be found")
  ENDIF()
ENDIF()


IF(FEATURE_THREADS_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${TBB_FOLDER}/include/tbb
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    )

  ADD_SUBDIRECTORY(${TBB_FOLDER}/src)
ENDIF()


IF(FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/UMFPACK/Source)
  ADD_SUBDIRECTORY(${UMFPACK_FOLDER}/AMD/Source)

  INSTALL(DIRECTORY
      ${UMFPACK_FOLDER}/UMFPACK/Include/
      ${UMFPACK_FOLDER}/AMD/Include/
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    )
ENDIF()


IF(FEATURE_MUPARSER_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${MUPARSER_FOLDER})
ENDIF()


IF(FEATURE_MESQUITE_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${MESQUITE_FOLDER})

  #
  # Mesquite writes some files into the build directory, so we need to include
  # those as well.
  #
  STRING(REGEX MATCH "mesquite-([0-9]+)\\.([0-9]+)([0-9]+)" MSQ_CONFIGURE_FOLDER  ${MESQUITE_FOLDER})
  
  INSTALL(DIRECTORY 
      ${CMAKE_BINARY_DIR}/bundled/${MSQ_CONFIGURE_FOLDER}/include
      ${MESQUITE_FOLDER}/src/include
      ${MESQUITE_FOLDER}/src/Control
      ${MESQUITE_FOLDER}/src/MappingFunction
      ${MESQUITE_FOLDER}/src/MappingFunction/Linear
      ${MESQUITE_FOLDER}/src/MappingFunction/Lagrange
      ${MESQUITE_FOLDER}/src/Mesh
      ${MESQUITE_FOLDER}/src/Misc
      ${MESQUITE_FOLDER}/src/ObjectiveFunction
      ${MESQUITE_FOLDER}/src/QualityAssessor
      ${MESQUITE_FOLDER}/src/QualityImprover
      ${MESQUITE_FOLDER}/src/QualityImprover/OptSolvers
      ${MESQUITE_FOLDER}/src/QualityImprover/Relaxation
      ${MESQUITE_FOLDER}/src/QualityMetric
      ${MESQUITE_FOLDER}/src/QualityMetric/Debug
      ${MESQUITE_FOLDER}/src/QualityMetric/Shape
      ${MESQUITE_FOLDER}/src/QualityMetric/Smoothness
      ${MESQUITE_FOLDER}/src/QualityMetric/Volume
      ${MESQUITE_FOLDER}/src/QualityMetric/Untangle
      ${MESQUITE_FOLDER}/src/QualityMetric/TMP
      ${MESQUITE_FOLDER}/src/TargetCalculator
      ${MESQUITE_FOLDER}/src/TargetMetric
      ${MESQUITE_FOLDER}/src/TargetMetric/Misc
      ${MESQUITE_FOLDER}/src/TargetMetric/Shape
      ${MESQUITE_FOLDER}/src/TargetMetric/ShapeOrient
      ${MESQUITE_FOLDER}/src/TargetMetric/ShapeSize
      ${MESQUITE_FOLDER}/src/TargetMetric/ShapeSizeOrient
      ${MESQUITE_FOLDER}/src/TargetMetric/Size
      ${MESQUITE_FOLDER}/src/TargetMetric/Untangle
      ${MESQUITE_FOLDER}/src/Wrappers
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
    )
ENDIF()



MESSAGE(STATUS "Setting up bundled features - Done")

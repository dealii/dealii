## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Configuration for the Mesquite library:
#

#
# Mesquite writes some files into the build directory, so we need to include
# those as well.
#
MACRO(FEATURE_MESQUITE_CONFIGURE_BUNDLED)
  STRING(REGEX MATCH "mesquite-([0-9]+)\\.([0-9]+)([0-9]+)" MSQ_CONFIGURE_FOLDER  ${MESQUITE_FOLDER})

  SET(MESQUITE_BUNDLED_INCLUDE_DIRS
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
    )
ENDMACRO()

MACRO(FEATURE_MESQUITE_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find Mesquite and supporting libraries!\n"
    "Please ensure that the libraries are installed on your computer.\n"
    "If the libraries are not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ MESQUITE_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DMESQUITE_DIR=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    "Alternatively you may choose to compile the bundled libraries\n"
    "by setting DEAL_II_ALLOW_BUNDLED=ON or DEAL_II_FORCE_BUNDLED_MESQUITE=ON.\n\n"
    )
ENDMACRO()

CONFIGURE_FEATURE(MESQUITE)
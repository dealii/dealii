## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

if(DEAL_II_COMPONENT_EXAMPLES)
  message(STATUS "Setting up examples")

  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${DEAL_II_EXAMPLES_RELDIR}
    COMPONENT examples
    FILES_MATCHING
    #
    # Exclude folder structures: doc, doxygen, CMakeFiles,...
    #
    PATTERN "CMakeFiles*" EXCLUDE
    PATTERN "doc*" EXCLUDE
    #
    # Glob Includes:
    #
    PATTERN "*.cu"
    PATTERN "*.cc"
    PATTERN "*.prm"
    PATTERN "*.inp"
    PATTERN "*.ipynb"
    PATTERN "step*/CMakeLists.txt"
    #
    # Special files:
    #
    PATTERN "output.reference.dat"              # step-39
    PATTERN "postprocess.pl"                    # step-39
    PATTERN "obstacle.pbm"                      # step-42
    PATTERN "example.geo"                       # step-49
    PATTERN "example.msh"                       # step-49
    PATTERN "topography.txt.gz"                 # step-53
    PATTERN "input/initial_mesh_3d.vtk"         # step-54
    PATTERN "input/DTMB-5415_bulbous_bow.iges"  # step-54
    PATTERN "sphere_r6.geo"                     # step-97
    PATTERN "sphere_r6.msh"                     # step-97
    PATTERN "sphere_r7.geo"                     # step-97
    PATTERN "sphere_r7.msh"                     # step-97
    PATTERN "sphere_r8.geo"                     # step-97
    PATTERN "sphere_r8.msh"                     # step-97
    PATTERN "sphere_r9.geo"                     # step-97
    PATTERN "sphere_r9.msh"                     # step-97
    )

  message(STATUS "Setting up examples - Done")

endif()

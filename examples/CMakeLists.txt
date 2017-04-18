## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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

IF(DEAL_II_COMPONENT_EXAMPLES)
  MESSAGE(STATUS "Setting up examples")

  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
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
    PATTERN "*.cc"
    PATTERN "*.prm"
    PATTERN "*.inp"
    PATTERN "step*/CMakeLists.txt"
    PATTERN "doxygen/CMakeLists.txt"
    #
    # Special files:
    #
    PATTERN "output.reference.dat"              # step-39
    PATTERN "postprocess.pl"                    # step-39
    PATTERN "obstacle.pbm"                      # step-42
    PATTERN "untitled.geo"                      # step-49
    PATTERN "untitled.msh"                      # step-49
    PATTERN "topography.txt.gz"                 # step-53
    PATTERN "input/initial_mesh_3d.vtk"         # step-54
    PATTERN "input/DTMB-5415_bulbous_bow.iges"  # step-54
    )

  IF(DEAL_II_COMPILE_EXAMPLES)
    #
    # Set up all executables:
    #
    FILE(GLOB _steps ${CMAKE_CURRENT_SOURCE_DIR}/step-*/*.cc)
    FOREACH(_step ${_steps})
      GET_FILENAME_COMPONENT(_name ${_step} NAME_WE)
      GET_FILENAME_COMPONENT(_directory ${_step} DIRECTORY)

      #
      # Extract dependency information from CMakeLists.txt file.
      #
      SET(_dependency_string TRUE)
      IF(EXISTS "${_directory}/CMakeLists.txt")
        FILE(STRINGS "${_directory}/CMakeLists.txt" _dependency_string
          REGEX "^IF.*DEAL_II.* # keep in one line$"
          )
        STRING(REPLACE "IF(" "" _dependency_string "${_dependency_string}")
        STRING(REPLACE ") # keep in one line" "" _dependency_string "${_dependency_string}")
        IF("${_dependency_string}" STREQUAL "")
          SET(_dependency_string "TRUE")
        ELSE()
          SET(_dependency_string "NOT (${_dependency_string})")
        ENDIF()
      ENDIF()

      #
      # Now evaluate the _dependency_string and store the result in _setup.
      # We have to be careful - ${_dependency_string} evaluating to true
      # means that not all dependencies are fulfilled.
      #
      SET(_setup FALSE)
      EVALUATE_EXPRESSION("
        IF(${_dependency_string})
          SET(_setup TRUE)
        ENDIF()")

      IF(_setup)
        FOREACH(_build ${DEAL_II_BUILD_TYPES})
          STRING(TOLOWER ${_build} _build_lowercase)
          ADD_EXECUTABLE(${_name}.${_build_lowercase} ${_step})
          DEAL_II_INSOURCE_SETUP_TARGET(${_name}.${_build_lowercase} ${_build})
          SET_TARGET_PROPERTIES(${_name}.${_build_lowercase}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
            )
          ADD_DEPENDENCIES(examples ${_name}.${_build_lowercase})
          INSTALL(TARGETS ${_name}.${_build_lowercase}
            DESTINATION ${DEAL_II_EXAMPLES_RELDIR}/${_name}
            )
        ENDFOREACH()

      ELSE()

        MESSAGE(STATUS "  ${_name} - dependencies not satisfied")
      ENDIF()

    ENDFOREACH()
  ENDIF()

  MESSAGE(STATUS "Setting up examples - Done")

ENDIF()

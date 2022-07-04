## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2021 by the deal.II authors
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
    )

  IF(DEAL_II_COMPILE_EXAMPLES)
    #
    # Make sure that there are no deprecated functions used in the tutorials.
    #
    STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

    #
    # Set up all executables:
    #
    FILE(GLOB _steps
      ${CMAKE_CURRENT_SOURCE_DIR}/step-*/step-*.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/step-*/step-*.cu)
    FOREACH(_step ${_steps})
      GET_FILENAME_COMPONENT(_name ${_step} NAME_WE)
      GET_FILENAME_COMPONENT(_directory ${_step} DIRECTORY)

      #
      # Extract dependency information from CMakeLists.txt file.
      #
      SET(_setup FALSE)
      IF(EXISTS "${_directory}/CMakeLists.txt")
        FILE(STRINGS "${_directory}/CMakeLists.txt" _dependency_string
          REGEX "^IF.*DEAL_II.* # keep in one line$"
          )
        STRING(REPLACE "IF(" "" _dependency_string "${_dependency_string}")
        STRING(REPLACE ") # keep in one line" "" _dependency_string "${_dependency_string}")
        IF("${_dependency_string}" STREQUAL "")
          SET(_setup TRUE)
        ELSE()
          # if the dependency string evaluates to TRUE then the example
          # CMakeLists.txt encounters a fatal error - we want the opposite logic
          # here so add a NOT.
          EVALUATE_EXPRESSION("
            IF(NOT (${_dependency_string}))
              SET(_setup TRUE)
            ENDIF()")
        ENDIF()
      ENDIF()

      IF(_setup)
        FOREACH(_build ${DEAL_II_BUILD_TYPES})
          STRING(TOLOWER ${_build} _build_lowercase)
          ADD_EXECUTABLE(${_name}.${_build_lowercase} ${_step})
          DEAL_II_INSOURCE_SETUP_TARGET(${_name}.${_build_lowercase} ${_build})
          SET_TARGET_PROPERTIES(${_name}.${_build_lowercase}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
            )
          #
          # In case CMake is instructed to add rpaths to the library and
          # exectuble on installation, make sure that we add an additional
          # rpath to the library location as well:
          #
          IF(CMAKE_INSTALL_RPATH_USE_LINK_PATH)
            SET_TARGET_PROPERTIES(${_name}.${_build_lowercase}
              PROPERTIES
              INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}"
              )
          ENDIF()

          ADD_DEPENDENCIES(examples ${_name}.${_build_lowercase})
          INSTALL(TARGETS ${_name}.${_build_lowercase}
            DESTINATION ${DEAL_II_EXAMPLES_RELDIR}/${_name}
            )
        ENDFOREACH()

      ELSE()

        MESSAGE(STATUS "  ${_name} - dependencies not satisfied")
      ENDIF()

    ENDFOREACH()

    # the same as above but for the examples folder
    FILE(GLOB _steps
      ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/*.cc)
    FOREACH(_step ${_steps})
      GET_FILENAME_COMPONENT(_name ${_step} NAME_WE)

        FOREACH(_build ${DEAL_II_BUILD_TYPES})
          STRING(TOLOWER ${_build} _build_lowercase)
          ADD_EXECUTABLE(${_name}.${_build_lowercase} ${_step})
          DEAL_II_INSOURCE_SETUP_TARGET(${_name}.${_build_lowercase} ${_build})

          SET_TARGET_PROPERTIES(${_name}.${_build_lowercase}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
            )
        ENDFOREACH()
    ENDFOREACH()
  ENDIF()

  MESSAGE(STATUS "Setting up examples - Done")

ENDIF()

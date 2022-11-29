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
    )

  if(DEAL_II_COMPILE_EXAMPLES)
    #
    # Make sure that there are no deprecated functions used in the tutorials.
    #
    strip_flag(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

    #
    # Set up all executables:
    #
    file(GLOB _steps
      ${CMAKE_CURRENT_SOURCE_DIR}/step-*/step-*.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/step-*/step-*.cu)
    foreach(_step ${_steps})
      get_filename_component(_name ${_step} NAME_WE)
      get_filename_component(_directory ${_step} DIRECTORY)

      #
      # Extract dependency information from CMakeLists.txt file.
      #
      set(_setup FALSE)
      if(EXISTS "${_directory}/CMakeLists.txt")
        file(STRINGS "${_directory}/CMakeLists.txt" _dependency_string
          REGEX "^if.*DEAL_II.* # keep in one line$"
          )
        string(REPLACE "if(" "" _dependency_string "${_dependency_string}")
        string(REPLACE ") # keep in one line" "" _dependency_string "${_dependency_string}")
        if("${_dependency_string}" STREQUAL "")
          set(_setup TRUE)
        else()
          # if the dependency string evaluates to TRUE then the example
          # CMakeLists.txt encounters a fatal error - we want the opposite logic
          # here so add a NOT.
          evaluate_expression("
            if(NOT (${_dependency_string}))
              set(_setup TRUE)
            endif()")
        endif()
      endif()

      if(_setup)
        foreach(_build ${DEAL_II_BUILD_TYPES})
          string(TOLOWER ${_build} _build_lowercase)
          set(_target "example_${_name}_${_build_lowercase}")
          string(REPLACE "-" "_" _target "${_target}")
          add_executable(${_target} ${_step})
          insource_setup_target(${_target} ${_build})
          set_target_properties(${_target}
            PROPERTIES
            RUNTIME_OUTPUT_NAME "${_name}.${_build_lowercase}"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
            )
          #
          # In case CMake is instructed to add rpaths to the library and
          # exectuble on installation, make sure that we add an additional
          # rpath to the library location as well:
          #
          if(CMAKE_INSTALL_RPATH_USE_LINK_PATH)
            set_target_properties(${_target}
              PROPERTIES
              INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}"
              )
          endif()

          add_dependencies(examples ${_target})
          install(TARGETS ${_target}
            DESTINATION ${DEAL_II_EXAMPLES_RELDIR}/${_name}
            )
        endforeach()

      else()

        message(STATUS "  ${_name} - dependencies not satisfied")
      endif()

    endforeach()

    #
    # Also compile some documentation examples:
    #

    file(GLOB _steps ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/*.cc)
    foreach(_step ${_steps})
      get_filename_component(_name ${_step} NAME_WE)

        foreach(_build ${DEAL_II_BUILD_TYPES})
          string(TOLOWER ${_build} _build_lowercase)
          set(_target "example_${_name}_${_build_lowercase}")
          string(REPLACE "-" "_" _target "${_target}")
          add_executable(${_target} ${_step})
          insource_setup_target(${_target} ${_build})
          set_target_properties(${_target}
            PROPERTIES
            RUNTIME_OUTPUT_NAME "${_name}.${_build_lowercase}"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
            )
        endforeach()
    endforeach()
  endif()

  message(STATUS "Setting up examples - Done")

endif()

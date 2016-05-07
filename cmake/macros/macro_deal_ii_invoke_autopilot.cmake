## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
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
# This file implements the DEAL_II_INVOKE_AUTOPILOT macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_INVOKE_AUTOPILOT()
#
# where it is assumed that the following variables are defined:
#
#       TARGET         -  a string used for the project and target name
#       TARGET_SRC     -  a list of source file to compile for target
#                         ${TARGET}
#       TARGET_RUN     -  (optional) the command line that should be
#                         invoked by "make run", will be set to default
#                         values if undefined. If no run target should be
#                         created, set it to an empty string.
#       CLEAN_UP_FILES -  (optional) a list of files (globs) that will be
#                         removed with "make runclean" and "make
#                         distclean", will be set to default values if
#                         empty
#

MACRO(DEAL_II_INVOKE_AUTOPILOT)

  # Set CMAKE_BUILD_TYPE=Debug if both 
  # Debug and Release mode are given
  IF("${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      FORCE)
  ENDIF()


  # Generator specific values:
  IF(CMAKE_GENERATOR MATCHES "Ninja")
    SET(_make_command "$ ninja")
  ELSE()
    SET(_make_command " $ make")
  ENDIF()

  # Define and setup a compilation target:
  ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})
  DEAL_II_SETUP_TARGET(${TARGET})

  MESSAGE(STATUS "Autopilot invoked")

  # Define a custom target to easily run the program:

  IF(NOT DEFINED TARGET_RUN)
    SET(TARGET_RUN ${TARGET})
  ENDIF()

  IF(CMAKE_SYSTEM_NAME MATCHES "(CYGWIN|Windows)")
    #
    # Hack for Cygwin and Windows targets: Export PATH to point to the
    # dynamic library.
    #
    SET(_delim ":")
    IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
      SET(_delim ";")
    ENDIF()
    FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/run_target.cmake
      "SET(ENV{PATH} \"${CMAKE_CURRENT_BINARY_DIR}${_delim}${DEAL_II_PATH}/${DEAL_II_EXECUTABLE_RELDIR}${_delim}\$ENV{PATH}\")\n"
      "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TYPE}\\\\${TARGET_RUN}\n"
      "  RESULT_VARIABLE _return_value\n"
      "  )\n"
      "IF(NOT \"\${_return_value}\" STREQUAL \"0\")\n"
      "  MESSAGE(SEND_ERROR \"\nProgram terminated with exit code: \${_return_value}\")\n"
      "ENDIF()\n"
      )
    SET(_command
      ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/run_target.cmake
      )

  ELSE()

    SET(_command ${TARGET_RUN})
  ENDIF()

  IF(NOT "${TARGET_RUN}" STREQUAL "")
    ADD_CUSTOM_TARGET(run
      COMMAND ${_command}
      DEPENDS ${TARGET}
      COMMENT "Run ${TARGET} with ${CMAKE_BUILD_TYPE} configuration"
      )
    SET(_run_targets
      "#      ${_make_command} run            - to (compile, link and) run the program\n"
      )
  ENDIF()


  #
  # Provide a target to sign the generated executable with a Mac OSX
  # developer key. This avoids problems with an enabled firewall and MPI
  # tasks that need networking.
  #

  IF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    IF(DEFINED OSX_CERTIFICATE_NAME)
      ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${TARGET}.signed
        COMMAND codesign -f -s ${OSX_CERTIFICATE_NAME} ${TARGET}
        COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${TARGET}.signed
        COMMENT "Digitally signing ${TARGET}"
        DEPENDS ${TARGET}
        VERBATIM
        )
      ADD_CUSTOM_TARGET(sign ALL
        DEPENDS ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${TARGET}.signed
        )
      ADD_DEPENDENCIES(run sign)
    ELSE()
      ADD_CUSTOM_TARGET(sign
        COMMAND
           ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo '**           Error: No Mac OSX developer certificate specified           **'
        && ${CMAKE_COMMAND} -E echo '**         Please reconfigure with -DOSX_CERTIFICATE_NAME="<...>"        **'
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo ''
        COMMENT "Digitally signing ${TARGET}"
        )
    ENDIF()

    SET(_run_targets
      "${_run_targets}#\n#      ${_make_command} sign           - to sign the executable with the supplied OSX developer key\n"
      )
  ENDIF()

  # Define custom targets to easily switch the build type:
  ADD_CUSTOM_TARGET(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
    )

  ADD_CUSTOM_TARGET(release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Release"
    )

  # Only mention release and debug targets if it is actually possible to
  # switch between them:
  IF(${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
    SET(_switch_targets
"#      ${_make_command} debug          - to switch the build type to 'Debug'
#      ${_make_command} release        - to switch the build type to 'Release'\n"
      )
  ENDIF()

  # And another custom target to clean up all files generated by the program:
  IF("${CLEAN_UP_FILES}" STREQUAL "")
    SET(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd *.d2)
  ENDIF()
  ADD_CUSTOM_TARGET(runclean
    COMMAND ${CMAKE_COMMAND} -E remove ${CLEAN_UP_FILES}
    COMMENT "runclean invoked"
    )

  # Define a distclean target to remove every generated file:
  ADD_CUSTOM_TARGET(distclean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target runclean
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
    COMMENT "distclean invoked"
    )

  # Define a strip-comments target:
  FIND_PACKAGE(Perl QUIET)
  IF(PERL_FOUND)
    ADD_CUSTOM_TARGET(strip_comments
      COMMAND ${PERL_EXECUTABLE} -pi -e 's\#^[ \\t]*//.*\\n\#\#g;' ${TARGET_SRC}
      COMMENT "strip comments"
      )
  ENDIF()


  # Print out some usage information to file:
  FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"MESSAGE(
\"###
#
#  Project  ${TARGET}  set up with  ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION}  found at
#      ${DEAL_II_PATH}
#
#  CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#  You can now run
#      ${_make_command}                - to compile and link the program
${_run_targets}#
${_switch_targets}#
")
  IF(NOT CMAKE_GENERATOR MATCHES "Ninja")
    FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"#      ${_make_command} edit_cache     - to change (cached) configuration variables
#                               and rerun the configure and generate phases of CMake
#
")
  ENDIF()
  IF(PERL_FOUND)
    FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"#      ${_make_command} strip_comments - to strip the source files in this
#                               directory off the documentation comments
")
  ENDIF()
  FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"#      ${_make_command} clean          - to remove the generated executable as well as
#                               all intermediate compilation files
#      ${_make_command} runclean       - to remove all output generated by the program
#      ${_make_command} distclean      - to clean the directory from _all_ generated
#                               files (includes clean, runclean and the removal
#                               of the generated build system)
#      ${_make_command} info           - to view this message again
#
#  Have a nice day!
#
###\")"
     )

  ADD_CUSTOM_TARGET(info
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
    )

  # Print this message once:
  IF(NOT USAGE_PRINTED)
    INCLUDE(${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake)
    SET(USAGE_PRINTED TRUE CACHE INTERNAL "")
  ELSE()
    MESSAGE(STATUS "Run  ${_make_command} info  to print a detailed help message")
  ENDIF()

ENDMACRO()


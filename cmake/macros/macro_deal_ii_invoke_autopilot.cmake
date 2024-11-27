## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This file implements the DEAL_II_INVOKE_AUTOPILOT macro, which is
# part of the deal.II library.
#
# Usage:
#       deal_ii_invoke_autopilot()
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

macro(deal_ii_invoke_autopilot)
  message(STATUS "Autopilot invoked")

  get_property(_generator_is_multi_config
    GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG
    )

  if(${_generator_is_multi_config})
    message(STATUS "Setting up multiple generator configurations")
    set(CMAKE_BUILD_TYPE "multiple generator")
    set(CMAKE_CONFIGURATION_TYPES "")
    if(${DEAL_II_BUILD_TYPE} MATCHES "Debug")
      list(APPEND CMAKE_CONFIGURATION_TYPES "Debug")
    endif()
    if(${DEAL_II_BUILD_TYPE} MATCHES "Release")
      list(APPEND CMAKE_CONFIGURATION_TYPES "Release")
    endif()
  endif()

  #
  # Generate compile_commands.json
  #

  set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

  #
  # Make sure that compile_commands.json is available in the source
  # directory (to help IDEs like VSCode with clangd for code analysis):
  #

  # In order to check whether we have an "in source" build we have to make
  # sure to remove all symlinks and resolve to absolute file paths. This is
  # necessary because with symlinks CMAKE_SOURCE_DIR and CMAKE_BINARY dir
  # might look different but are in fact the same (fully resolved)
  # directory...
  get_filename_component(_source "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(_build  "${CMAKE_BINARY_DIR}" REALPATH)
  if(NOT "${_source}" STREQUAL "${_build}")
    # Ensure that we do not override "compile_commands.json" in the source
    # directory if already present and then create a symlink using the
    # (unresolved) CMAKE_BINARY_DIR variable.
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/compile_commands.json")
      file(CREATE_LINK
        "${CMAKE_BINARY_DIR}/compile_commands.json"
        "${CMAKE_SOURCE_DIR}/compile_commands.json"
        SYMBOLIC
        )
    endif()
  endif()

  # Define and set up a compilation target:
  add_executable(${TARGET} ${TARGET_SRC})

  #
  # To ensure maximal compatibility with existing user codes we use the
  # deal_ii_setup_target() macro when setting up the target for a standard
  # single-configuration generator und switch to the dealii::dealii target
  # only for a multiple generator configuration.
  #
  # This hopefully ensures that subtle differences in how cmake expands
  # compiler and linker flags coming from either the target directly or
  # from the interface target dealii::dealii do not propagate to user
  # project configurations.
  #
  if(${_generator_is_multi_config})
    if(NOT DEAL_II_TARGET_CONFIG_INCLUDED)
      include(${DEAL_II_TARGET_CONFIG})
      set(DEAL_II_TARGET_CONFIG_INCLUDED TRUE)
    endif()
    target_link_libraries(${TARGET} dealii::dealii)
  else()
    deal_ii_setup_target(${TARGET})
  endif()

  # Generator specific values:
  if(CMAKE_GENERATOR MATCHES "Ninja")
    set(_make_command "$ ninja")
  else()
    set(_make_command " $ make")
  endif()

  # Define a custom target to easily run the program:
  if(NOT DEFINED TARGET_RUN)
    set(TARGET_RUN ${TARGET})
  endif()

  if(CMAKE_SYSTEM_NAME MATCHES "(CYGWIN|Windows)")
    #
    # Hack for Cygwin and Windows targets: Export PATH to point to the
    # dynamic library.
    #
    set(_delim ":")
    if(CMAKE_SYSTEM_NAME MATCHES "Windows")
      set(_delim ";")
    endif()
    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/run_target.cmake
      "set(ENV{PATH} \"${CMAKE_CURRENT_BINARY_DIR}${_delim}${DEAL_II_PATH}/${DEAL_II_EXECUTABLE_RELDIR}${_delim}\$ENV{PATH}\")\n"
      "execute_process(COMMAND ${CMAKE_BUILD_TYPE}\\\\${TARGET_RUN}\n"
      "  RESULT_VARIABLE _return_value\n"
      "  )\n"
      "if(NOT \"\${_return_value}\" STREQUAL \"0\")\n"
      "  message(SEND_ERROR \"\nProgram terminated with exit code: \${_return_value}\")\n"
      "endif()\n"
      )
    set(_command
      ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/run_target.cmake
      )
  else()
    set(_command ${TARGET_RUN})
  endif()

  if(NOT "${TARGET_RUN}" STREQUAL "")
    add_custom_target(run
      COMMAND ${_command}
      DEPENDS ${TARGET}
      COMMENT "Run ${TARGET} with $<$<CONFIG:Release>:Release>$<$<CONFIG:Debug>:Debug> configuration"
      )
    set(_run_targets
      "#      ${_make_command} run            - to (compile, link and) run the program\n"
      )
  endif()

  #
  # Provide a target to sign the generated executable with a Mac OSX
  # developer key. This avoids problems with an enabled firewall and MPI
  # tasks that need networking.
  #

  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    if(DEFINED OSX_CERTIFICATE_NAME)
      add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${TARGET}.signed
        COMMAND codesign -f -s ${OSX_CERTIFICATE_NAME} ${TARGET}
        COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${TARGET}.signed
        COMMENT "Digitally signing ${TARGET}"
        DEPENDS ${TARGET}
        VERBATIM
        )
      add_custom_target(sign ALL
        DEPENDS ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${TARGET}.signed
        )
      add_dependencies(run sign)
    else()
      add_custom_target(sign
        COMMAND
           ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo '**           Error: No Mac OSX developer certificate specified           **'
        && ${CMAKE_COMMAND} -E echo '**         Please reconfigure with -DOSX_CERTIFICATE_NAME="<...>"        **'
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo ''
        COMMENT "Digitally signing ${TARGET}"
        )
    endif()

    set(_run_targets
      "${_run_targets}#\n#      ${_make_command} sign           - to sign the executable with the supplied OSX developer key\n"
      )
  endif()

  if(NOT ${_generator_is_multi_config})
    #
    # Define custom targets to easily switch the build type:
    #

    if(${DEAL_II_BUILD_TYPE} MATCHES "Debug")
      add_custom_target(debug
        COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug mode. Now recompile with: ${_make_command}"
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMENT "Switching CMAKE_BUILD_TYPE to Debug"
        VERBATIM
        )
    endif()

    if(${DEAL_II_BUILD_TYPE} MATCHES "Release")
      add_custom_target(release
        COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Release mode. Now recompile with: ${_make_command}"
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMENT "Switching CMAKE_BUILD_TYPE to Release"
        VERBATIM
        )
    endif()
  endif()

  #
  # And another custom target to clean up all files generated by the program:
  #

  if("${CLEAN_UP_FILES}" STREQUAL "")
    set(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd *.d2)
  endif()
  add_custom_target(runclean
    COMMAND ${CMAKE_COMMAND} -E remove ${CLEAN_UP_FILES}
    COMMENT "runclean invoked"
    )

  #
  # Define a distclean target to remove every generated file:
  #

  add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target runclean
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove
      CMakeCache.txt cmake_install.cmake compile_commands.json Makefile
      build.ninja rules.ninja .ninja_deps .ninja_log
    COMMENT "distclean invoked"
    )

  #
  # Define a strip_comments target:
  #

  find_package(Perl QUIET)
  if(PERL_FOUND)
    add_custom_target(strip_comments
      COMMAND ${PERL_EXECUTABLE} -pi -e 's\#^[ \\t]*//.*\\n\#\#g;' ${CMAKE_SOURCE_DIR}/${TARGET_SRC}
      COMMENT "strip comments"
      )
  endif()

  #
  # Construct help message:
  #

  # Only mention release and debug targets if it is actually possible to
  # switch between them:
  if(NOT ${_generator_is_multi_config} AND ${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
    set(_switch_targets
"#      ${_make_command} debug          - to switch the build type to 'Debug'
#      ${_make_command} release        - to switch the build type to 'Release'\n"
      )
  endif()

  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"message(
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
  if(NOT CMAKE_GENERATOR MATCHES "Ninja")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"#      ${_make_command} edit_cache     - to change (cached) configuration variables
#                               and rerun the configure and generate phases of CMake
#
")
  endif()
  if(PERL_FOUND)
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"#      ${_make_command} strip_comments - to strip the source files in this
#                               directory off their comments; this is irreversible
")
  endif()
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
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

  add_custom_target(info
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
    )

  # Print this message once:
  if(NOT USAGE_PRINTED)
    include(${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake)
    set(USAGE_PRINTED TRUE CACHE INTERNAL "")
  else()
    message(STATUS "Run  ${_make_command} info  to print a detailed help message")
  endif()

endmacro()

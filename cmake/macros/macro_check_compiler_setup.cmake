## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2016 - 2023 by the deal.II authors
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
# Usage:
#   check_compiler_setup("compiler flag string" "linker flag string" _var
#     [libraries]
#     )
#
# This macro tries to compile and link a simple "int main(){ return 0; }
# with the given set of provided compiler and linker flags  and an optional
# list of libraries to link against. If the test is successful the variable
# ${_var} is set to true, otherwise it is set to false.
#

macro(check_compiler_setup _compiler_flags _linker_flags _var)
  message(STATUS "Performing Test ${_var}")

  #
  # We used to do some fancy caching here, but for now simply rerun this
  # check every time and do not bother with an elaborate caching strategy.
  # We simply keep the build directory for the check around and rely on the
  # build system to skip unnecessary recompilation.
  #

  #
  # Prepare compile and link options:
  #

  separate_arguments(_compile_options UNIX_COMMAND "${_compiler_flags}")
  shell_escape_option_groups(_compile_options)
  separate_arguments(_link_options UNIX_COMMAND "${_linker_flags}")
  shell_escape_option_groups(_link_options)

  #
  # Ideally, we would like to simply list our internal interface_* targets,
  # which might recursively list imported targets... But unfortunately,
  # CMake as of version 3.25 does not support this operation.
  #
  # Remark: What we try to accomplish here is somewhat fundamentally
  # incompatible with how CMake envisions import targets to function.
  # Namely, the final link line is only constructed after the "configure
  # phase" during the "generator phase".
  #
  # As a workaround, at least for now, let us make the assumption that we
  # only ever encounter an interface_* target and that imported targets
  # have been expanded.
  #

  set(_libraries)
  foreach(_entry ${ARGN})
    if(TARGET ${_entry})
      get_target_property(_value ${_entry} INTERFACE_LINK_LIBRARIES)
      if(NOT "${_value}" MATCHES "-NOTFOUND")
        list(APPEND _libraries ${_value})
      endif()
      get_target_property(_values ${_entry} INTERFACE_COMPILE_OPTIONS)
      if(NOT "${_values}" MATCHES "-NOTFOUND")
        list(APPEND _compile_options ${_values})
      endif()
      get_target_property(_values ${_entry} INTERFACE_LINK_OPTIONS)
      if(NOT "${_values}" MATCHES "-NOTFOUND")
        list(APPEND _link_options ${_values})
      endif()
    else()
      list(APPEND _libraries ${_entry})
    endif()
  endforeach()

  try_compile(
    ${_var}
    ${CMAKE_CURRENT_BINARY_DIR}/check_compiler_setup/${_var}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros/check_compiler_setup
    CheckCompilerSetup${_var}
    CMAKE_FLAGS "-DTEST_COMPILE_OPTIONS=${_compile_options}"
                "-DTEST_LINK_OPTIONS=${_link_options}"
                "-DTEST_LINK_LIBRARIES=${_libraries}"
                "-DCMAKE_VERBOSE_MAKEFILE=ON"
    OUTPUT_VARIABLE _output
    )

  if(${_var})
    message(STATUS "Performing Test ${_var} - Success")
  else()
    message(STATUS "Performing Test ${_var} - Failed")
    message(STATUS "Compiler output:\n\n${_output}")
  endif()
endmacro()

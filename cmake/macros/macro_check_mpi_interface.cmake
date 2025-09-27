## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2024 by the deal.II authors
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
# Check whether a feature is compiled against the same MPI library as the
# one deal.II picked up
#
# Usage:
#     check_mpi_interface(_feature _var),
#

macro(check_mpi_interface _feature _var)
  if(DEAL_II_WITH_MPI AND MPI_LIBRARIES)

    set(_nope FALSE)

    foreach(_library ${${_feature}_LIBRARIES})
      if( _library MATCHES "/libmpi\\.(a|so)[^/]*$")

        get_filename_component(_file1 ${_library} REALPATH)

        set(_not_found TRUE)
        foreach(_mpi_library ${MPI_LIBRARIES})
          get_filename_component(_file2 ${_mpi_library} REALPATH)
          if("${_file1}" STREQUAL "${_file2}")
            set(_not_found FALSE)
            break()
          endif()
        endforeach()

        if(_not_found)
          set(_nope TRUE)
          set(_spurious_library ${_library})
          break()
        endif()
      endif()
    endforeach()

    if(_nope)
      message(STATUS "Could not find a sufficient ${_feature} installation: "
        "${_feature} is compiled against a different MPI library than the one "
        "deal.II picked up."
        )
      to_string(_str ${MPI_LIBRARIES})
      set(${_feature}_ADDITIONAL_ERROR_STRING
        ${${_feature}_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ${_feature} installation:\n"
        "${_feature} has to be compiled against the same MPI library as deal.II "
        "but the link line of ${_feature} contains:\n"
        "  ${_spurious_library}\n"
        "which is not listed in MPI_LIBRARIES:\n"
        "  MPI_LIBRARIES = \"${_str}\"\n"
        )
      set(${_var} FALSE)
    endif()
  endif()
endmacro()

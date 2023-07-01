## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# define_object_target(<library>)
#
# A small wrapper around add_library that will define a target for each
# build type specified in DEAL_II_BUILD_TYPES. The specified library name
# must start with either object_ or bundled_.
#
# The function automatically "links" all object targets to the interface
# targets specified in DEAL_II_TARGETS and DEAL_II_TARGETS_(BUILD|RELEASE).
# The scope is PUBLIC so that properties of the interface targets propagate
# to the final shared library targets.
#
# In addition, if the library name is of the form "object_*" then all
# bundled targets are added to the "link" interface as well. The scope is
# private.
#

function(define_object_library _library)

  if(NOT "${_library}" MATCHES "^(object|bundled)_")
    message(FATAL_ERROR
      "Internal error: The specified target name must begin with object_ "
      "or bundled_. Encountered: ${_library}"
      )
  endif()

  foreach(_build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)
    set(_target "${_library}_${_build_lowercase}")

    add_library(${_target} ${ARGN})

    #
    # Add standard properties defined in DEAL_II_* variables:
    #

    populate_target_properties(${_target} ${_build})

    if("${_library}" MATCHES "^bundled_")
      #
      # Record all bundled object libraries in the global property
      # DEAL_II_BUNDLED_TARGETS_${_build}
      #
      set_property(GLOBAL APPEND PROPERTY DEAL_II_BUNDLED_TARGETS_${_build} ${_target})
    else()
      get_property(_bundled_object_targets
        GLOBAL PROPERTY DEAL_II_BUNDLED_TARGETS_${build}
        )
      target_link_libraries(${_target} PRIVATE ${_bundled_object_targets})
    endif()

    set_property(GLOBAL APPEND PROPERTY DEAL_II_OBJECT_TARGETS_${_build} ${_target})
  endforeach()

endfunction()

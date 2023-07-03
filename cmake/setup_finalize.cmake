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


########################################################################
#                                                                      #
#                      Finalize the configuration:                     #
#                                                                      #
########################################################################

#
# Hide some cmake specific cached variables. This is annoying...
#
mark_as_advanced(file_cmd)

#
# Append the saved initial (cached) variables ${flags}_SAVED at the end of
# ${flags}, see setup_cached_compiler_flags.cmake and the main
# CMakeLists.txt for details.
#
foreach(_flags ${DEAL_II_USED_FLAGS})
  # Strip leading and trailing whitespace:
  string(STRIP "${${_flags}} ${${_flags}_SAVED}" ${_flags})
endforeach()

#
# Sanity check: The variables defined in DEAL_II_REMOVED_FLAGS must not be
# used during the configuration stage:
#
foreach(_flag ${DEAL_II_REMOVED_FLAGS})
  if(NOT "${${_flag}}" STREQUAL "")
    message(FATAL_ERROR
      "\nInternal configuration error: The variable ${_flag} was set to a "
      "non empty value during the configuration! (The corresponding "
      "DEAL_II_* variable should have been used.)\n"
      "${_flag}=\"${${_flag}}\"\n"
      )
  endif()
endforeach()

########################################################################
#                                                                      #
#                Query for git repository information:                 #
#                                                                      #
########################################################################

deal_ii_query_git_information("DEAL_II")

file(WRITE ${CMAKE_BINARY_DIR}/revision.log
"###
#
#  Git information:
#        Branch:    ${DEAL_II_GIT_BRANCH}
#        Revision:  ${DEAL_II_GIT_REVISION}
#        Timestamp: ${DEAL_II_GIT_TIMESTAMP}
#
###"
  )


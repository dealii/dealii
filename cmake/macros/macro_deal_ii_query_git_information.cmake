## ---------------------------------------------------------------------
##
## Copyright (C) 2015 by the deal.II authors
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
# This file implements the DEAL_II_QUERY_GIT_INFORMATION macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_QUERY_GIT_INFORMATION()
#
# This will try to gather information about current branch, as well as
# short and long revision. If ${CMAKE_SOURCE_DIR} is the root of a git
# repository the following variables will be populated:
#
#       GIT_BRANCH
#       GIT_REVISION
#       GIT_SHORTREV
#
# If this macro is called within the deal.II build system the variables are
# prefixed with DEAL_II_:
#
#       DEAL_II_GIT_BRANCH
#       DEAL_II_GIT_REVISION
#       DEAL_II_GIT_SHORTREV
#

MACRO(DEAL_II_QUERY_GIT_INFORMATION)

  MESSAGE(STATUS "Query git repository information.")

  #
  # If DEAL_II_BASE_NAME is defined and DEAL_II_PROJECT_CONFIG_INCLUDED was
  # not set, we assume that we are called from within the deal.II build
  # system. In this case we prepend all variables by "DEAL_II_"
  #
  IF( DEFINED DEAL_II_BASE_NAME AND
      NOT DEFINED DEAL_II_PROJECT_CONFIG_INCLUDED )
    SET(_prefix "DEAL_II_")
  ELSE()
    SET(_prefix "")
  ENDIF()

  FIND_PACKAGE(Git)

  #
  # Only run the following if we have git and the source directory seems to
  # be under version control.
  #
  IF(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git/HEAD)
    #
    # Bogus configure_file calls to trigger a reconfigure, and thus an
    # update of branch and commit information every time HEAD has changed.
    #
    CONFIGURE_FILE(
      ${CMAKE_SOURCE_DIR}/.git/HEAD
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/HEAD
      )
    FILE(STRINGS ${CMAKE_SOURCE_DIR}/.git/HEAD _head_ref LIMIT_COUNT 1)
    STRING(REPLACE "ref: " "" _head_ref ${_head_ref})
    IF(EXISTS ${CMAKE_SOURCE_DIR}/.git/${_head_ref})
      CONFIGURE_FILE(
        ${CMAKE_SOURCE_DIR}/.git/${_head_ref}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/HEAD_REF
        )
    ENDIF()

    #
    # Query for revision:
    #

    EXECUTE_PROCESS(
       COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:"%H %h"
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
       OUTPUT_VARIABLE _info
       RESULT_VARIABLE _result
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
    IF(${_result} EQUAL 0)
      STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+)\"$"
        "\\1" ${_prefix}GIT_REVISION "${_info}")
      STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+)\"$"
        "\\2" ${_prefix}GIT_SHORTREV "${_info}")
    ENDIF()

    #
    # Query for branch:
    #

    EXECUTE_PROCESS(
       COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
       OUTPUT_VARIABLE _branch
       RESULT_VARIABLE _result
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
    IF(${_result} EQUAL 0)
      STRING(REGEX REPLACE "refs/heads/" "" ${_prefix}GIT_BRANCH "${_branch}")
    ENDIF()
  ENDIF()

ENDMACRO()

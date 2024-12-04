## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2015 - 2022 by the deal.II authors
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
# This file implements the DEAL_II_QUERY_GIT_INFORMATION macro, which is
# part of the deal.II library.
#
# Usage:
#       deal_ii_query_git_information()
#       deal_ii_query_git_information("CUSTOM_PREFIX")
#
# This will try to gather information about current branch, as well as
# short and long revision. If ${DEAL_II_SOURCE_DIR} is the root of a git
# repository the following variables will be populated:
#
#       GIT_BRANCH
#       GIT_REVISION
#       GIT_SHORTREV
#       GIT_TAG
#
# The macro can be called with an optional PREFIX argument to prefix the
# variables:
#
#       PREFIX_GIT_BRANCH
#       PREFIX_GIT_REVISION
#       PREFIX_GIT_SHORTREV
#       PREFIX_GIT_TAG
#       PREFIX_GIT_TIMESTAMP
#

macro(deal_ii_query_git_information)

  message(STATUS "Query git repository information.")

  # Set prefix.
  set(_prefix "")
  if(NOT "${ARGN}" STREQUAL "")
    set(_prefix "${ARGN}_")
  endif()

  find_package(Git)

  #
  # Only run the following if we have git and the source directory seems to
  # be under version control.
  #
  if(GIT_FOUND AND EXISTS ${DEAL_II_SOURCE_DIR}/.git/HEAD)
    #
    # Bogus configure_file calls to trigger a reconfigure, and thus an
    # update of branch and commit information every time HEAD has changed.
    #
    configure_file(
      ${DEAL_II_SOURCE_DIR}/.git/HEAD
      ${DEAL_II_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/HEAD
      )
    file(STRINGS ${DEAL_II_SOURCE_DIR}/.git/HEAD _head_ref LIMIT_COUNT 1)
    string(REPLACE "ref: " "" _head_ref ${_head_ref})
    if(EXISTS ${DEAL_II_SOURCE_DIR}/.git/${_head_ref})
      configure_file(
        ${DEAL_II_SOURCE_DIR}/.git/${_head_ref}
        ${DEAL_II_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/HEAD_REF
        )
    endif()

    #
    # Query for revision:
    #
    # Format options for git log:
    #   %H   - the full sha1 hash of the commit, aka "revision"
    #   %h   - the abbreviated sha1 hash of the commit, aka "shortrev"
    #   %cd  - the commit date (for which we use the strict iso date format)
    #

    # --date=iso-strict has been introduced in git 2.2 (released Dec 2014)
    set(_date_format)
    if(NOT ${GIT_VERSION_STRING} VERSION_LESS 2.2)
      set(_date_format "--date=iso-strict")
    endif()

    execute_process(
       COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:"revision=%H, shortrev=%h, date=%cd" ${_date_format}
       WORKING_DIRECTORY ${DEAL_II_SOURCE_DIR}
       OUTPUT_VARIABLE _info
       RESULT_VARIABLE _result
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
    if(${_result} EQUAL 0)
      string(REGEX REPLACE "^\"revision=(.+), shortrev=(.+), date=(.+)\"$"
        "\\1" ${_prefix}GIT_REVISION "${_info}")
      string(REGEX REPLACE "^\"revision=(.+), shortrev=(.+), date=(.+)\"$"
        "\\2" ${_prefix}GIT_SHORTREV "${_info}")
      string(REGEX REPLACE "^\"revision=(.+), shortrev=(.+), date=(.+)\"$"
        "\\3" ${_prefix}GIT_TIMESTAMP "${_info}")

      #
      # Replace "T" by a space in order to have the same output as
      #   $ date --rfc-3339=seconds"
      #
      string(REPLACE "T" " " ${_prefix}GIT_TIMESTAMP "${${_prefix}GIT_TIMESTAMP}")
    endif()

    #
    # Query for branch:
    #

    execute_process(
       COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
       WORKING_DIRECTORY ${DEAL_II_SOURCE_DIR}
       OUTPUT_VARIABLE _branch
       RESULT_VARIABLE _result
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
    if(${_result} EQUAL 0)
      string(REGEX REPLACE "refs/heads/" "" ${_prefix}GIT_BRANCH "${_branch}")
    endif()

    #
    # Query for tag:
    #
    set(_script "")
    if(EXISTS     ${DEAL_II_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/scripts/get_latest_tag.sh)
      set(_script ${DEAL_II_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/scripts/get_latest_tag.sh)
    elseif(EXISTS ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/get_latest_tag.sh)
      set(_script ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/get_latest_tag.sh)
    endif()
    if(NOT "${_script}" STREQUAL "")
       execute_process(
          COMMAND ${_script}
          WORKING_DIRECTORY ${DEAL_II_SOURCE_DIR}
          OUTPUT_VARIABLE _tag
          RESULT_VARIABLE _result
          OUTPUT_STRIP_TRAILING_WHITESPACE
          )
       if(${_result} EQUAL 0)
         set(${_prefix}GIT_TAG ${_tag})
       endif()
    else()
       message(STATUS "Could not locate get_latest_tag.sh. " ${_prefix}GIT_TAG " will not be set.")
    endif()

  endif()

endmacro()

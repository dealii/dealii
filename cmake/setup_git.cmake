## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal.II authors
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
# Query git commit and branch information. If git was found and the source
# directory is under version control the following variables are populated:
#
#     DEAL_II_GIT_BRANCH
#     DEAL_II_GIT_REVISION
#     DEAL_II_GIT_SHORTREV
#

FIND_PACKAGE(Git)

#
# Only run the following if we have git and the source directory seems to
# be under version control.
#
IF(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git/HEAD)
  #
  # Bogus configure_file calls to trigger a reconfigure, and thus an update
  # of branch and commit information every time HEAD has changed.
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
  # Query for branch and commit information:
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
      "\\1" DEAL_II_GIT_REVISION "${_info}")
    STRING(REGEX REPLACE "^\"([^ ]+) ([^ ]+)\"$"
      "\\2" DEAL_II_GIT_SHORTREV "${_info}")
  ENDIF()

  EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
     OUTPUT_VARIABLE _branch
     RESULT_VARIABLE _result
     OUTPUT_STRIP_TRAILING_WHITESPACE
     )
  IF(${_result} EQUAL 0)
    STRING(REGEX REPLACE "refs/heads/" "" DEAL_II_GIT_BRANCH "${_branch}")
  ENDIF()
ENDIF()

FILE(WRITE ${CMAKE_BINARY_DIR}/revision.log
"###
#
#  Git information:
#        Branch:   ${DEAL_II_GIT_BRANCH}
#        Revision: ${DEAL_II_GIT_REVISION}
#
###"
  )

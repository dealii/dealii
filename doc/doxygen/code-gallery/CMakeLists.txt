## ---------------------------------------------------------------------
##
## Copyright (C) 2015 by the deal.II Authors
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
# A target for the preparation of all the stuff happening in here...
#

ADD_CUSTOM_TARGET(code-gallery)


#
# Check whether someone has either specified a global variable that
# points to a code gallery directory, or whether it has been put into
# a top-level code-gallery/ directory alongside the tutorial/
# directory
#
SET_IF_EMPTY(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
IF (EXISTS ${DEAL_II_CODE_GALLERY_DIRECTORY}/README.md)

  MESSAGE(STATUS "Setting up code gallery documentation from ${DEAL_II_CODE_GALLERY_DIRECTORY}.")

  # Collect the names of all code gallery projects. To
  # do so, find all 'author' files, then strip the last two
  # levels of these paths.
  #
  # For unclear reasons, the glob returns these files as
  # "/a/b/c/name//doc//author", so make sure we eat the
  # double slashes in the second step
  FILE(GLOB _code_gallery_names
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
  STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

  # Now set up targets for each of the code gallery programs

  # ...TODO...

ELSE()
  MESSAGE(STATUS "Generating the code-gallery documentation is only possible if the code gallery exists in ${DEAL_II_CODE_GALLERY_DIRECTORY}.")
ENDIF()


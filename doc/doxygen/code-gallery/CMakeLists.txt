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
# Define ways to get the code-gallery repository downloaded at
# the time that we call cmake, so that we can subsequently
# use the files that we find during download as inputs for
# doxygen. This has to happen at the time we run cmake because
# doing so at build time is too late for us to set up the
# proper dependencies for the build targets.
#
# To do so, we need to execute git on a command line; this
# (and consequently everything else below) will not work
# on Windows. It will also not work if we can't find a version
# of git
#
FIND_PACKAGE(Git)
IF( (NOT CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
    AND
    (NOT CMAKE_SYSTEM_NAME MATCHES "Windows")
    AND GIT_FOUND )

  # First update or download the code-gallery repository
  IF (IS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/code-gallery)
    MESSAGE(STATUS "Updating code-gallery git repository")
    EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} checkout master
                    COMMAND git pull origin master
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/code-gallery
                    TIMEOUT 30
                    OUTPUT_QUIET
                    ERROR_QUIET
		    RESULT_VARIABLE _result)
  ELSE()
    MESSAGE(STATUS "Cloning code-gallery git repository")
    EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} clone https://github.com/bangerth/code-gallery.git
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    TIMEOUT 30
                    OUTPUT_QUIET
                    ERROR_QUIET
		    RESULT_VARIABLE _result)
  ENDIF()

  # See if getting the git repository update succeeded
  IF (NOT ${_result})
    MESSAGE(STATUS "Updating code-gallery git repository - Done")

    # Now collect the names of all code gallery projects. To
    # do so, find all 'author' files, then strip the last two
    # levels of these paths.
    #
    # For unclear reasons, the glob returns these files as
    # "/a/b/c/name//doc//author", so make sure we eat the
    # double slashes in the second step
    FILE(GLOB _code_gallery_names
         "${CMAKE_CURRENT_BINARY_DIR}/code-gallery/*/doc/author")
    STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

    # Now set up targets for each of the code gallery programs

    # ...TODO...

  ELSE()
    MESSAGE(WARNING "Could not clone or update the code-gallery repository ${_result}")
  ENDIF()
ELSE()
  MESSAGE(WARNING "Generating the code-gallery documentation is only possible on unix-like systems with 'git' installed.")
ENDIF()


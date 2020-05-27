## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
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
# Configuration for thread support in deal.II with the help of the CPP
# Taskflow library:
#


MACRO(FEATURE_CPP_TASKFLOW_FIND_EXTERNAL var)
  FIND_PACKAGE(CPP_TASKFLOW)

  IF(CPP_TASKFLOW_FOUND)
    SET(${var} TRUE)
  ENDIF()

  IF(CPP_TASKFLOW_VERSION VERSION_LESS "2.4")
    # Clear the previously determined version numbers to avoid confusion
    SET(CPP_TASKFLOW_VERSION "bundled")
    SET(CPP_TASKFLOW_VERSION_MAJOR "")
    SET(CPP_TASKFLOW_VERSION_MINOR "")

    MESSAGE(STATUS
      "The externally provided Cpp Taskflow library is older than version 2.4, "
      "which cannot be used with deal.II."
      )
    SET(CPP_TASKFLOW_ADDITIONAL_ERROR_STRING
      "The externally provided Cpp Taskflow library is older than version\n"
      "2.4, which is the oldest version compatible with deal.II."
      )
    SET(${var} FALSE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_CPP_TASKFLOW_CONFIGURE_BUNDLED)
  LIST(APPEND CPP_TASKFLOW_BUNDLED_INCLUDE_DIRS ${CPP_TASKFLOW_FOLDER}/include)
ENDMACRO()


CONFIGURE_FEATURE(CPP_TASKFLOW)


## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Configuration for the SLEPC library:
#

SET(FEATURE_SLEPC_DEPENDS PETSC)


MACRO(FEATURE_SLEPC_FIND_EXTERNAL var)
  FIND_PACKAGE(SLEPC)

  IF(SLEPC_FOUND)
    #
    # Check whether SLEPc and PETSc are compatible according to
    # SLEPc's rules: This is equivalent to asking if the VERSION_MAJOR
    # and VERSION_MINOR of PETSc and SLEPc are
    # equivalent; and where VERSION_SUBMINORs are allowed to differ.
    #
    IF( ("${SLEPC_VERSION_MAJOR}" STREQUAL "${PETSC_VERSION_MAJOR}")
       AND
       ("${SLEPC_VERSION_MINOR}" STREQUAL "${PETSC_VERSION_MINOR}"))
      SET(${var} TRUE)
    ELSE()

      MESSAGE(STATUS "Could not find a sufficient SLEPc installation: "
        "The SLEPc library must have the same version as the PETSc library."
        )
      SET(SLEPC_ADDITIONAL_ERROR_STRING
        "Could not find a sufficient SLEPc installation: "
        "The SLEPc library must have the same version as the PETSc library.\n"
        )

      UNSET(SLEPC_INCLUDE_DIR_ARCH CACHE)
      UNSET(SLEPC_INCLUDE_DIR_COMMON CACHE)
      UNSET(SLEPC_LIBRARY CACHE)
      SET(SLEPC_DIR "" CACHE PATH
        "An optional hint to a SLEPc directory"
        )
      MARK_AS_ADVANCED(CLEAR SLEPC_DIR)

      SET(${var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()


MACRO(FEATURE_SLEPC_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find the SLEPc library!\n"
    ${SLEPC_ADDITIONAL_ERROR_STRING}
    "Please ensure that the SLEPc library version 3.0.0 or newer is installed on your computer\n"
    "and the version is the same as the one of the installed PETSc library.\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "SLEPc installed with --prefix=<...> to a destination:\n"
    "    $ SLEPC_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DSLEPC_DIR=\"...\" <...>\n"
    "SLEPc compiled in source tree:\n"
    "    $ SLEPC_DIR=\"...\"\n"
    "    $ cmake -DSLEPC_DIR=\"...\"\n"
    "or set the relevant variables by hand in ccmake.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(SLEPC)

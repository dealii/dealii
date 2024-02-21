## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# Configuration for the SLEPC library:
#

set(FEATURE_SLEPC_DEPENDS PETSC)


macro(feature_slepc_find_external var)
  find_package(DEAL_II_SLEPC)

  if(SLEPC_FOUND)
    #
    # Check whether SLEPc and PETSc are compatible according to
    # SLEPc's rules: This is equivalent to asking if the VERSION_MAJOR
    # and VERSION_MINOR of PETSc and SLEPc are
    # equivalent; and where VERSION_SUBMINORs are allowed to differ.
    #
    if( ("${SLEPC_VERSION_MAJOR}" STREQUAL "${PETSC_VERSION_MAJOR}")
       AND
       ("${SLEPC_VERSION_MINOR}" STREQUAL "${PETSC_VERSION_MINOR}"))
      set(${var} TRUE)
    else()
      set(SLEPC_VERSION_STR "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}")
      set(PETSC_VERSION_STR "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}")
      message(STATUS "Could not find a sufficient SLEPc installation: "
        "The SLEPc library "
        ${SLEPC_VERSION_STR}
        " must have the same version as the PETSc library "
        ${PETSC_VERSION_STR}
        )
      set(SLEPC_ADDITIONAL_ERROR_STRING
        "Could not find a sufficient SLEPc installation: "
        "The SLEPc library must have the same version as the PETSc library.\n"
        )

      set(${var} FALSE)
    endif()
  endif()
endmacro()


macro(feature_slepc_error_message)
  message(FATAL_ERROR "\n"
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
endmacro()


configure_feature(SLEPC)

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
# Configuration for the trilinos library:
#


MACRO(FEATURE_TRILINOS_FIND_EXTERNAL var)
  FIND_PACKAGE(TRILINOS)

  IF(TRILINOS_FOUND)
    #
    # So, we have a library. Let's see whether we can use it:
    #
    SET(${var} TRUE)

    #
    # Set TRILINOS_DIR to something meaningful if empty
    #
    IF("${TRILINOS_DIR}" STREQUAL "")
      SET(TRILINOS_DIR "<system location>")
    ENDIF()

    #
    # Check whether all required modules of trilinos are installed:
    #
    MESSAGE(STATUS
      "Check whether the found trilinos package contains all required modules:"
      )

    FOREACH(_module
      Amesos Epetra Ifpack AztecOO Sacado Teuchos
      )
      ITEM_MATCHES(_module_found ${_module} ${Trilinos_PACKAGE_LIST})
      IF(_module_found)
        MESSAGE(STATUS "Found ${_module}")
      ELSE()
        MESSAGE(STATUS "Module ${_module} not found!")
        SET(_modules_missing "${_modules_missing} ${_module}")
        SET(${var} FALSE)
      ENDIF()
    ENDFOREACH()

    IF(NOT ${var})
      MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
        "Missing ${_modules_missing}"
        )
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "is missing one or more modules necessary for the deal.II Trilinos interfaces:\n"
        "  ${_modules_missing}\n\n"
        "Please re-install Trilinos with the missing Trilinos subpackages enabled.\n\n"
        )
    ENDIF()

    #
    # Trilinos 10.6 had quite a number of bugs we ran into, see
    # for example
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5062
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5319
    #
    # The same is unfortunately true for 10.8.[01]:
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5370
    #
    IF((TRILINOS_VERSION_MAJOR EQUAL 10 AND
        TRILINOS_VERSION_MINOR EQUAL 6)
       OR
       (TRILINOS_VERSION_MAJOR EQUAL 10 AND
        TRILINOS_VERSION_MINOR EQUAL 8 AND
        TRILINOS_VERSION_SUBMINOR LESS 2))

      MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
        "Version ${TRILINOS_VERSION_MAJOR}.${TRILINOS_VERSION_MINOR}.${TRILINOS_VERSION_SUBMINOR} has bugs that make "
        "it incompatible with deal.II. Please use versions before 10.6 or after 10.8.1"
        )
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "with version ${TRILINOS_VERSION_MAJOR}.${TRILINOS_VERSION_MINOR}.${TRILINOS_VERSION_SUBMINOR} has bugs that make\n"
        "it incompatible with deal.II. Please use versions before 10.6 or after\n"
        "10.8.1.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Trilinos has to be configured with the same MPI configuration as
    # deal.II.
    #
    IF( (TRILINOS_WITH_MPI AND NOT DEAL_II_WITH_MPI)
         OR
         (NOT TRILINOS_WITH_MPI AND DEAL_II_WITH_MPI))
      MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
        "Trilinos has to be configured with the same MPI configuration as deal.II."
        )
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "has to be configured with the same MPI configuration as deal.II, but found:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
        "  TRILINOS_WITH_MPI = ${TRILINOS_WITH_MPI}\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Trilinos has to be configured with 32bit indices if deal.II uses unsigned int.
    #
    IF(TRILINOS_WITH_NO_32BIT_INDICES AND NOT DEAL_II_WITH_64BIT_INDICES)
      MESSAGE(STATUS "deal.II was configured to use 32bit global indices but "
        "Trilinos was not."
        ) 
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "has to be configured to use the same number of bits as deal.II, but "
        "found:\n"
        "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
        "  TRILINOS_WITH_NO_32BIT_INDICES = ${TRILINOS_WITH_NO_32_BIT_INDICES}\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Trilinos has to be configured with 64bit indices if deal.II uses unsigned long
    # long int.
    #
    IF(TRILINOS_WITH_NO_64BIT_INDICES AND DEAL_II_WITH_64BIT_INDICES)
      MESSAGE(STATUS "deal.II was configured to use 64bit global indices but "
        "Trilinos was not."
        )
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "has to be configured to use the same number of bits as deal.II, but "
        "found:\n"
        "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
        "  TRILINOS_WITH_NO_64BIT_INDICES = ${TRILINOS_WITH_NO_64_BIT_INDICES}\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Some versions of Sacado_cmath.hpp do things that aren't compatible
    # with the -std=c++0x flag of GCC, see deal.II FAQ.
    # Test whether that is indeed the case
    #
    IF(DEAL_II_WITH_CXX11 AND NOT TRILINOS_SUPPORTS_CPP11)

      IF(TRILINOS_HAS_C99_TR1_WORKAROUND)
        LIST(APPEND TRILINOS_DEFINITIONS "HAS_C99_TR1_CMATH")
        LIST(APPEND TRILINOS_USER_DEFINITIONS "HAS_C99_TR1_CMATH")
      ELSE()
        MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
          "The installation is not compatible with the C++ standard selected for "
          "this compiler."
          )
        SET(TRILINOS_ADDITIONAL_ERROR_STRING
          ${TRILINOS_ADDITIONAL_ERROR_STRING}
          "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
          "is not compatible with the C++ standard selected for\n"
          "this compiler. See the deal.II FAQ page for a solution.\n\n"
          )
        SET(${var} FALSE)
      ENDIF()
    ENDIF()

    CHECK_MPI_INTERFACE(TRILINOS ${var})
  ENDIF()
ENDMACRO()


MACRO(FEATURE_TRILINOS_CONFIGURE_EXTERNAL)
  SET(DEAL_II_EXPAND_TRILINOS_VECTOR "TrilinosWrappers::Vector")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCKVECTOR "TrilinosWrappers::BlockVector")
  SET(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN "TrilinosWrappers::SparsityPattern")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN "TrilinosWrappers::BlockSparsityPattern")
  SET(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR "TrilinosWrappers::MPI::BlockVector")
  SET(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR "TrilinosWrappers::MPI::Vector")

  #
  # Disable a bunch of warnings caused by Trilinos headers:
  #
  ENABLE_IF_SUPPORTED(TRILINOS_CXX_FLAGS "-Wno-unused")
  ENABLE_IF_SUPPORTED(TRILINOS_CXX_FLAGS "-Wno-extra")
  ENABLE_IF_SUPPORTED(TRILINOS_CXX_FLAGS "-Wno-overloaded-virtual")
ENDMACRO()


CONFIGURE_FEATURE(TRILINOS)

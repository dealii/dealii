## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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
      Amesos Epetra Ifpack AztecOO Sacado Teuchos ML MueLu
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
    # We require at least Trilinos 12.4
    #
    IF(TRILINOS_VERSION VERSION_LESS 12.4)

      MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II requires at least version 12.4, but version ${TRILINOS_VERSION} was found."
        )
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "with version ${TRILINOS_VERSION} is too old.\n"
        "deal.II requires at least version 12.4.\n\n"
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
    # Trilinos has to be configured with 32bit indices if deal.II uses
    # unsigned int.
    #
    IF(TRILINOS_WITH_NO_32BIT_INDICES AND NOT DEAL_II_WITH_64BIT_INDICES)
      MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II was configured to use 32bit global indices but "
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
    # Trilinos has to be configured with 64bit indices if deal.II uses
    # unsigned long long int.
    #
    IF(TRILINOS_WITH_NO_64BIT_INDICES AND DEAL_II_WITH_64BIT_INDICES)
      MESSAGE(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II was configured to use 64bit global indices but "
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
    IF(NOT TRILINOS_SUPPORTS_CPP11)

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

    IF (${var})
      FOREACH(_optional_module ROL Zoltan)
      ITEM_MATCHES(_module_found ${_optional_module} ${Trilinos_PACKAGE_LIST})
      IF(_module_found)
          MESSAGE(STATUS "Found ${_optional_module}")
          STRING(TOUPPER "${_optional_module}" _optional_module_upper)
          SET(DEAL_II_TRILINOS_WITH_${_optional_module_upper} ON)
      ELSE()
          MESSAGE(STATUS "Module ${_optional_module} not found!")
      ENDIF()
      ENDFOREACH()
    ENDIF()
  ENDIF()
ENDMACRO()


MACRO(FEATURE_TRILINOS_CONFIGURE_EXTERNAL)
  SET(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN "TrilinosWrappers::SparsityPattern")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN "TrilinosWrappers::BlockSparsityPattern")
  SET(DEAL_II_EXPAND_TRILINOS_SPARSE_MATRICES 
      "TrilinosWrappers::SparseMatrix"
      "TrilinosWrappers::BlockSparseMatrix")
  SET(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR "TrilinosWrappers::MPI::BlockVector")
  SET(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR "TrilinosWrappers::MPI::Vector")
  # Note: Only CMake 3.0 and greater support line continuation with the "\" character
  #       Elements of string lists are naturally separated by a ";"
  SET(DEAL_II_EXPAND_TRILINOS_SACADO_TYPES_FAD
      "Sacado::Fad::DFad<double>"
      "Sacado::Fad::DFad<float>"
      "Sacado::Fad::DFad<Sacado::Fad::DFad<double>>"
      "Sacado::Fad::DFad<Sacado::Fad::DFad<float>>")
  SET(DEAL_II_EXPAND_TRILINOS_SACADO_TYPES_RAD
      "Sacado::Rad::ADvar<double>"
      "Sacado::Rad::ADvar<float>"
      "Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>"
      "Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>")

  IF (TRILINOS_WITH_MPI)
    SET(DEAL_II_EXPAND_EPETRA_VECTOR "LinearAlgebra::EpetraWrappers::Vector")
  ENDIF()
  IF (TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD)
    SET(DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD ${TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD})
  ENDIF()
ENDMACRO()


CONFIGURE_FEATURE(TRILINOS)

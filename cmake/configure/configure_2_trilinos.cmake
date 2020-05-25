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
        Amesos Epetra Ifpack AztecOO Teuchos ML
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
    # deal.II has to be configured with MPI if both Trilinos and PETSc are
    # enabled.
    #
    IF(DEAL_II_WITH_TRILINOS AND DEAL_II_WITH_PETSC AND NOT DEAL_II_WITH_MPI)
      MESSAGE(STATUS "Incompatible configuration settings: "
        "MPI must be enabled to use both Trilinos and PETSc, as both libraries "
        "provide mutually incompatible MPI stubs."
        )
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "Incompatible Trilinos and PETSc libraries found. Both libraries were "
        "configured without MPI support and cannot be used at the same time due "
        "to incompatible MPI stub files. Either reconfigure deal.II, Trilinos, "
        "and PETSc with MPI support, or disable one of the libraries.\n"
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

    CHECK_MPI_INTERFACE(TRILINOS ${var})

    #
    # Check which optional features of trilinos are installed.
    #
    IF (${var})
      #
      # Check for modules.
      #
      FOREACH(_optional_module EpetraExt ROL Sacado Tpetra MueLu Zoltan)
        ITEM_MATCHES(_module_found ${_optional_module} ${Trilinos_PACKAGE_LIST})
        IF(_module_found)
          MESSAGE(STATUS "Found ${_optional_module}")
          STRING(TOUPPER "${_optional_module}" _optional_module_upper)
          SET(DEAL_II_TRILINOS_WITH_${_optional_module_upper} ON)
        ELSE()
          MESSAGE(STATUS "Module ${_optional_module} not found!")
        ENDIF()
      ENDFOREACH()

      #
      # Check for third-party libraries (tpl).
      #
      FOREACH(_optional_tpl MUMPS)
        ITEM_MATCHES(_tpl_found ${_optional_tpl} ${Trilinos_TPL_LIST})
        IF(_tpl_found)
          MESSAGE(STATUS "Found ${_optional_tpl}")
          STRING(TOUPPER "${_optional_tpl}" _optional_tpl_upper)
          SET(DEAL_II_TRILINOS_WITH_${_optional_tpl_upper} ON)
        ELSE()
          MESSAGE(STATUS "Module ${_optional_tpl} not found!")
        ENDIF()
      ENDFOREACH()
    ENDIF()

    IF(DEAL_II_TRILINOS_WITH_TPETRA)
      #
      # Check if Tpetra is usable in fact.
      #
      LIST(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      LIST(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
      ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

      CHECK_SYMBOL_EXISTS(
        "KOKKOS_ENABLE_CUDA_LAMBDA"
        "Kokkos_Macros.hpp"
        DEAL_II_KOKKOS_LAMBDA_EXISTS
        )
      IF(DEAL_II_KOKKOS_LAMBDA_EXISTS)
        ADD_FLAGS(CMAKE_REQUIRED_FLAGS "--expt-extended-lambda")
      ENDIF()

      LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <Tpetra_Vector.hpp>
        int
        main()
        {
          using LO       = int;
          using GO       = unsigned int;
          using Node     = Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial>;
          using map_type = Tpetra::Map<LO, GO, Node>;
          Teuchos::RCP<const map_type>         dummy_map = Teuchos::rcp(new map_type());
          Tpetra::Vector<double, LO, GO, Node> dummy_vector(dummy_map);
          (void)dummy_vector;
          return 0;
        }
        "
        TRILINOS_TPETRA_IS_FUNCTIONAL
        )

      RESET_CMAKE_REQUIRED()

      IF(NOT TRILINOS_TPETRA_IS_FUNCTIONAL)
        MESSAGE(
          STATUS
          "Tpetra was found but is not usable! Disabling Tpetra support."
          )
        SET(DEAL_II_TRILINOS_WITH_TPETRA OFF)
      ENDIF()
    ENDIF()

    IF(DEAL_II_TRILINOS_WITH_MUELU)
      #
      # Check if MueLu is actually usable.
      #
      LIST(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      LIST(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
      ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

      LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <MueLu_CreateEpetraPreconditioner.hpp>
        int
        main()
        {
          Epetra_CrsMatrix *matrix;
          const auto teuchos_wrapped_matrix = Teuchos::rcp(matrix, false);	
          Teuchos::ParameterList parameters;
          MueLu::CreateEpetraPreconditioner(teuchos_wrapped_matrix, parameters);
          return 0;
        }
        "
        TRILINOS_MUELU_IS_FUNCTIONAL
        )

      RESET_CMAKE_REQUIRED()

      IF(NOT TRILINOS_MUELU_IS_FUNCTIONAL)
        MESSAGE(
          STATUS
          "MueLu was found but is not usable through Epetra! Disabling MueLu support."
          )
        SET(DEAL_II_TRILINOS_WITH_MUELU OFF)
      ENDIF()
    ENDIF()

    IF(${DEAL_II_TRILINOS_WITH_SACADO})
      #
      # Look for Sacado_config.h - we'll query it to determine C++11 support:
      #
      DEAL_II_FIND_FILE(SACADO_CONFIG_H Sacado_config.h
        HINTS ${Trilinos_INCLUDE_DIRS}
        NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
        NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
        )

      #
      # GCC 6.3.0 has a bug that prevents the creation of complex
      # numbers templated on Sacado::Rad::ADvar types:
      #
      # include/c++/6.3.0/complex: In instantiation of
      # ‘struct std::complex<Sacado::Rad::ADvar<double> >’:
      # include/c++/6.3.0/complex:206:16: error: ‘std::complex<_Tp>& std::complex<_Tp>::operator=(const std::complex<_Tp>&) [with _Tp = Sacado::Rad::ADvar<double>]’ declared to take const reference, but implicit declaration would take non-const
      #
      # Test whether the compiler hits this issue
      #
      DEAL_II_FIND_FILE(SACADO_TRAD_HPP Sacado_trad.hpp
        HINTS ${Trilinos_INCLUDE_DIRS}
        NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
        NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
        )

      IF(EXISTS ${SACADO_TRAD_HPP})
        LIST(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
        ADD_FLAGS(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_VERSION_FLAG}")

        CHECK_CXX_SOURCE_COMPILES(
          "
          #include <Sacado_trad.hpp>
          #include <complex>
          int main ()
          {
            Sacado::Rad::ADvar<double> sacado_rad_double; // Works
            std::complex<Sacado::Rad::ADvar<double> > complex_sacado_rad_double; // Doesn't work
            return 0;
          }
          "
          TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
          )
        RESET_CMAKE_REQUIRED()
      ENDIF()

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
  IF (TRILINOS_WITH_MPI)
    SET(DEAL_II_EXPAND_EPETRA_VECTOR "LinearAlgebra::EpetraWrappers::Vector")
    IF (${DEAL_II_TRILINOS_WITH_TPETRA})
      SET(DEAL_II_EXPAND_TPETRA_VECTOR_DOUBLE "LinearAlgebra::TpetraWrappers::Vector<double>")
      SET(DEAL_II_EXPAND_TPETRA_VECTOR_FLOAT "LinearAlgebra::TpetraWrappers::Vector<float>")
      IF (${DEAL_II_WITH_COMPLEX_NUMBERS})
        SET(DEAL_II_EXPAND_TPETRA_VECTOR_COMPLEX_DOUBLE "LinearAlgebra::TpetraWrappers::Vector<std::complex<double>>")
        SET(DEAL_II_EXPAND_TPETRA_VECTOR_COMPLEX_FLOAT "LinearAlgebra::TpetraWrappers::Vector<std::complex<float>>")
      ENDIF()
    ENDIF()
  ENDIF()
  IF(${DEAL_II_TRILINOS_WITH_SACADO})
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

    IF (TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD)
      SET(DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD ${TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD})
    ENDIF()
  ENDIF()
ENDMACRO()


CONFIGURE_FEATURE(TRILINOS)

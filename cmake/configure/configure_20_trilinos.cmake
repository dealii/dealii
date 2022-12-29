## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2021 by the deal.II authors
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

set(FEATURE_TRILINOS_DEPENDS MPI)

macro(feature_trilinos_find_external var)
  find_package(DEAL_II_TRILINOS)

  if(TRILINOS_FOUND)
    #
    # So, we have a library. Let's see whether we can use it:
    #
    set(${var} TRUE)

    #
    # Set TRILINOS_DIR to something meaningful if empty
    #
    if("${TRILINOS_DIR}" STREQUAL "")
      set(TRILINOS_DIR "<system location>")
    endif()

    #
    # Check whether all required modules of trilinos are installed:
    #
    message(STATUS
      "Checking whether the found trilinos package contains all required modules:"
      )

    foreach(_module
        Amesos Epetra Ifpack AztecOO Teuchos ML
      )
      item_matches(_module_found ${_module} ${Trilinos_PACKAGE_LIST})
      if(_module_found)
        message(STATUS "  Found ${_module}")
      else()
        message(STATUS "  Module ${_module} not found!")
        set(_modules_missing "${_modules_missing} ${_module}")
        set(${var} FALSE)
      endif()
    endforeach()

    if(NOT ${var})
      message(STATUS "Could not find a sufficient Trilinos installation: "
        "Missing ${_modules_missing}"
        )
      set(TRILINOS_ADDITIONAL_ERROR_STRING
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "is missing one or more modules necessary for the deal.II Trilinos interfaces:\n"
        "  ${_modules_missing}\n\n"
        "Please re-install Trilinos with the missing Trilinos subpackages enabled.\n\n"
        )
    endif()

    #
    # We require at least Trilinos 12.4
    #
    if(TRILINOS_VERSION VERSION_LESS 12.4)

      message(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II requires at least version 12.4, but version ${TRILINOS_VERSION} was found."
        )
      set(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "with version ${TRILINOS_VERSION} is too old.\n"
        "deal.II requires at least version 12.4.\n\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Trilinos has to be configured with the same MPI configuration as
    # deal.II.
    #
    if(NOT TRILINOS_WITH_MPI)
      message(STATUS "Could not find a sufficient Trilinos installation: "
        "Trilinos has to have MPI support enabled."
        )
      set(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "has to be configured with MPI support, but found:\n"
        "  TRILINOS_WITH_MPI = ${TRILINOS_WITH_MPI}\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Trilinos has to be configured with 32bit indices if deal.II uses
    # unsigned int.
    #
    if(TRILINOS_WITH_NO_32BIT_INDICES AND NOT DEAL_II_WITH_64BIT_INDICES)
      message(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II was configured to use 32bit global indices but "
        "Trilinos was not."
        )
      set(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "has to be configured to use the same number of bits as deal.II, but "
        "found:\n"
        "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
        "  TRILINOS_WITH_NO_32BIT_INDICES = ${TRILINOS_WITH_NO_32_BIT_INDICES}\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Trilinos has to be configured with 64bit indices if deal.II uses
    # unsigned long long int.
    #
    if(TRILINOS_WITH_NO_64BIT_INDICES AND DEAL_II_WITH_64BIT_INDICES)
      message(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II was configured to use 64bit global indices but "
        "Trilinos was not."
        )
      set(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "has to be configured to use the same number of bits as deal.II, but "
        "found:\n"
        "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
        "  TRILINOS_WITH_NO_64BIT_INDICES = ${TRILINOS_WITH_NO_64_BIT_INDICES}\n"
        )
      set(${var} FALSE)
    endif()

    check_mpi_interface(TRILINOS ${var})

    #
    # Check which optional features of trilinos are installed.
    #
    if (${var})
      #
      # Check for modules.
      #
      foreach(_optional_module Belos EpetraExt Kokkos MueLu NOX ROL Sacado SEACAS Tpetra Zoltan)
        item_matches(_module_found ${_optional_module} ${Trilinos_PACKAGE_LIST})
        if(_module_found)
          message(STATUS "  Found ${_optional_module}")
          string(TOUPPER "${_optional_module}" _optional_module_upper)
          set(DEAL_II_TRILINOS_WITH_${_optional_module_upper} ON)
        else()
          message(STATUS "  Module ${_optional_module} not found!")
        endif()
      endforeach()

      #
      # Check for third-party libraries (tpl).
      #
      foreach(_optional_tpl MUMPS)
        item_matches(_tpl_found ${_optional_tpl} ${Trilinos_TPL_LIST})
        if(_tpl_found)
          message(STATUS "  Found ${_optional_tpl}")
          string(TOUPPER "${_optional_tpl}" _optional_tpl_upper)
          set(DEAL_II_TRILINOS_WITH_${_optional_tpl_upper} ON)
        else()
          message(STATUS "  Module ${_optional_tpl} not found!")
        endif()
      endforeach()
    endif()

    if(DEAL_II_TRILINOS_WITH_KOKKOS)
      if(DEAL_II_FORCE_BUNDLED_KOKKOS)
        set(TRILINOS_ADDITIONAL_ERROR_STRING
          ${TRILINOS_ADDITIONAL_ERROR_STRING}
          "The Trilinos installation (found at \"${TRILINOS_DIR}\")"
          "includes Kokkos, but DEAL_II_FORCE_BUNDLED_KOKKOS=ON!\n")
        set(${var} FALSE)
      endif()

      if(Kokkos_ENABLE_CUDA)
        # We need to disable SIMD vectorization for CUDA device code.
        # Otherwise, nvcc compilers from version 9 on will emit an error message like:
        # "[...] contains a vector, which is not supported in device code". We
        # would like to set the variable in check_01_cpu_feature but at that point
        # we don't know if CUDA support is enabled in Kokkos
        set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 0)
      endif()

      # We need a recent version of Trilinos to use kokkos_check. We want to use
      # VERSION_GREATER_EQUAL 13.0 but this requires CMake 3.7
      if(TRILINOS_VERSION VERSION_GREATER 12.99 AND Kokkos_ENABLE_CUDA)
        KOKKOS_CHECK(OPTIONS CUDA_LAMBDA)
      endif()
    endif()

    if(DEAL_II_TRILINOS_WITH_TPETRA)
      #
      # Check if Tpetra is usable in fact.
      #
      list(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      list(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})

      list(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

      # For the case of Trilinos being compiled with openmp support the
      # following Tpetra test needs -fopenmp to succeed. Make sure that we
      # supply the correct compiler and linker flags:
      add_flags(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_LINKER_FLAGS}")

      if(DEAL_II_WITH_64BIT_INDICES)
        set(_global_index_type "std::uint64_t")
      else()
        set(_global_index_type "unsigned int")
      endif()

      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <cstdint>
        #include <Tpetra_Vector.hpp>
        int
        main()
        {
          using LO       = int;
          using GO       = ${_global_index_type};
          using map_type = Tpetra::Map<LO, GO>;
          Teuchos::RCP<const map_type>   dummy_map = Teuchos::rcp(new map_type());
          Tpetra::Vector<double, LO, GO> dummy_vector(dummy_map);
          (void)dummy_vector;
          return 0;
        }
        "
        TRILINOS_TPETRA_IS_FUNCTIONAL
        )

      reset_cmake_required()

      if(NOT TRILINOS_TPETRA_IS_FUNCTIONAL)
        message(
          STATUS
          "Tpetra was found but is not usable! Disabling Tpetra support."
          )
        set(DEAL_II_TRILINOS_WITH_TPETRA OFF)
      endif()
    endif()

    if(DEAL_II_TRILINOS_WITH_MUELU)
      #
      # Check if MueLu is actually usable.
      #
      list(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      list(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})

      list(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

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

      reset_cmake_required()

      if(NOT TRILINOS_MUELU_IS_FUNCTIONAL)
        message(
          STATUS
          "MueLu was found but is not usable through Epetra! Disabling MueLu support."
          )
        set(DEAL_II_TRILINOS_WITH_MUELU OFF)
      endif()
    endif()

    # the only thing we use from SEACAS right now is ExodusII, so just check
    # that it works
    if(${DEAL_II_TRILINOS_WITH_SEACAS})
      list(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      list(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES})
      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <exodusII.h>
        int
        main()
        {
          int component_word_size = sizeof(double);
          int floating_point_word_size = 0;
          float ex_version = 0;
          const int ex_id = ex_open(\"test.ex\",
                                    EX_READ,
                                    &component_word_size,
                                    &floating_point_word_size,
                                    &ex_version);
          ex_close(ex_id);
          return 0;
        }
        "
        TRILINOS_SEACAS_IS_FUNCTIONAL
        )

      reset_cmake_required()

      if(NOT TRILINOS_SEACAS_IS_FUNCTIONAL)
        message(
          STATUS
          "SEACAS was found but doesn't seem to include ExodusII. Disabling SEACAS support."
          )
        set(DEAL_II_TRILINOS_WITH_SEACAS OFF)
      endif()
    endif()

    if(${DEAL_II_TRILINOS_WITH_SACADO})
      #
      # Look for Sacado_config.h - we'll query it to determine C++11 support:
      #
      deal_ii_find_file(SACADO_CONFIG_H Sacado_config.h
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
      deal_ii_find_file(SACADO_TRAD_HPP Sacado_trad.hpp
        HINTS ${Trilinos_INCLUDE_DIRS}
        NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
        NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
        )

      if(EXISTS ${SACADO_TRAD_HPP})
        list(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})

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
        reset_cmake_required()
      endif()

    endif()
  endif()
endmacro()


macro(feature_trilinos_configure_external)
  set(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN "TrilinosWrappers::SparsityPattern")
  set(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN "TrilinosWrappers::BlockSparsityPattern")
  set(DEAL_II_EXPAND_TRILINOS_SPARSE_MATRICES 
      "TrilinosWrappers::SparseMatrix"
      "TrilinosWrappers::BlockSparseMatrix")
  set(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR "TrilinosWrappers::MPI::BlockVector")
  set(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR "TrilinosWrappers::MPI::Vector")
  set(DEAL_II_EXPAND_EPETRA_VECTOR "LinearAlgebra::EpetraWrappers::Vector")
  if (${DEAL_II_TRILINOS_WITH_TPETRA})
    set(DEAL_II_EXPAND_TPETRA_VECTOR_DOUBLE "LinearAlgebra::TpetraWrappers::Vector<double>")
    set(DEAL_II_EXPAND_TPETRA_VECTOR_FLOAT "LinearAlgebra::TpetraWrappers::Vector<float>")
    if (${DEAL_II_WITH_COMPLEX_NUMBERS})
      set(DEAL_II_EXPAND_TPETRA_VECTOR_COMPLEX_DOUBLE "LinearAlgebra::TpetraWrappers::Vector<std::complex<double>>")
      set(DEAL_II_EXPAND_TPETRA_VECTOR_COMPLEX_FLOAT "LinearAlgebra::TpetraWrappers::Vector<std::complex<float>>")
    endif()
  endif()
  if(${DEAL_II_TRILINOS_WITH_SACADO})
    # Note: Only CMake 3.0 and greater support line continuation with the "\" character
    #       Elements of string lists are naturally separated by a ";"
    set(DEAL_II_EXPAND_TRILINOS_SACADO_TYPES_FAD
        "Sacado::Fad::DFad<double>"
        "Sacado::Fad::DFad<float>"
        "Sacado::Fad::DFad<Sacado::Fad::DFad<double>>"
        "Sacado::Fad::DFad<Sacado::Fad::DFad<float>>")
    set(DEAL_II_EXPAND_TRILINOS_SACADO_TYPES_RAD
        "Sacado::Rad::ADvar<double>"
        "Sacado::Rad::ADvar<float>"
        "Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>"
        "Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>")

    if (TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD)
      set(DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD ${TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD})
    endif()
  endif()
endmacro()


configure_feature(TRILINOS)

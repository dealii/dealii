## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# Configuration for the trilinos library:
#

set(FEATURE_TRILINOS_DEPENDS MPI)

#
# A list of optional Trilinos modules we use:
#
set(_deal_ii_trilinos_optional_modules
  Amesos2 Belos EpetraExt Ifpack2 Kokkos MueLu NOX ROL Sacado SEACAS Tpetra Zoltan 
)

#
# A list of optional Trilinos TPLs we use:
#
set(_deal_ii_trilinos_optional_tpls MUMPS)

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
    # We require at least Trilinos 13.2
    #
    if(TRILINOS_VERSION VERSION_LESS 13.2)
      message(STATUS "Could not find a sufficient Trilinos installation: "
        "deal.II requires at least version 13.2, but version ${TRILINOS_VERSION} was found."
      )
      set(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation (found at \"${TRILINOS_DIR}\")\n"
        "with version ${TRILINOS_VERSION} is too old.\n"
        "deal.II requires at least version 13.2.\n\n"
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
    if(${var})
      #
      # Check for modules.
      #
      foreach(_optional_module ${_deal_ii_trilinos_optional_modules})
        item_matches(_module_found ${_optional_module} ${Trilinos_PACKAGE_LIST})
        if(_module_found)
          message(STATUS "  Found ${_optional_module}")
          string(TOUPPER "${_optional_module}" _optional_module_upper)
          set(TRILINOS_WITH_${_optional_module_upper} ON)
        else()
          message(STATUS "  Module ${_optional_module} not found!")
        endif()
      endforeach()

      #
      # Check for third-party libraries (tpl).
      #
      foreach(_optional_tpl ${_deal_ii_trilinos_optional_tpls})
        item_matches(_tpl_found ${_optional_tpl} ${Trilinos_TPL_LIST})
        if(_tpl_found)
          message(STATUS "  Found ${_optional_tpl}")
          string(TOUPPER "${_optional_tpl}" _optional_tpl_upper)
          set(TRILINOS_WITH_${_optional_tpl_upper} ON)
        else()
          message(STATUS "  Module ${_optional_tpl} not found!")
        endif()
      endforeach()
    endif()

    if(TRILINOS_WITH_KOKKOS)
      if(DEAL_II_FORCE_BUNDLED_KOKKOS)
        set(TRILINOS_ADDITIONAL_ERROR_STRING
          ${TRILINOS_ADDITIONAL_ERROR_STRING}
          "The Trilinos installation (found at \"${TRILINOS_DIR}\")"
          "includes Kokkos, but DEAL_II_FORCE_BUNDLED_KOKKOS=ON!\n")
        set(${var} FALSE)
      endif()

      #
      # When configuring Kokkos we have to ensure that we actually pick up the
      # correct Kokkos installation coming from Trilinos.
      #
      # FIXME: this logic should probably be refactored into
      # FindDEAL_II_TRILINOS.cmake...
      #
      set(TRILINOS_KOKKOS_DIR "${TRILINOS_CONFIG_DIR}/..")
    endif()

    if(TRILINOS_WITH_TPETRA)
      #
      # Check if Tpetra is usable in fact.
      #
      list(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      list(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})

      list(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})
      list(APPEND CMAKE_REQUIRED_FLAGS ${TRILINOS_CXX_FLAGS})

      # For the case of Trilinos being compiled with openmp support the
      # following Tpetra test needs -fopenmp to succeed. Make sure that we
      # supply the correct compiler and linker flags:
      add_flags(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_LINKER_FLAGS}")

      if(DEAL_II_WITH_64BIT_INDICES)
        set(_global_index_type "long long")
      else()
        set(_global_index_type "int")
      endif()

      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <cstdint>
        #include <Tpetra_Vector.hpp>
        int main()
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

      if(TRILINOS_TPETRA_IS_FUNCTIONAL)
        #
        # We need to figure out what instantiations are used in Tpetra so
        # that we can populate our DEAL_II_EXPAND_TPETRA correctly. We need
        # to dof this here prior to the call to reset_cmake_required().
        #
        check_cxx_symbol_exists(HAVE_TPETRA_INST_DOUBLE "TpetraCore_config.h" _tpetra_inst_double)
        check_cxx_symbol_exists(HAVE_TPETRA_INST_FLOAT "TpetraCore_config.h" _tpetra_inst_float)
        check_cxx_symbol_exists(HAVE_TPETRA_INST_COMPLEX_DOUBLE "TpetraCore_config.h" _tpetra_inst_complex_double)
        check_cxx_symbol_exists(HAVE_TPETRA_INST_COMPLEX_FLOAT "TpetraCore_config.h" _tpetra_inst_complex_float)
      else()
        message(STATUS "Tpetra was found but is not usable due to a mismatch in ordinal number types.")
        set(TRILINOS_WITH_TPETRA OFF)

        check_cxx_symbol_exists(HAVE_TPETRA_INT_INT "TpetraCore_config.h" _tpetra_int_int)
        check_cxx_symbol_exists(HAVE_TPETRA_INT_LONG_LONG "TpetraCore_config.h" _tpetra_int_long_long)

        if(NOT _tpetra_int_long_long AND DEAL_II_WITH_64BIT_INDICES)
          message( STATUS
            "  Tpetra was configured *without* support for 64-bit global indices"
            " but deal.II is configured to use 64-bit global indices."
            )
          message(STATUS
            "  Either reconfigure deal.II with -DDEAL_II_WITH_64BIT_INDICES=OFF"
            " or rebuild Trilinos with -DTpetra_INST_INT_LONG_LONG=ON"
            )

        elseif(NOT _tpetra_int_int AND NOT DEAL_II_WITH_64BIT_INDICES)
          message( STATUS
            "  Tpetra was configured *without* support for 32-bit global indices"
            " but deal.II is configured to use 32-bit global indices."
            )
          message(STATUS
            "  Either reconfigure deal.II with -DDEAL_II_WITH_64BIT_INDICES=ON"
            " or rebuild Trilinos with -DTpetra_INST_INT_INT=ON"
            )
        endif()

        reset_cmake_required()
      endif()
    endif()

    if(TRILINOS_WITH_MUELU)
      #
      # Check if MueLu is actually usable.
      #
      list(APPEND CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
      list(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})

      list(APPEND CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})
      list(APPEND CMAKE_REQUIRED_FLAGS ${TRILINOS_CXX_FLAGS})

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
        set(TRILINOS_WITH_MUELU OFF)
      endif()
    endif()

    # the only thing we use from SEACAS right now is ExodusII, so just check
    # that it works
    if(${TRILINOS_WITH_SEACAS})
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
        set(TRILINOS_WITH_SEACAS OFF)
      endif()
    endif()

    if(${TRILINOS_WITH_SACADO})
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
  #
  # Propagate optional Trilinos modules and TPLs into DEAL_II namespace:
  #

  foreach(_module ${_deal_ii_trilinos_optional_modules} ${_deal_ii_trilinos_optional_tpls})
    string(TOUPPER "${_module}" _module_upper)
    if(${TRILINOS_WITH_${_module_upper}})
      set(DEAL_II_TRILINOS_WITH_${_module_upper} ON)
    endif()
  endforeach()

  #
  # Figure out all the possible instantiations we need:
  #

  set(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN "TrilinosWrappers::SparsityPattern")
  set(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN "TrilinosWrappers::BlockSparsityPattern")
  set(DEAL_II_EXPAND_TRILINOS_SPARSE_MATRICES "TrilinosWrappers::SparseMatrix" "TrilinosWrappers::BlockSparseMatrix")
  set(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR "TrilinosWrappers::MPI::BlockVector")
  set(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR "TrilinosWrappers::MPI::Vector")
  set(DEAL_II_EXPAND_EPETRA_VECTOR "LinearAlgebra::EpetraWrappers::Vector")

  if(${DEAL_II_TRILINOS_WITH_TPETRA})
    if(_tpetra_inst_double)
      set(DEAL_II_EXPAND_TPETRA_TYPES "double")
      set(DEAL_II_EXPAND_TPETRA_VECTOR_DOUBLE
        "LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host>"
        "LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>")
      set(DEAL_II_EXPAND_TPETRA_BLOCKVECTOR_DOUBLE
        "LinearAlgebra::TpetraWrappers::BlockVector<double, MemorySpace::Host>"
        "LinearAlgebra::TpetraWrappers::BlockVector<double, MemorySpace::Default>")
    endif()

    if(_tpetra_inst_float)
      set(DEAL_II_EXPAND_TPETRA_VECTOR_FLOAT
        "LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Host>"
        "LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default>")
      set(DEAL_II_EXPAND_TPETRA_BLOCKVECTOR_FLOAT
        "LinearAlgebra::TpetraWrappers::BlockVector<float, MemorySpace::Host>"
        "LinearAlgebra::TpetraWrappers::BlockVector<float, MemorySpace::Default>")
    endif()

    if(${DEAL_II_WITH_COMPLEX_NUMBERS})
      if(_tpetra_inst_complex_double)
        set(DEAL_II_EXPAND_TPETRA_VECTOR_COMPLEX_DOUBLE
          "LinearAlgebra::TpetraWrappers::Vector<std::complex<double>, MemorySpace::Host>"
          "LinearAlgebra::TpetraWrappers::Vector<std::complex<double>, MemorySpace::Default>")
        set(DEAL_II_EXPAND_TPETRA_BLOCKVECTOR_COMPLEX_DOUBLE
          "LinearAlgebra::TpetraWrappers::BlockVector<std::complex<double>, MemorySpace::Host>"
          "LinearAlgebra::TpetraWrappers::BlockVector<std::complex<double>, MemorySpace::Default>")
      endif()

      if(_tpetra_inst_complex_float)
        set(DEAL_II_EXPAND_TPETRA_VECTOR_COMPLEX_FLOAT
          "LinearAlgebra::TpetraWrappers::Vector<std::complex<float>, MemorySpace::Host>"
          "LinearAlgebra::TpetraWrappers::Vector<std::complex<float>, MemorySpace::Default>")
        set(DEAL_II_EXPAND_TPETRA_BLOCKVECTOR_COMPLEX_FLOAT
          "LinearAlgebra::TpetraWrappers::BlockVector<std::complex<float>, MemorySpace::Host>"
          "LinearAlgebra::TpetraWrappers::BlockVector<std::complex<float>, MemorySpace::Default>")
      endif()
    endif()
  endif()

  if(${DEAL_II_TRILINOS_WITH_SACADO})
    set(DEAL_II_EXPAND_TRILINOS_SACADO_TYPES_FAD
      "Sacado::Fad::DFad<double>"
      "Sacado::Fad::DFad<float>"
      "Sacado::Fad::DFad<Sacado::Fad::DFad<double>>"
      "Sacado::Fad::DFad<Sacado::Fad::DFad<float>>"
      )
    set(DEAL_II_EXPAND_TRILINOS_SACADO_TYPES_RAD
      "Sacado::Rad::ADvar<double>"
      "Sacado::Rad::ADvar<float>"
      "Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>"
      "Sacado::Rad::ADvar<Sacado::Fad::DFad<float>>"
      )
    if(TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD)
      set(DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD ${TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD})
    endif()
  endif()
endmacro()


configure_feature(TRILINOS)

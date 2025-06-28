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
# Configuration for the petsc library:
#

set(FEATURE_PETSC_DEPENDS MPI)


macro(feature_petsc_find_external var)
  find_package(DEAL_II_PETSC)

  if(PETSC_FOUND)
    #
    # So, we have found a petsc library. Let's check whether we can use it.
    #
    set(${var} TRUE)

    #
    # We support petsc from version 3.7.x onwards
    #
    if(${PETSC_VERSION} VERSION_LESS 3.7.0)
      message(STATUS "Could not find a sufficiently modern PETSc installation: "
        "Version >=3.7.0 required!"
        )
      set(PETSC_ADDITIONAL_ERROR_STRING
        "Could not find a sufficiently modern PETSc installation: "
        "Version >=3.7.0 required!\n"
        )
      set(${var} FALSE)
    endif()

    #
    # PETSc has to be configured with the same MPI configuration as
    # deal.II.
    #
    # petscconf.h should export PETSC_HAVE_MPIUNI 1 in case mpi support is
    # _NOT_ enabled.
    # So we check for this:
    #
    if(PETSC_WITH_MPIUNI)
      message(STATUS "Could not find a sufficient PETSc installation: "
        "PETSc has to be configured with MPI support."
        )
      set(PETSC_ADDITIONAL_ERROR_STRING
        ${PETSC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient PETSc installation:\n"
        "PETSc has to be configured with the MPI support, but found:\n"
        "  PETSC_WITH_MPI   = (NOT ${PETSC_WITH_MPIUNI})\n"
        )
      set(${var} FALSE)
    endif()

    # If PETSc is compiled with complex scalar type we need to have support
    # for complex values within deal.II as well.
    #
    if(PETSC_WITH_COMPLEX AND NOT DEAL_II_WITH_COMPLEX_VALUES)
        message(STATUS "The PETSc configuration is incompatible with the deal.II configuration: "
          "PETSc is compiled with complex scalar type. "
          "This requires support for complex values in deal.II as well."
          )
      set(PETSC_ADDITIONAL_ERROR_STRING
        ${PETSC_ADDITIONAL_ERROR_STRING}
        "The PETSc configuration is incompatible with the deal.II configuration:\n"
        "PETSc is compiled with complex scalar type. "
        "This requires support for complex values in deal.II as well.\n"
        "  DEAL_II_WITH_COMPLEX_VALUES = ${DEAL_II_WITH_COMPLEX_VALUES}\n"
        "  PETSC_WITH_COMPLEX = (${PETSC_WITH_COMPLEX})\n"
        )
      set(${var} FALSE)
    endif()

    if(PETSC_WITH_KOKKOS)
      if(DEAL_II_FORCE_BUNDLED_KOKKOS)
        set(PETSC_ADDITIONAL_ERROR_STRING
          ${PETSC_ADDITIONAL_ERROR_STRING}
          "The PETSc installation (found at \"${PETSC_DIR}\")"
          "includes Kokkos, but DEAL_II_FORCE_BUNDLED_KOKKOS=ON!\n")
        set(${var} FALSE)
      endif()
    endif()


    check_mpi_interface(PETSC ${var})
  endif()
endmacro()


macro(feature_petsc_configure_external)
  #
  # Propagate some PETSc configuration variables into DEAL_II namespace:
  #

  set(DEAL_II_PETSC_WITH_COMPLEX ${PETSC_WITH_COMPLEX})
  set(DEAL_II_PETSC_WITH_HYPRE ${PETSC_WITH_HYPRE})
  set(DEAL_II_PETSC_WITH_MUMPS ${PETSC_WITH_MUMPS})
  set(DEAL_II_PETSC_WITH_KOKKOS ${PETSC_WITH_KOKKOS})

  #
  # Figure out all the possible instantiations we need:
  #

  set(DEAL_II_EXPAND_PETSC_MPI_VECTOR "PETScWrappers::MPI::Vector")
  set(DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR "PETScWrappers::MPI::BlockVector")
  set(DEAL_II_EXPAND_PETSC_SPARSE_MATRICES
      "PETScWrappers::SparseMatrix"
      "PETScWrappers::MPI::SparseMatrix"
      "PETScWrappers::MPI::BlockSparseMatrix")
  #
  # FIXME:
  # temporary variable until deal.II fully support complex-valued PETSc
  if( NOT PETSC_WITH_COMPLEX )
    set(DEAL_II_EXPAND_PETSC_MPI_VECTOR_REAL "PETScWrappers::MPI::Vector")
    set(DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR_REAL "PETScWrappers::MPI::BlockVector")
  else()
    message(STATUS "Compiling with complex-valued algebra")
  endif()
endmacro()


macro(feature_petsc_error_message)
  message(FATAL_ERROR "\n"
    "Could not find the petsc library!\n"
    ${PETSC_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that the petsc library version 3.7.0 or newer is "
    "installed on your computer and is configured with the same mpi options "
    "as deal.II\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "PETSc installed with --prefix=<...> to a destination:\n"
    "    $ PETSC_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DPETSC_DIR=\"...\" <...>\n"
    "PETSc compiled in source tree:\n"
    "    $ PETSC_DIR=\"...\"  PETSC_ARCH=\"...\" cmake <...>\n"
    "    $ cmake -DPETSC_DIR=\"...\" -DPETSC_ARCH=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n\n"
    )
endmacro()


configure_feature(PETSC)

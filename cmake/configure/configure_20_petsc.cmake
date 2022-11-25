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
# Configuration for the petsc library:
#

set(FEATURE_PETSC_AFTER MPI)


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
    # Petsc has to be configured with the same MPI configuration as
    # deal.II.
    #
    # petscconf.h should export PETSC_HAVE_MPIUNI 1 in case mpi support is
    # _NOT_ enabled.
    # So we check for this:
    #
    if( (PETSC_WITH_MPIUNI AND DEAL_II_WITH_MPI)
         OR
         (NOT PETSC_WITH_MPIUNI AND NOT DEAL_II_WITH_MPI))
      message(STATUS "Could not find a sufficient PETSc installation: "
        "PETSc has to be configured with the same MPI configuration as deal.II."
        )
      set(PETSC_ADDITIONAL_ERROR_STRING
        ${PETSC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient PETSc installation:\n"
        "PETSc has to be configured with the same MPI configuration as deal.II, but found:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
        "  PETSC_WITH_MPI   = (NOT ${PETSC_WITH_MPIUNI})\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Petsc has to be configured with the same number of bits for indices as
    # deal.II.
    #
    # petscconf.h should export PETSC_WITH_64BIT_INDICES 1 in case 64bits
    # indices support is enabled.
    # So we check for this:
    #
    if( (NOT PETSC_WITH_64BIT_INDICES AND DEAL_II_WITH_64BIT_INDICES)
         OR
         (PETSC_WITH_64BIT_INDICES AND NOT DEAL_II_WITH_64BIT_INDICES))
      message(STATUS "Could not find a sufficient PETSc installation: "
        "PETSc has to be configured to use the same number of bits for the "
        "global indices as deal.II."
        )
      set(PETSC_ADDITIONAL_ERROR_STRING
        ${PETSC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient PETSc installation:\n"
        "PETSc has to be configured to use the same number of bits for the "
        "global indices as deal.II, but found:\n"
        "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
        "  PETSC_WITH_64BIT_INDICES = (${PETSC_WITH_64BIT_INDICES})\n"
        )
      set(${var} FALSE)
    endif()

    # If PETSc is compiled with complex scalar type we need to have support
    # for complex values within deal.II as well.
    #
    if( PETSC_WITH_COMPLEX AND NOT DEAL_II_WITH_COMPLEX_VALUES )
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

    if(DEAL_II_PETSC_WITH_KOKKOS)
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
set(DEAL_II_PETSC_WITH_COMPLEX ${PETSC_WITH_COMPLEX})
set(DEAL_II_PETSC_WITH_HYPRE ${PETSC_WITH_HYPRE})
set(DEAL_II_PETSC_WITH_MUMPS ${PETSC_WITH_MUMPS})
set(DEAL_II_PETSC_WITH_KOKKOS ${PETSC_WITH_KOKKOS})

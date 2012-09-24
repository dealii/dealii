#
# Configuration for the petsc library:
#

MACRO(FEATURE_PETSC_FIND_EXTERNAL var)

  FIND_PACKAGE(PETSC)

  IF(PETSC_FOUND)
    #
    # So, we have found a petsc library. Let's check whether we can use it.
    #
    SET(${var} TRUE)

    #
    # We support petsc from version 3.x.x onwards
    #
    IF(PETSC_VERSION_MAJOR LESS 3)
      MESSAGE(STATUS
        "Could not find a sufficient modern petsc installation: "
        "Version >=3.0.0 required!"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Petsc has to be configured with the same MPI configuration as
    # deal.II.
    #
    # petscconf.h should export PETSC_HAVE_MPIUNI 1 in case  mpi support is
    # _NOT_ enabled.
    # So we check for this:
    #
    IF( (PETSC_WITH_MPIUNI AND DEAL_II_COMPILER_SUPPORTS_MPI)
         OR
         (NOT PETSC_WITH_MPIUNI AND NOT DEAL_II_COMPILER_SUPPORTS_MPI))
      MESSAGE(STATUS
        "Could not find a sufficient petsc installation: "
        "Petsc has to be configured with the same MPI configuration as deal.II."
        )
      SET(${var} FALSE)
    ENDIF()

  ENDIF()

ENDMACRO()


MACRO(FEATURE_PETSC_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})

  # The user has to know the location of the petsc headers as well: # TODO
  # LIST(APPEND DEAL_II_EXTERNAL_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS})

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    ${PETSC_LIBRARIES}
    )

  SET(DEAL_II_USE_PETSC TRUE)

  #
  # Work around a stupidity in PETSc that makes sure it interferes in
  # a completely obnoxious way with boost.
  # TODO: Obosolete?
  #
  ADD_DEFINITIONS(-DPETSC_SKIP_UNDERSCORE_CHKERR)

  #
  # Set some definitions for config.h:
  #

  IF(NOT PETSC_RELEASE)
    SET(DEAL_II_USE_PETSC_DEV TRUE)
  ENDIF()

  IF(PETSC_COMPLEX)
    SET(DEAL_II_USE_PETSC_COMPLEX TRUE)
  ENDIF()

  SET(DEAL_II_EXPAND_PETSC_VECTOR "PETScWrappers::Vector")
  SET(DEAL_II_EXPAND_PETSC_MPI_VECTOR "PETScWrappers::MPI::Vector")
  SET(DEAL_II_EXPAND_PETSC_BLOCKVECTOR "PETScWrappers::BlockVector")
  SET(DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR "PETScWrappers::MPI::BlockVector")

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(PETSC)


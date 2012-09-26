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
      MESSAGE(WARNING "\n"
        "Could not find a sufficient modern petsc installation: "
        "Version >=3.0.0 required!\n\n"
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
      MESSAGE(WARNING "\n"
        "Could not find a sufficient petsc installation: "
        "Petsc has to be configured with the same MPI configuration as deal.II.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()

  ENDIF()

ENDMACRO()


MACRO(FEATURE_PETSC_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})

  # The user has to know the location of the petsc headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS})

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    ${PETSC_LIBRARIES}
    )

  SET(DEAL_II_USE_PETSC TRUE)

  #
  # Disable a bunch of warnings when compiling with petsc:
  #
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-long-long")


  #
  # Work around a stupidity in PETSc that makes sure it interferes in
  # a completely obnoxious way with boost.
  # TODO: Obsolete?
  #
  #SET(PETSC_SKIP_UNDERSCORE_CHKERR TRUE)

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
  SET(DEAL_II_EXPAND_PETSC_BLOCKVECTOR "PETScWrappers::BlockVector")

  IF(DEAL_II_WITH_MPI)
    SET(DEAL_II_EXPAND_PETSC_MPI_VECTOR "PETScWrappers::MPI::Vector")
    SET(DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR "PETScWrappers::MPI::BlockVector")
  ENDIF()

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_PETSC_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_PETSC_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "\n"
    "Could not find the petsc library!\n\n"
    "Please ensure that the petsc library version 3.0.0 or newer is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "PETSc installed with --prefix=<...> to a destination:\n"
    "    $ PETSC_DIR=\"...\" cmake <...>\n"
    "    $ ccmake -DPETSC_DIR=\"...\" cmake <...>\n"
    "PETSc compiled in source tree:\n"
    "    $ PETSC_DIR=\"...\"  PETSC_ARCH=\"...\" cmake <...>\n"
    "    $ ccmake -DPETSC_DIR=\"...\" -DPETSC_ARCH=\"...\" cmake <...>\n"
    "or set the relevant variables by hand in ccmake.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(PETSC)


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
    IF(PETSC_MAJOR LESS 3)
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
    FILE(STRINGS "${PETSC_PETSCCONF_H}" PETSC_RELEASE_STRING
      REGEX "#define.*PETSC_HAVE_MPIUNI 1")
    IF("${PETSC_RELEASE_STRING}" STREQUAL "")
      SET(PETSC_WITH_MPIUNI FALSE)
    ELSE()
      SET(PETSC_WITH_MPIUNI TRUE)
    ENDIF()

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
  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(PETSC)


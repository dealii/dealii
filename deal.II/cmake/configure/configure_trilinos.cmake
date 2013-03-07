#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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
    # Check whether all required modules of trilinos are installed:
    #
    MESSAGE(STATUS
      "Check whether the found trilinos package contains all required modules:"
      )


    FOREACH(_module
      amesos epetra ifpack aztecoo sacado teuchos
      )
      LIST_CONTAINS(_module_found ${_module} ${Trilinos_LIBRARIES})
      IF(_module_found)
        MESSAGE(STATUS "Found ${_module}")
      ELSE()
        MESSAGE(STATUS "Module ${_module} not found!")
        SET(_modules_missing "${_modules_missing} ${_module}")
        SET(${var} FALSE)
      ENDIF()
    ENDFOREACH()

    IF(NOT ${var})
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        "The Trilinos installation found at\n"
        "  ${TRILINOS_DIR}\n"
        "is missing one or more modules necessary for the deal.II Trilinos interfaces:\n"
        "  ${_modules_missing}\n\n"
        "Please re-install Trilinos with the missing Trilinos subpackages
        enabled.\n\n"
        )
      MESSAGE(WARNING "\n" ${TRILINOS_ADDITIONAL_ERROR_STRING} "\n")
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

      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation found at\n"
        "  ${TRILINOS_DIR}\n"
        "with version ${TRILINOS_VERSION_MAJOR}.${TRILINOS_VERSION_MINOR}.${TRILINOS_VERSION_SUBMINOR} has bugs that make\n"
        "it incompatible with deal.II. Please use versions before 10.6 or after\n"
        "10.8.1.\n\n"
        )
      MESSAGE(WARNING "\n" ${TRILINOS_ADDITIONAL_ERROR_STRING} "\n")
      SET(${var} FALSE)
    ENDIF()

    #
    # Trilinos has to be configured with the same MPI configuration as
    # deal.II.
    #
    IF( (TRILINOS_WITH_MPI AND NOT DEAL_II_WITH_MPI)
         OR
         (NOT TRILINOS_WITH_MPI AND DEAL_II_WITH_MPI))
      SET(TRILINOS_ADDITIONAL_ERROR_STRING
        ${TRILINOS_ADDITIONAL_ERROR_STRING}
        "The Trilinos installation found at\n"
        "  ${TRILINOS_DIR}\n"
        "has to be configured with the same MPI configuration as deal.II, but found:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
        "  TRILINOS_WITH_MPI = ${TRILINOS_WITH_MPI}\n"
        )
      MESSAGE(WARNING "\n" ${TRILINOS_ADDITIONAL_ERROR_STRING} "\n")
      SET(${var} FALSE)
    ENDIF()

    #
    # Some verions of Sacado_cmath.hpp does things that aren't compatible
    # with the -std=c++0x flag of GCC, see deal.II FAQ.
    # Test whether that is indeed the case
    #
    SET(CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIRS})
    PUSH_TEST_FLAG("${DEAL_II_CXX11_FLAG}")

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <Sacado_cmath.hpp>
      int main(){ return 0; }
      "
      TRILINOS_SUPPORTS_CPP11)

    IF(DEAL_II_CAN_USE_CXX11 AND NOT TRILINOS_SUPPORTS_CPP11)
      #
      # Try whether exporting HAS_C99_TR1_CMATH helps:
      #
      PUSH_TEST_FLAG("-DHAS_C99_TR1_CMATH")
      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <Sacado_cmath.hpp>
        int main(){ return 0; }
        "
        TRILINOS_HAS_C99_TR1_WORKAROUND)
      POP_TEST_FLAG()

      IF(TRILINOS_HAS_C99_TR1_WORKAROUND)
        LIST(APPEND DEAL_II_DEFINITIONS "HAS_C99_TR1_CMATH")
        LIST(APPEND DEAL_II_USER_DEFINITIONS "HAS_C99_TR1_CMATH")
      ELSE()
        SET(TRILINOS_ADDITIONAL_ERROR_STRING
          ${TRILINOS_ADDITIONAL_ERROR_STRING}
          "The Trilinos installation found at\n"
          "  ${TRILINOS_DIR}\n"
          "is not compatible with the C++ standard selected for\n"
          "this compiler. See the deal.II FAQ page for a solution.\n\n"
          )
        MESSAGE(WARNING "\n" ${TRILINOS_ADDITIONAL_ERROR_STRING} "\n")
        SET(${var} FALSE)
      ENDIF()
    ENDIF()

    POP_TEST_FLAG()
    SET(CMAKE_REQUIRED_INCLUDES)

    #
    # Remove the following variables from the cache to force a recheck:
    #
    UNSET(TRILINOS_SUPPORTS_CPP11 CACHE)
    UNSET(TRILINOS_HAS_C99_TR1_WORKAROUND CACHE)

  ENDIF(TRILINOS_FOUND)
ENDMACRO()


MACRO(FEATURE_TRILINOS_CONFIGURE_EXTERNAL)
  INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIRS})

  # The user has to know the location of the trilinos headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${TRILINOS_INCLUDE_DIRS})


  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    # The Trilinos libraries:
    ${TRILINOS_LIBRARIES}
    # All external libraries necessary for the Trilinos libraries. Nice and
    # easy :-)
    ${Trilinos_TPL_LIBRARIES}
    )

  SET(DEAL_II_EXPAND_TRILINOS_VECTOR "TrilinosWrappers::Vector")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCKVECTOR "TrilinosWrappers::BlockVector")
  SET(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN "TrilinosWrappers::SparsityPattern")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN "TrilinosWrappers::BlockSparsityPattern")

  IF(DEAL_II_WITH_MPI)
    SET(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR "TrilinosWrappers::MPI::BlockVector")
    SET(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR "TrilinosWrappers::MPI::Vector")
  ENDIF()

  #
  #  used with -W -Wall (which includes -Wunused). Regrettable
  #  though it may be, these warnings pretty much drown everything
  #  else and we better disable some of the warnings to enable us
  #  to see through the clutter.
  #
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-extra")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-overloaded-virtual")
ENDMACRO()


CONFIGURE_FEATURE(TRILINOS)

#
# Configuration for the trilinos library:
#

#
# TODO: Fix up WARNING/STATUS of the failure messages.
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

    SET(macro_modules_list amesos epetra ifpack aztecoo sacado teuchos)

    FOREACH(macro_module ${macro_modules_list})
      LIST_CONTAINS(macro_module_found ${macro_module} ${Trilinos_LIBRARIES})
      IF(macro_module_found)
        MESSAGE(STATUS "Found ${macro_module}")
      ELSE()
        MESSAGE(STATUS "Module ${macro_module} not found!")
        SET(macro_modules_missing "${macro_modules_missing} ${macro_module}")
        SET(${var} FALSE)
      ENDIF()
    ENDFOREACH()

    IF(NOT ${var})
      MESSAGE(WARNING "\n"
        "The Trilinos installation is missing one or more modules necessary for\n"
        "the deal.II Trilinos interfaces:\n"
        "${macro_modules_missing}\n\n"
        "Please re-install Trilinos with the missing Trilinos subpackages enabled.\n\n"
        )
    ENDIF()


    #
    # Trilinos 10.6 had quite a number of bugs we ran into, see
    # for example
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5062
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5319
    #
    IF(TRILINOS_VERSION_MAJOR EQUAL 10 AND TRILINOS_VERSION_MINOR EQUAL 6)
      MESSAGE(WARNING "\n"
        "Trilinos versions ${TRILINOS_VERSION_MAJOR}.${TRILINOS_VERSION_MINOR}.x have bugs that make\n"
        "it incompatible with deal.II. Please use versions before 10.6 or after\n"
        "10.8.1.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()


    #
    # The same is unfortunately true for 10.8.[01]:
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5370
    #
    IF( TRILINOS_VERSION_MAJOR EQUAL 10 AND
        TRILINOS_VERSION_MINOR EQUAL 8 AND
        TRILINOS_VERSION_SUBMINOR LESS 2 )
      MESSAGE(WARNING "\n"
        "Trilinos versions 10.8.0 and 10.8.1 have bugs that make\n"
        "it incompatible with deal.II. Please use versions before 10.6 or after\n"
        "10.8.1.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()


    #
    # Trilinos has to be configured with the same MPI configuration as
    # deal.II.
    #
    # Epetra installs Epetra_MpiComm.h if configured trilinos was
    # configured with mpi. We use this as a check for the mpi configuration
    # of Epetra.
    #
    IF(EXISTS "${TRILINOS_INCLUDE_DIRS}/Epetra_MpiComm.h")
      SET(TRILINOS_WITH_MPI TRUE)
    ENDIF()

    IF( (TRILINOS_WITH_MPI AND NOT DEAL_II_COMPILER_SUPPORTS_MPI)
         OR
         (NOT TRILINOS_WITH_MPI AND DEAL_II_COMPILER_SUPPORTS_MPI))
      MESSAGE(WARNING "\n"
        "Trilinos has to be configured with the same MPI configuration as deal.II.\n\n"
        )
      SET(${var} FALSE)
    ENDIF()


    #
    # Some verions of Sacado_cmath.hpp does things that aren't compatible
    # with the -std=c++0x flag of GCC, see deal.II FAQ.
    # Test whether that is indeed the case
    #
    LIST(APPEND CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIRS})
    ADD_FLAGS(CMAKE_REQUIRED_FLAGS "-std=c++0x")
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
      ADD_FLAGS(CMAKE_REQUIRED_FLAGS "-DHAS_C99_TR1_CMATH")
      CHECK_CXX_SOURCE_COMPILES(
        "
        #include <Sacado_cmath.hpp>
        int main(){ return 0; }
        "
        TRILINOS_HAS_C99_TR1_WORKAROUND)
      STRIP_FLAG(CMAKE_REQUIRED_FLAGS "-DHAS_C99_TR1_CMATH")

      IF(TRILINOS_HAS_C99_TR1_WORKAROUND)
        LIST(APPEND DEAL_II_DEFINITIONS "HAS_C99_TR1_CMATH")
        LIST(APPEND DEAL_II_USER_DEFINITIONS "HAS_C99_TR1_CMATH")
      ELSE()
        MESSAGE(WARNING "\n"
          "Your Trilinos installation is not compatible with the C++ standard selected for\n"
          "this compiler. See the deal.II FAQ page for a solution.\n\n"
          )
        SET(${var} FALSE)
      ENDIF()
    ENDIF()
    STRIP_FLAG(CMAKE_REQUIRED_FLAGS "-std=c++0x")
    LIST(REMOVE_ITEM CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIRS})

  ENDIF(TRILINOS_FOUND)
ENDMACRO()


MACRO(FEATURE_TRILINOS_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIRS})

  # The user has to know the location of the trilinos headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${TRILINOS_INCLUDE_DIRS})


  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    ${TRILINOS_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
    )

  SET(DEAL_II_USE_TRILINOS TRUE)

  SET(DEAL_II_EXPAND_TRILINOS_VECTOR "TrilinosWrappers::Vector")
  SET(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR "TrilinosWrappers::MPI::Vector")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCKVECTOR "TrilinosWrappers::BlockVector")
  SET(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR "TrilinosWrappers::MPI::BlockVector")
  SET(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN "TrilinosWrappers::SparsityPattern")
  SET(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN "TrilinosWrappers::BlockSparsityPattern")


  #
  #  used with -W -Wall (which includes -Wunused). Regrettable
  #  though it may be, these warnings pretty much drown everything
  #  else and we better disable some of the warnings to enable us
  #  to see through the clutter.
  #
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-extra")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-overloaded-virtual")

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(TRILINOS)


#
# Configuration for the trilinos library:
#

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckIncludeFile)


MACRO(FEATURE_TRILINOS_FIND_EXTERNAL var)

  FIND_PACKAGE(TRILINOS)

  IF(TRILINOS_FOUND)
    SET(${var} TRUE)
    #
    # So, we have a library. Let's see whether we can use it:
    #


    #
    # Trilinos 10.6 had quite a number of bugs we ran into, see
    # for example
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5062
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5319
    #
    IF(TRILINOS_MAJOR EQUAL 10 AND TRILINOS_MINOR EQUAL 6)
      MESSAGE(WARNING "
Trilinos versions ${TRILINOS_MAJOR}.${TRILINOS_MINOR}.x have bugs that make
it incompatible with deal.II. Please use versions before 10.6 or after
10.8.1.

")
      SET(${var} FALSE)
    ENDIF()


    #
    # The same is unfortunately true for 10.8.[01]:
    #   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5370
    #
    IF( TRILINOS_MAJOR EQUAL 10 AND
        TRILINOS_MINOR EQUAL 8 AND
        TRILINOS_SUBMINOR LESS 2 )
      MESSAGE(WARNING "
Trilinos versions 10.8.0 and 10.8.1 have bugs that make
it incompatible with deal.II. Please use versions before 10.6 or after
10.8.1.

")
      SET(${var} FALSE)
    ENDIF()


    #
    # Trilinos has to be configured with the same MPI configuration as
    # deal.II. So check this:
    #
    # TODO: Refine this check...
    #
    IF("${TRILINOS_MPI_LIBRARIES}" EQUAL "")
      SET(TRILINOS_USE_MPI TRUE)
    ENDIF()
    IF( (TRILINOS_USE_MPI AND NOT DEAL_II_COMPILER_SUPPORTS_MPI) OR
        (NOT TRILINOS_USE_MPI AND DEAL_II_COMPILER_SUPPORTS_MPI))
      MESSAGE(WARNING "
Trilinos has to be configured with the same MPI configuration as deal.II.

")
      SET(${var} FALSE)
    ENDIF()


    #
    # Some verions of Sacado_cmath.hpp does things that aren't compatible
    # with the -std=c++0x flag of GCC, see deal.II FAQ.
    # Test whether that is indeed the case
    #
    LIST(APPEND CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIR})
    LIST(APPEND CMAKE_REQUIRED_FLAGS "-std=c++0x")
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <Sacado_cmath.hpp>
      int main(){ return 0; }
      "
      TRILINOS_SUPPORTS_CPP11)
    LIST(REMOVE_ITEM CMAKE_REQUIRED_FLAGS "-std=c++0x")
    LIST(REMOVE_ITEM CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIR})

    IF(DEAL_II_CAN_USE_CXX11 AND NOT TRILINOS_SUPPORTS_CPP11)
      MESSAGE(WARNING "
Your Trilinos installation is not compatible with the C++ standard selected for this compiler.
See the deal.II FAQ page for a solution.

")
      SET(${var} FALSE)
    ENDIF()


    #
    # Check whether we can find all header files used in deal.II:
    #
    LIST(APPEND CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIR})
    CHECK_INCLUDE_FILE("Amesos.h" TRILINOS_HAVE_AMESOS_H)
    CHECK_INCLUDE_FILE("Epetra_CrsGraph.h" TRILINOS_HAVE_EPETRA_CRSGRAPH_H)
    CHECK_INCLUDE_FILE("Epetra_CrsMatrix.h" TRILINOS_HAVE_EPETRA_CRSMATRIX_H)
    CHECK_INCLUDE_FILE("Epetra_Import.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Epetra_LinearProblem.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Epetra_Map.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Epetra_MultiVector.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Epetra_Operator.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Epetra_SerialComm.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Epetra_Vector.h" TRILINOS_HAVE_EPETRA_H)
    CHECK_INCLUDE_FILE("Ifpack.h" TRILINOS_HAVE_IFPACK_H)
    CHECK_INCLUDE_FILE("ml_MultiLevelPreconditioner.h" TRILINOS_HAVE_ML_MULTILEVELPRECONDITIONER_H)
    CHECK_INCLUDE_FILE("AztecOO.h" TRILINOS_HAVE_AZTECOO_H)
    CHECK_INCLUDE_FILE("AztecOO_Operator.h" TRILINOS_HAVE_AZTECOO_OPERATOR_H)
    CHECK_INCLUDE_FILE("Sacado.hpp" TRILINOS_HAVE_SACADO_HPP)
    CHECK_INCLUDE_FILE("Teuchos_ParameterList.hpp" TRILINOS_HAVE_TEUCHOS_PARAMETERLIST_HPP)
    CHECK_INCLUDE_FILE("Teuchos_RCP.hpp" TRILINOS_HAVE_TEUCHOS_RCP_HPP)
    CHECK_INCLUDE_FILE("Teuchos_RefCountPtr.hpp" TRILINOS_HAVE_TEUCHOS_REFCOUNTPTR_HPP)
    LIST(REMOVE_ITEM CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIR})

    #
    # TODO:
    #
    IF(NOT
        (TRILINOS_HAVE_AMESOS_H AND
        TRILINOS_HAVE_EPETRA_CRSGRAPH_H AND
        TRILINOS_HAVE_EPETRA_CRSMATRIX_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_EPETRA_H AND
        TRILINOS_HAVE_IFPACK_H AND
        TRILINOS_HAVE_ML_MULTILEVELPRECONDITIONER_H AND
        TRILINOS_HAVE_AZTECOO_H AND
        TRILINOS_HAVE_AZTECOO_OPERATOR_H AND
        TRILINOS_HAVE_SACADO_HPP AND
        TRILINOS_HAVE_TEUCHOS_PARAMETERLIST_HPP AND
        TRILINOS_HAVE_TEUCHOS_RCP_HPP AND
        TRILINOS_HAVE_TEUCHOS_REFCOUNTPTR_HPP))
#      MESSAGE(WARNING "
#The Trilinos installation is missing one or more header files necessary for
#the deal.II Trilinos interfaces. Please re-install Trilinos with the missing
#Trilinos subpackages enabled.
#
#")
#     SET(${var} FALSE)
    ENDIF()

  ENDIF(TRILINOS_FOUND)

ENDMACRO()


MACRO(FEATURE_TRILINOS_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})

  LIST(APPEND deal_ii_external_libraries
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
  CHECK_CXX_COMPILER_FLAG(
    "-Wno-unused"
    DEAL_II_HAVE_WNO_UNUSED_FLAG
    )
  IF(DEAL_II_HAVE_WNO_UNUSED_FLAG)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused")
  ENDIF()

  CHECK_CXX_COMPILER_FLAG(
    "-Wno-overloaded-virtual"
    DEAL_II_HAVE_WNO_OVERLOADED_VIRTUAL_FLAG
    )
  IF(DEAL_II_HAVE_WNO_OVERLOADED_VIRTUAL_FLAG)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overloaded-virtual")
  ENDIF()

  CHECK_CXX_COMPILER_FLAG(
    "-Wno-extra"
    DEAL_II_HAVE_WNO_EXTRA_FLAG
    )
  IF(DEAL_II_HAVE_WNO_EXTRA_FLAG)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-extra")
  ENDIF()


  SET(${var} TRUE)

ENDMACRO()


CONFIGURE_FEATURE(TRILINOS)

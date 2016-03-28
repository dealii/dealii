CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)
INCLUDE(../setup_testsubproject.cmake)
PROJECT(testsuite CXX)
# FIXME: lac/matrix_out_03 fails with complex-valued PETSc
IF(NOT DEAL_II_PETSC_WITH_COMPLEX)
  DEAL_II_PICKUP_TESTS()
ENDIF()

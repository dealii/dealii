## ---------------------------------------------------------------------
##
## Copyright (C) 2020 - 2021 by the deal.II authors
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
# Make sure that the tests are picked up by a global CTest call
#
file(APPEND ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "subdirs(contrib/python-bindings/tests)\n")
enable_testing()

file(GLOB _tests "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
foreach(_test_path ${_tests})
 get_filename_component(_test ${_test_path} NAME_WE)
 get_filename_component(_test_directory ${_test_path} DIRECTORY)
 add_test(NAME python-bindings/${_test} COMMAND ${PYTHON_EXECUTABLE} ${_test_path})
 set_tests_properties(python-bindings/${_test} PROPERTIES ENVIRONMENT
   "PYTHONPATH=${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}/../:$ENV{PYTHONPATH};DEAL_II_PYTHON_TESTPATH=${_test_directory}"
   )
endforeach()

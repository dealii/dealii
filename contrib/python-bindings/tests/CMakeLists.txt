## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2020 - 2022 by the deal.II authors
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

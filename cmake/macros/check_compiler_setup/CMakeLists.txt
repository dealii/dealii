## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2023 - 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

# use a "Custom" build type to avoid auto populated compiler flags
set(CMAKE_BUILD_TYPE "Custom")
set(CMAKE_CONFIGURATION_TYPES "DEBUG;RELEASE")
cmake_minimum_required(VERSION 3.16.0)
project(CheckCompilerSetup)
add_executable(CheckCompilerSetupExec dummy.cpp)

target_compile_options(CheckCompilerSetupExec PRIVATE ${TEST_COMPILE_OPTIONS})
target_link_options(CheckCompilerSetupExec PRIVATE ${TEST_LINK_OPTIONS})
target_link_libraries(CheckCompilerSetupExec ${TEST_LINK_LIBRARIES})

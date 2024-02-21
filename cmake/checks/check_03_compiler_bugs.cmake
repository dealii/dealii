## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

########################################################################
#                                                                      #
#                   Check for various compiler bugs:                   #
#                                                                      #
########################################################################


#
# In intel (at least 13.1 and 14), vectorization causes
# wrong code. See https://code.google.com/p/dealii/issues/detail?id=156
# or tests/hp/solution_transfer.cc
# A work-around is to disable all vectorization.
#
# - Timo Heister, 2013, 2015
#
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.3" )
  enable_if_supported(DEAL_II_CXX_FLAGS_RELEASE "-no-vec")
endif()


#
# Intel 16.0.1 produces wrong code that creates a race condition in
# tests/fe/curl_curl_01.debug but 16.0.2 is known to work. Blacklist this
# version. Also see github.com/dealii/dealii/issues/2203
#
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "16.0.1" )
  message(FATAL_ERROR "Intel compiler version 16.0.1 is not supported, please update to 16.0.2 or newer!")
endif()


#
# Check for a regression in gcc-11.1.0 where a deleted move constructor
# prevents templated constructor from being used. For details see
#
#   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100644
#   https://github.com/dealii/dealii/issues/12244
#   https://github.com/dealii/dealii/pull/12246
#
# - Mathias Anselmann, Matthias Maier, David Wells, 2021
#
check_cxx_compiler_bug(
  "
  struct NonMovable {
    NonMovable() = default;
    NonMovable(NonMovable &&) = delete;
  };
  template <class T> struct Maybe {
    NonMovable mMember;
    template <typename U> Maybe(Maybe<U> &&) : mMember() {}
  };
  void unlucky(Maybe<int> &&x) { Maybe<int> var{(Maybe<int> &&) x}; }
  int main() { return 0; }
  "
  DEAL_II_DELETED_MOVE_CONSTRUCTOR_BUG)

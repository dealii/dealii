## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

########################################################################
#                                                                      #
#                   Check for various compiler bugs:                   #
#                                                                      #
########################################################################


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

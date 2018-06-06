## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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

########################################################################
#                                                                      #
#                   Check for various compiler bugs:                   #
#                                                                      #
########################################################################

#
# In some cases, we would like to name partial specializations
# as friends. However, the standard forbids us to do so. But
# then, we can declare the general template as a friend, and
# at least gcc extends the friendship to all specializations
# of the templates, which is not what the standard says.
#
# With other compilers, most notably cxx, this does not work.
# In this case, we can make individual specializations friends,
# which in turn gcc rejects. So check, whether this is possible.
#
# The respective clause in the standard is 14.5.3.1, which gives
# this example:
#   template<class T> class task {
#     friend class task<int>;
#   };
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  template <int N, typename T> class X;
  template <typename T>        class X<1,T>;

  template <typename P> class Y {
      static int i;
      template <int N, typename T> friend class X;
      friend class X<1,P>;
  };

  template <typename T> class X<1,T> {
      int f () { return Y<T>::i; };     // access private field
  };
  int main() { return 0; }
  "
  DEAL_II_TEMPL_SPEC_FRIEND_BUG)


#
# This is a variant of the previous test. Some icc 11.0
# builds (sub-releases) on Windows apparently don't allow
# the declaration of an explicit specialization of member
# arrays of templates:
#
# template <int dim>
# struct X
# {
#    static const int N = 2*dim;
#    static const int x[N];
# };
# template <> const int X<2>::x[N];
#
# That version of icc requests that there be an initialization,
# i.e. it thinks that this is the *definition*, not merely a
# *declaration* of an explicit specialization. This is wrong,
# however.
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  template <int dim>
  struct X
  {
    static const int N = 2*dim;
    static const int x[N];
  };
  template <> const int X<2>::x[N];
  int main() { return 0; }
  "
  DEAL_II_MEMBER_ARRAY_SPECIALIZATION_BUG
  )


#
# Some older versions of gcc deduce pointers to const functions in
# template contexts to pointer-to-function of const objects.
# This is not correct
#
# Check for this misfeature.
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  template <typename T> struct identity { typedef T type; };
  template <typename C> void new_thread (void (C::*fun_ptr)(),
                typename identity<C>::type &c) {}
  template <typename C> void new_thread (void (C::*fun_ptr)() const,
                const typename identity<C>::type &c) {}
  struct X { void f() const{} };

  int main()
  {
    X x;
    new_thread (&X::f, x);
  }
  "
  DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  )



#
# Microsoft Visual C++ has a bug where the resulting object
# from calling std::bind does not have a const operator(),
# so we cannot pass such objects as const references as we
# usually do with input arguments of other functions.
#
# - Wolfgang Bangerth, 2014
#
PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_FLAGS}")
CHECK_CXX_COMPILER_BUG(
  "
  #include <functional>

  void f(int, int) {}

  template <typename F>
  void g(const F &func)
  {
    func(1);
  }

  int main ()
  {
    g (std::bind(&f, std::placeholders::_1, 1));
  }
  "
  DEAL_II_BIND_NO_CONST_OP_PARENTHESES
  )
RESET_CMAKE_REQUIRED()


#
# In intel (at least 13.1 and 14), vectorization causes
# wrong code. See https://code.google.com/p/dealii/issues/detail?id=156
# or tests/hp/solution_transfer.cc
# A work-around is to disable all vectorization.
#
# - Timo Heister, 2013, 2015
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.3" )
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-no-vec")
ENDIF()


#
# gcc-4.8.1 has some problems with the constexpr "vertices_per_cell" in the
# definition of alternating_form_at_vertices.
#
# TODO: Write a unit test.
#
# For now, just enable the workaround for Windows targets
#
# - Matthias Maier, 2013
#
IF( CMAKE_SYSTEM_NAME MATCHES "CYGWIN"
    OR CMAKE_SYSTEM_NAME MATCHES "Windows" )
  SET(DEAL_II_CONSTEXPR_BUG TRUE)
ENDIF()


#
# Intel 16.0.1 produces wrong code that creates a race condition in
# tests/fe/curl_curl_01.debug but 16.0.2 is known to work. Blacklist this
# version. Also see github.com/dealii/dealii/issues/2203
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "16.0.1" )
  MESSAGE(FATAL_ERROR "Intel compiler version 16.0.1 is not supported, please update to 16.0.2 or newer!")
ENDIF()

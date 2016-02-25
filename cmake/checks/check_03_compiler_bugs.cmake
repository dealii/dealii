## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

########################################################################
#                                                                      #
#                   Check for various compiler bugs:                   #
#                                                                      #
########################################################################


#
# On some gcc 4.3 snapshots, a 'const' qualifier on a return type triggers a
# warning. This is unfortunate, since we happen to stumble on this
# in some of our template trickery with iterator classes. If necessary,
# do not use the relevant warning flag
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
PUSH_CMAKE_REQUIRED("-Wreturn-type")
PUSH_CMAKE_REQUIRED("-Werror")
CHECK_CXX_COMPILER_BUG(
  "
  const double foo() { return 1.; }
  int main() { return 0; }
  "
  DEAL_II_WRETURN_TYPE_CONST_QUALIFIER_BUG
  )
RESET_CMAKE_REQUIRED()

IF(DEAL_II_WRETURN_TYPE_CONST_QUALIFIER_BUG)
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-return-type)
ENDIF()


#
# gcc 4.4 has an interesting problem in that it doesn't
# care for one of BOOST signals2's header files and produces
# dozens of pages of error messages of the form
#   warning: invoking macro BOOST_PP_CAT argument 1: \
#   empty macro arguments are undefined in ISO C90 and ISO C++98
# This can be avoided by not using -pedantic for this compiler.
# For all other versions, we use this flag, however.
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
   CMAKE_CXX_COMPILER_VERSION MATCHES "4.4.")
  STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")
ENDIF()


#
# Newer gcc versions generate a large number of warnings inside boost if we
# are compiling without cxx11 but with -pedantic and there is no way to
# silence them.
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
   NOT DEAL_II_WITH_CXX11)
  STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")
ENDIF()


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
# Many compilers get this wrong (see Section 14.7.3.1, number (4)):
#
#   template <int dim> struct T {
#     static const int i;
#   };
#
#   template <> const int T<1>::i;
#   template <> const int T<1>::i = 1;
#
# First, by Section 14.7.3.14 of the standard, the first template<>
# line must necessarily be the _declaration_ of a specialization,
# and the second is then its definition. There is therefore no
# reason to report a doubly defined variable (Intel ICC 6.0), or
# to choke on these lines at all (Sun Forte)
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  template <int dim> struct T
  {
    static const int i;
  };
  template <> const int T<1>::i;
  template <> const int T<1>::i = 1;
  int main() {return 0;}
  "
  DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
  )


#
# Some older versions of gcc compile this, despite the 'explicit'
# keyword:
#
# struct X {
#     template <typename T>
#     explicit X(T);
# };
# void f(X);
# int main () { f(1); }
#
# Check for this misfeature.
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_SOURCE_COMPILES(
  "
  struct X {
    template <typename T>
    explicit X(T) {}
  };
  void f(X) {}
  int main() { f(1); }
  "
  DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
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
# Check for GCC bug 36052, see
#   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=36052
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#

CHECK_CXX_COMPILER_BUG(
  "
  struct S {
      typedef double value_type;
  };

  template <typename T> struct Traits {
      typedef const typename T::value_type dereference_type;
  };

  template <class BlockVectorType> struct ConstIterator {
      typedef typename Traits<BlockVectorType>::dereference_type dereference_type;

      dereference_type operator * () const  { return 0; }
  };
  template class ConstIterator<S>;
  int main(){return 0;}
  "
  DEAL_II_TYPE_QUALIFIER_BUG)

IF(DEAL_II_TYPE_QUALIFIER_BUG)
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS -Wno-ignored-qualifiers)
ENDIF()


#
# On Mac OS X, gcc appears to have a bug that prevents us from
# compiling a bit of code that involves boost::bind. Check for
# that.
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
IF(DEAL_II_HAVE_BUNDLED_DIRECTORY)
  CHECK_CXX_COMPILER_BUG(
    "
    #include <complex>
    #include <iostream>
    #include \"${BOOST_FOLDER}/include/boost/bind.hpp\"

    template<typename number>
    void bug_function (number test)
    {
      std::cout << test << std::endl;
    }
    int main()
    {
      std::complex<float> float_val (1., 2.);
      boost::bind(&bug_function<std::complex<float> >,
                  float_val)();
      return 0;
    }
    "
    DEAL_II_BOOST_BIND_COMPILER_BUG
    )
ENDIF()


#
# Microsoft Visual C++ has a bug where the resulting object
# from calling std::bind does not have a const operator(),
# so we cannot pass such objects as const references as we
# usually do with input arguments of other functions.
#
# - Wolfgang Bangerth, 2014
#
IF(DEAL_II_WITH_CXX11)
  PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
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
ELSE()
  CHECK_CXX_COMPILER_BUG(
    "
    #include <functional>
    #include \"${BOOST_FOLDER}/include/boost/bind.hpp\"

    void f(int, int) {}

    template <typename F>
    void g(const F &func)
    {
      func(1);
    }

    int main ()
    {
      using boost::bind;
      using boost::reference_wrapper;
  
      // now also import the _1, _2 placeholders from the global namespace
      // into the current one as suggested above
      using ::_1;

      g (boost::bind(&f, boost::_1, 1));
    }
    "
    DEAL_II_BIND_NO_CONST_OP_PARENTHESES
    )
ENDIF()


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
# Intel (at least 14, 15) has a bug where it incorrectly detects multiple
# matching function candidates and dies during type resolution in a
# perfectly valid SFINAE scenario. This seems to happen because the templated
# variant is not discarded (where it should be):
#
# error: more than one instance of overloaded function
#     "has_vmult_add<Range, T>::test [with Range=double, T=MyMatrix]"
# matches the argument list:
#     function template "void has_vmult_add<Range, T>::test<C>(decltype((<expression>))) [with Range=double, T=MyMatrix]"
#     function template "void has_vmult_add<Range, T>::test<C>(decltype((&C::vmult_add<double>))) [with Range=double, T=MyMatrix]"
# [...]
#
# - Matthias Maier, 2015
#

IF(DEAL_II_WITH_CXX11)
  PUSH_CMAKE_REQUIRED("${DEAL_II_CXX_VERSION_FLAG}")
  CHECK_CXX_COMPILER_BUG(
    "
    template <typename Range, typename T> struct has_vmult_add
    {
      template <typename C>
      static void test(decltype(&C::vmult_add));

      template <typename C>
      static void test(decltype(&C::template vmult_add<Range>));

      typedef decltype(test<T>(0)) type;
    };

    struct MyMatrix
    {
      void vmult_add() const;
    };

    int main()
    {
      typedef has_vmult_add<double, MyMatrix>::type test;
    }
    "
    DEAL_II_ICC_SFINAE_BUG
    )
  RESET_CMAKE_REQUIRED()
ENDIF()

#
# Intel 16.0.1 produces wrong code that creates a race condition in
# tests/fe/curl_curl_01.debug but 16.0.2 is known to work. Blacklist this
# version. Also see github.com/dealii/dealii/issues/2203
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "16.0.1" )
  MESSAGE(FATAL_ERROR "Intel compiler version 16.0.1 is not supported, please update to 16.0.2 or newer!")
ENDIF()

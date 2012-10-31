#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
#####

#
# Check for various compiler bugs:
#


#
# Some compiler versions, notably ICC, have trouble with the
# following code in which we explicitly call a destructor.
# This has to be worked around with a typedef. The problem is
# that the workaround fails with some other compilers, so that
# we can not unconditionally use the workaround...
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  namespace dealii
  {
    namespace FEValuesViews
    {
      template <int dim, int spacedim> struct Scalar {};
    }

    template <int dim, int spacedim>
    struct X
    {
        FEValuesViews::Scalar<dim,spacedim> scalars[dim*spacedim];

        void f()
          {
            scalars[0].dealii::FEValuesViews::Scalar<dim,spacedim>::~Scalar ();
          }
    };

    template struct X<2,2>;
  }
  int main() { return 0; }
  "
  DEAL_II_EXPLICIT_DESTRUCTOR_BUG
  )


#
# On some gcc 4.3 snapshots, a 'const' qualifier on a return type triggers a
# warning. This is unfortunate, since we happen to stumble on this
# in some of our template trickery with iterator classes. If necessary,
# do not use the relevant warning flag
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
PUSH_TEST_FLAG("-Wreturn-type")
PUSH_TEST_FLAG("-Werror")
CHECK_CXX_COMPILER_BUG(
  "
  const double foo() { return 1.; }
  int main() { return 0; }
  "
  DEAL_II_WRETURN_TYPE_CONST_QUALIFIER_BUG
  )
POP_TEST_FLAG()
POP_TEST_FLAG()

IF(DEAL_II_WRETURN_TYPE_CONST_QUALIFIER_BUG)
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS -Wno-return-type)
  ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS -Wno-return-type)
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
  STRIP_FLAG(CMAKE_CXX_FLAGS "-pedantic")
  STRIP_FLAG(CMAKE_C_FLAGS "-pedantic")
ENDIF()


#
# Some gcc compiler versions have a problem when using an unsigned count
# in the std::advance function. Unfortunately, this also happens
# occasionally from within the standard library, so we can't prevent the
# warning messages. Since this is annoying, switch of the flag -W which
# causes this.
#
# - Matthias Maier, rewritten 2012
#

# TODO: We use the mpi.h header file for this check. We should test this
# bug with another header file than mpi.h ...
CHECK_INCLUDE_FILE_CXX("mpi.h" HAVE_MPI_H)

IF(HAVE_MPI_H)
  PUSH_TEST_FLAG("-Wunused-parameter")
  PUSH_TEST_FLAG("-Werror")
  CHECK_CXX_COMPILER_BUG(
    "
    #include <mpi.h>
    int main() { return 0; }
    "
    DEAL_II_ADVANCE_WARNING_BUG)
  POP_TEST_FLAG()
  POP_TEST_FLAG()

  IF(DEAL_II_ADVANCE_WARNING_BUG)
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused-parameter")
    ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS "-Wno-unused-parameter")
  ENDIF()
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
# Some versions of gcc get this example wrong:
#
# struct X
# {
#     template <typename T> void operator << (T);
# };
# int main()
# {
#   X x;
#   x.operator << <double> (1);
# }
#
# They want to see a "template" for disambiguation in
#    x.template operator << <double> (1);
# which shouldn't be necessary since the left hand side of the
# dot operator is not template dependent. Surprisingly, this is
# only the case for operators, not if operator<< were a regular
# function. Annoyingly, other compilers barf on seeing the
# disambiguating "template" keyword.
#
# - Wolfgang Bangerth, Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  struct X
  {
      template <typename T> void operator << (T) {}
  };
  int main()
  {
    X x;
    x.operator << <double> (1);
    return 0;
  }
  "
  DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG
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
  int main() { f(1) };
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
  DEAL_II_TYPE_QUALIFIER_BUG
  )
IF(DEAL_II_TYPE_QUALIFIER_BUG)
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS -Wno-ignored-qualifiers)
  ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS -Wno-ignored-qualifiers)
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


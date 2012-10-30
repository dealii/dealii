dnl -------------------------------------------------------------
dnl
dnl Check for a problem with dynamic cast on dynamic libs on
dnl Mac OS X Snow Leopard. The test is only called on Mac OS X.
dnl The test is due to Scott Miller <scott.miller@psu.edu>
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_DYNAMIC_CAST_BUG,
[
  AC_MSG_CHECKING(for dynamic_cast problem with shared libs)

  if (cd contrib/config/tests/darwin-dynamic-cast ; \
      echo $CXX > compile_line; \
      $CXX -dynamiclib BaseClass.cpp -o libDynamicCastTestLib.dylib ; \
      $CXX -L. -lDynamicCastTestLib main.cc -o main ; \
      ./main) ; then
    AC_MSG_RESULT(no) ;
  else
    AC_MSG_RESULT(yes) ;
    AC_DEFINE(DEAL_II_HAVE_DARWIN_DYNACAST_BUG, 1,
              [Defined if the compiler has a bug with dynamic casting
               and dynamic libraries])
    if(test "`sw_vers -productVersion`" != "10.8");then
	CXXFLAGSG="$CXXFLAGSG -mmacosx-version-min=10.4"
	CXXFLAGSO="$CXXFLAGSO -mmacosx-version-min=10.4"
	LDFLAGS="$LDFLAGS -mmacosx-version-min=10.4"
    fi
  fi
  rm -f contrib/config/tests/darwin-dynamic-cast/libDynamicCastTestLib.dylib
  rm -f contrib/config/tests/darwin-dynamic-cast/main.o
  rm -f contrib/config/tests/darwin-dynamic-cast/main
])

dnl -------------------------------------------------------------
dnl In some cases, we would like to name partial specializations
dnl as friends. However, the standard forbids us to do so. But
dnl then, we can declare the general template as a friend, and
dnl at least gcc extends the friendship to all specializations
dnl of the templates, which is not what the standard says.
dnl
dnl With other compilers, most notably cxx, this does not work.
dnl In this case, we can make individual specializations friends,
dnl which in turn gcc rejects. So check, whether this is possible.
dnl
dnl The respective clause in the standard is 14.5.3.1, which gives
dnl this example:
dnl   template<class T> class task {
dnl     friend class task<int>;
dnl   };
dnl
dnl
dnl Usage: DEAL_II_CHECK_TEMPL_SPEC_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPL_SPEC_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for template specialization friend bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
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

template class X<1,int>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE_UNQUOTED(DEAL_II_TEMPL_SPEC_FRIEND_BUG, 1,
                         [Define if we have to work around a bug with some
                          compilers that will not allow us to specify a
                          fully specialized class of a template as a friend.
                          See the aclocal.m4 file in the top-level directory
                          for a description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl This is a variant of the previous test. Some icc 11.0
dnl builds (sub-releases) on Windows apparently don't allow
dnl the declaration of an explicit specialization of member
dnl arrays of templates:
dnl ---------------------------------
dnl template <int dim>
dnl struct X
dnl {
dnl    static const int N = 2*dim;
dnl    static const int x[N];
dnl };
dnl template <> const int X<2>::x[N];
dnl ---------------------------------
dnl That version of icc requests that there be an initialization,
dnl i.e. it thinks that this is the *definition*, not merely a
dnl *declaration* of an explicit specialization. This is wrong,
dnl however.
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_ARRAY_SPECIALIZATION_SPEC_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_ARRAY_SPECIALIZATION_BUG, dnl
[
  AC_MSG_CHECKING(for member array specialization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <int dim>
      struct X
      {
        static const int N = 2*dim;
        static const int x[N];
      };

      template <> const int X<2>::x[N];
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_MEMBER_ARRAY_SPECIALIZATION_BUG, 1,
                     [Defined if the compiler refuses to allow the
                      explicit specialization of static member
                      arrays. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
    ])
])


dnl -------------------------------------------------------------
dnl We compile all the files in the deal.II subdirectory multiple
dnl times, for the various space dimensions. When a program
dnl needs the libs for more than one dimension, the symbols in
dnl these libs must be either static or weak. Unfortunately, some
dnl compilers don't mark functions in anonymous namespaces as
dnl static, and also mangle their names the same way each time,
dnl so we get into trouble with duplicate symbols. Check
dnl whether this is so.
dnl
dnl The check is a little more complicated than the usual one,
dnl since we can't use AC_TRY_COMPILE (we have to compile and link
dnl more than one file).
dnl
dnl Usage: DEAL_II_CHECK_ANON_NAMESPACE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ANON_NAMESPACE_BUG, dnl
[
  AC_MSG_CHECKING(for anonymous namespace and name mangling bug)

  dnl Create the testfile
  echo "namespace { int SYMBOL() {return 1;}; }" >  conftest.cc
  echo "static int f() { return SYMBOL(); }"     >> conftest.cc

  dnl Then compile it twice...
  $CXX -c conftest.cc -o conftest.1.$ac_objext
  $CXX -c conftest.cc -o conftest.2.$ac_objext

  dnl Create a file with main() and also compile it
  echo "int main () {}" > conftest.cc
  $CXX -c conftest.cc -o conftest.3.$ac_objext

  dnl Then try to link everything
  if (($CXX conftest.1.$ac_objext conftest.2.$ac_objext \
            conftest.3.$ac_objext -o conftest  2>&1 ) > /dev/null );\
  then
      AC_MSG_RESULT(no)
  else
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_ANON_NAMESPACE_BUG, 1,
                     [Defined if the compiler needs to see the static
                      keyword even for functions in anonymous namespaces,
                      to avoid duplicate symbol errors when linking.
                      For the details, look at aclocal.m4 in the
                      top-level directory.])
  fi
  rm -f conftest.1.$ac_objext conftest.2.$ac_objext conftest
])



dnl -------------------------------------------------------------
dnl A second test in this direction: if the name of a function is
dnl not mangled differently for each compiler invokation, then
dnl it should at least result in a weak symbol. Test this.
dnl
dnl Note that this is not a problem in itself if the name is
dnl mangled differently each time a file is compiled.
dnl
dnl Usage: DEAL_II_CHECK_ANON_NAMESPACE_BUG2
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ANON_NAMESPACE_BUG2, dnl
[
  AC_MSG_CHECKING(for anonymous namespace and weak linkage bug)

  dnl Create the testfile
  echo "namespace { int SYMBOL() {return 1;}; }" >  conftest.cc
  echo "static int f() { return SYMBOL(); }"     >> conftest.cc

  dnl Compile it
  $CXX -c conftest.cc -o conftest.$ac_objext

  dnl Then look for lines in the output of "nm" that have the name of
  dnl SYMBOL in them. Then make sure that we don't find a line with
  dnl " T " in it, i.e. a text symbol with strong linkage
  check="`nm conftest.$ac_objext | grep SYMBOL | grep ' T '`"

  dnl Then try to link everything
  if test "x$check" = "x" ;
  then
      AC_MSG_RESULT(no)
  else
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG, 1,
                     [Another test if the compiler needs to see the static
                      keyword even for functions in anonymous namespaces,
                      to avoid duplicate symbol errors when linking.
                      For the details, look at aclocal.m4 in the
                      top-level directory.])
  fi
  rm -f conftest.$ac_objext
])


dnl -------------------------------------------------------------
dnl We have so many templates in deal.II that sometimes we need
dnl to make it clear with which types a template parameter can
dnl be instantiated. There is a neat trick to do this: SFINAE
dnl (substitution failure is not an error). The idea is this: the
dnl C++ standard prescribes that a template function is only
dnl considered in a call, if all parts of its signature can be
dnl instantiated with the template parameter replaced by the
dnl respective types/values in this particular call. Example:
dnl   template <typename T>
dnl   typename T::type  foo(T) {...};
dnl   ...
dnl   foo(1);
dnl
dnl The compiler should detect that in this call, the template
dnl parameter T must be identified with the type "int". However,
dnl the return type T::type does not exist. The trick now is
dnl that this is not considered an error: this template is simply
dnl not considered, the compiler keeps on looking for another
dnl possible function foo.
dnl
dnl That allows for a neat trick to rule out a template for certain
dnl template arguments without changing the function signature at
dnl all: Make the return type un-instantiable if the template
dnl type presently considered for instantiation does not qualify.
dnl An example of this is shown below.
dnl
dnl Unfortunately, older compilers do not support this trick: they
dnl issue an error if the return type cannot be installed. In this
dnl case, we back out and just drop the constraint and allow for
dnl all template types. In that case, you'll simply get compiler
dnl error when it tries to compile this template (in case it's in
dnl the .h file), or linker errors for missing functions in case
dnl it's in the .cc file and explicitly instantiated.
dnl
dnl Usage: DEAL_II_CHECK_SFINAE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_SFINAE_BUG, dnl
[
  AC_MSG_CHECKING(for SFINAE bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <typename> struct int_or_double { static const bool value = false;};
      template <> struct int_or_double<int>    { static const bool value = true; };
      template <> struct int_or_double<double> { static const bool value = true; };

      template <bool, typename> struct constraint_and_return_value {};

      template <typename T> struct constraint_and_return_value<true,T>
      {
          typedef T type;
      };

      // deduction for T=char should file, since return type cannot be
      // instantiated...
      template <typename T>
      typename constraint_and_return_value<int_or_double<T>::value,void>::type
      f (T);

      // ...however, this is not an error. rather, the compiler should
      // instead choose the following function:
      void f(int);
    ],
    [
      f('c');
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. disabling template constraints)
      AC_DEFINE(DEAL_II_SFINAE_BUG, 1,
                     [Defined if the compiler does not support the
                      substitution-failure-is-not-an-error paradigm.
                      For the details, look at aclocal.m4 in the
                      top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl Some versions of gcc get this example wrong:
dnl ---------------------------------
dnl struct X
dnl {
dnl     template <typename T> void operator << (T);
dnl };
dnl
dnl void f()
dnl {
dnl   X x;
dnl   x.operator << <double> (1);
dnl }
dnl ---------------------------------
dnl They want to see a "template" for disambiguation in
dnl    x.template operator << <double> (1);
dnl which shouldn't be necessary since the left hand side of the
dnl dot operator is not template dependent. Surprisingly, this is
dnl only the case for operators, not if operator<< were a regular
dnl function. Annoyingly, other compilers barf on seeing the
dnl disambiguating "template" keyword.
dnl
dnl Usage: DEAL_II_CHECK_TEMPL_OP_DISAMBIGUATION_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPL_OP_DISAMBIGUATION_BUG, dnl
[
  AC_MSG_CHECKING(for template operator disambiguation bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      struct X
      {
          template <typename T> void operator << (T);
      };

      void f()
      {
        X x;
        x.operator << <double> (1);
      }
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG, 1,
                     [Defined if the compiler requires the use of the
                      template keyword for disambiguation keyword in
                      certain contexts in which it is not supposed to
                      do so. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl Some older versions of gcc compile this, despite the 'explicit'
dnl keyword:
dnl ---------------------------------
dnl struct X {
dnl     template <typename T>
dnl     explicit X(T);
dnl };
dnl
dnl void f(X);
dnl
dnl int main () { f(1); }
dnl ---------------------------------
dnl
dnl Check for this misfeature.
dnl
dnl Usage: DEAL_II_CHECK_EXPLICIT_CONSTRUCTOR_BUG
dnl
dnl --------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_EXPLICIT_CONSTRUCTOR_BUG, dnl
[
  AC_MSG_CHECKING(for explicit template constructor bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      struct X {
        template <typename T>
        explicit X(T);
      };

      void f(X);
    ],
    [
      f(1);
    ],
    [
      AC_MSG_RESULT(yes. disabling some functions)
      AC_DEFINE(DEAL_II_EXPLICIT_CONSTRUCTOR_BUG, 1,
                     [Defined if the compiler does not honor the explicit
                      keyword on template constructors.])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Some older versions of gcc deduce pointers to const functions in
dnl template contexts to pointer-to-function of const objects.
dnl This is not correct
dnl
dnl Check for this misfeature.
dnl
dnl Usage: DEAL_II_CHECK_CONST_MEMBER_DEDUCTION_BUG
dnl
dnl --------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CONST_MEMBER_DEDUCTION_BUG, dnl
[
  AC_MSG_CHECKING(for const member deduction bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <typename T> struct identity { typedef T type; };

      template <typename C> void new_thread (void (C::*fun_ptr)(),
                                       typename identity<C>::type &c);
      template <typename C> void new_thread (void (C::*fun_ptr)() const,
                                       const typename identity<C>::type &c);
      struct X { void f() const; };
    ],
    [
      X x;
      new_thread (&X::f, x);
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_CONST_MEMBER_DEDUCTION_BUG, 1,
                     [Defined if the compiler has a bug in deducing
                      the type of pointers to const member functions.])
    ])
])



dnl -------------------------------------------------------------
dnl Check for GCC bug 36052, see
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=36052
dnl
dnl Usage: DEAL_II_CHECK_TYPE_QUALIFIER_BUG
dnl
dnl --------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TYPE_QUALIFIER_BUG, dnl
[
  case "$GXX_VERSION" in
    gcc*)
      AC_MSG_CHECKING(for warning bug with type qualifiers)
      AC_LANG(C++)
      CXXFLAGS="$CXXFLAGSG -Werror"
      AC_TRY_COMPILE(
        [
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
        ],
        [
        ],
        [
          AC_MSG_RESULT(no)
        ],
        [
          AC_MSG_RESULT(yes)
          CXXFLAGSG="$CXXFLAGSG -Wno-ignored-qualifiers"
        ])
    ;;
  esac
])




dnl -------------------------------------------------------------
dnl In the gcc libstdc++ headers for std::complex, there is
dnl no defined default copy constructor, but a templated copy
dnl constructor. So when using using a normal assignment
dnl between identical types, the compiler synthesizes the
dnl default operator, rather than using the template.
dnl The code will still be ok, though.
dnl
dnl With -Wsynth in gcc we then get a warning. So if we find that
dnl this is still the case, disable -Wsynth, i.e. remove it from
dnl the list of warning flags.
dnl
dnl This is gcc bug 18644:
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=18644
dnl
dnl Usage: DEAL_II_CHECK_WSYNTH_AND_STD_COMPLEX
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_WSYNTH_AND_STD_COMPLEX, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for problem with -Wsynth and std::complex)
    AC_LANG(C++)
    CXXFLAGS="-Wsynth -Werror"
    AC_TRY_COMPILE(
      [
#       include <complex>
      ],
      [
        std::complex<double> x;
        x = std::complex<double>(1,0);
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
        CXXFLAGSG="`echo $CXXFLAGSG | perl -pi -e 's/-Wsynth\s*//g;'`"
      ])
  fi
])



dnl -------------------------------------------------------------
dnl Older gcc versions appear to frown upon the way we write the
dnl IsBlockMatrix<MatrixType> template. If that's the case,
dnl remove the -Wctor-dtor-privacy flag.
dnl
dnl Usage: DEAL_II_CHECK_CTOR_DTOR_PRIVACY
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CTOR_DTOR_PRIVACY, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for problem with -Wctor-dtor-privacy)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSG -Werror"
    AC_TRY_COMPILE(
      [
        template <typename T>
        struct IsInt
        {
          private:
            struct yes_type { char c[1]; };
            struct no_type  { char c[2]; };

            static yes_type check_for_int (const int *);

            static no_type check_for_int (...);

          public:
            static const bool value = (sizeof(check_for_int((T*)0))
                                       ==
                                       sizeof(yes_type));
        };

        const bool x = IsInt<double>::value;
      ],
      [
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
        CXXFLAGSG="$CXXFLAGSG -Wno-ctor-dtor-privacy"
      ])
  fi
])



dnl -------------------------------------------------------------
dnl On Mac OS X, gcc appears to have a bug that prevents us from
dnl compiling a bit of code that involves boost::bind. Check for
dnl that.
dnl
dnl Usage: DEAL_II_CHECK_BOOST_BIND_COMPILER_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BOOST_BIND_COMPILER_BUG, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for boost::bind compiler internal error)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSO $BOOST_INCLUDE_DIR"
    AC_TRY_COMPILE(
      [
#include <complex>
#include <iostream>
#include <boost/bind.hpp>

template<typename number>
void bug_function (number test)
{
  std::cout << test << std::endl;
}
      ],
      [
  std::complex<float> float_val (1., 2.);
  boost::bind(&bug_function<std::complex<float> >,
              float_val)();
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(DEAL_II_BOOST_BIND_COMPILER_BUG, 1,
                  [Defined if the compiler gets an internal error compiling
                   some code that involves boost::bind])
      ])
  fi
])



dnl -------------------------------------------------------------
dnl Some versions of icc on some platforms issue a lot of warnings
dnl about the unreliability of floating point comparisons. Check
dnl whether we can switch that off by checking whether the compiler
dnl allows -wdXXXX for this warning:
dnl #1572: `floating-point equality and inequality comparisons are
dnl        unreliable'
dnl        (while true, this warning also triggers on comparisons
dnl        with zero, or comparing two variables for which one is
dnl        greater; there is about no way to write numeric code
dnl        without triggering this warning many times over)
dnl
dnl Usage: DEAL_II_ICC_WD_1572
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_ICC_WD_1572, dnl
[
  AC_MSG_CHECKING(whether -wd1572 is allowed)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -wd1572"
  AC_TRY_COMPILE( [], [],
      [
        AC_MSG_RESULT(yes)
        CXXFLAGSG="$CXXFLAGSG -wd1572"
      ],
      [
        AC_MSG_RESULT(no)
      ])
])



dnl -------------------------------------------------------------
dnl Check whether socket functions (i.e. functions using the
dnl network) affect the numerical accuracy. To spice up things,
dnl on some versions of cygwin (<=1.5.13), calling the windows
dnl socket function (e.g. the functions defined in unistd.h)
dnl will set the internal accuracy of the FPU to 64 bits (double)
dnl instead of 80 bits. Consequently, all long double operations
dnl are performed with double accuracy only. The compiler will
dnl not notice this fact and still report long double limits in
dnl the linits<>-class. We have to check for this, otherwise
dnl some functions like the QGauss constructor will go into an
dnl infinite loop.
dnl
dnl Usage: DEAL_II_CHECK_BROKEN_SOCKETS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BROKEN_SOCKETS, dnl
[
  AC_MSG_CHECKING(for bad socket functions/FPU interaction)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_RUN(
    [
#include <unistd.h>
#include <limits>
#include <iostream>

int main()
{
   char buf[100];
   gethostname(buf,99);

   volatile long double x=1.0;
   x += std::numeric_limits<long double>::epsilon();

   if (x == 1.0)
     return 1;  // this shouldn't happen...
   else
     return 0;
}
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT([yes. disabling socket functions])
      AC_DEFINE(DEAL_II_BROKEN_SOCKETS, 1,
                [Define if the use of socket functionality leads to strange
                 results with floating point computations on cygwin
                 systems.])
    ])
])

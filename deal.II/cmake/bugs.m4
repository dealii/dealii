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
dnl On some systems, notably AIX and SUN Solaris, using threads
dnl leads to warnings since the POSIX_MUTEX_INITIALIZER preprocessor
dnl variable used to initialize POSIX mutex objects does not contain
dnl initializers for all elements of the mutex. This is not wrong,
dnl but leads to the message "warning: aggregate has a partly
dnl bracketed initializer", which is annoying since it shows up
dnl _very_ often in our files, although this is something that
dnl happens inside gcc systems headers. So avoid the warning if
dnl necessary
dnl
dnl Usage:
dnl   DEAL_II_CHECK_CHECK_PARTLY_BRACKETED_INITIALIZER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_PARTLY_BRACKETED_INITIALIZER, dnl
[
  if test "$enablethreads" = yes ; then
    case "$GXX_VERSION" in
      gcc*)
        AC_MSG_CHECKING(for only partly bracketed mutex initializer)
        AC_LANG(C++)
        CXXFLAGS="$CXXFLAGSG -Werror"
        AC_TRY_COMPILE(
        [
#       include <vector>
        ],
        [;],
        [
          AC_MSG_RESULT(no)
        ],
        [
          AC_MSG_RESULT(yes)
          CXXFLAGSG="$CXXFLAGSG -Wno-missing-braces"
          CXXFLAGSO="$CXXFLAGSO -Wno-missing-braces"
        ])
        ;;
      *)
        ;;
    esac
  fi
])



dnl -------------------------------------------------------------
dnl Check whether some backward compatibility features are disabled
dnl Usage:
dnl   DEAL_II_CHECK_COMPAT_BLOCKER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_COMPAT_BLOCKER, dnl
[
  AC_ARG_ENABLE(compat-blocker,
              [AS_HELP_STRING([--enable-compat-blocker=mapping],
              [Block functions that implicitely assume a Q1 mapping])],
      enable_compat_blocker="$enableval",
      enable_compat_blocker="")

  dnl Replace the comma-separated list by a space-separated one
  disable_compat=`echo $enable_compat_blocker | perl -pi -e 's/,/ /g;'`

  dnl Check that each entry is an allowed one
  for i in $disable_compat ; do
    case $i in
      mapping)
        AC_MSG_RESULT(Disabling backward compatibility feature: "$i")
        ;;
      *)
        AC_MSG_ERROR(Backward compatibility feature "$i" unknown)
        ;;
    esac
  done

  dnl Now for each known feature, either disable it or enable it.
  dnl Default is to enable. In order to have these flags in config.h,
  dnl it is necessary to AC_DEFINE them actually by name, rather than
  dnl by some loop variable, since otherwise autoheader can't generate
  dnl an entry for config.h for this variable
  for i in mapping ; do
    uppercase=`echo $i | perl -pi -e 'tr/a-z/A-Z/;'`
    if test -n "`echo $disable_compat | grep $i`"  ; then
      compat_value=false
    else
      compat_value=true
    fi

    case $i in
      mapping)
        AC_DEFINE_UNQUOTED(DEAL_II_COMPAT_MAPPING,$compat_value,
                           [Backward compatibility support for functions
                              and classes that do not take an explicit
                              mapping variable, but rather use a default
                              Q1 mapping instead])
        ;;

      *)
        AC_MSG_ERROR(Backward compatibility feature "$i" unknown)
        ;;
    esac
  done
])



dnl -------------------------------------------------------------
dnl On SunOS 4.x, the `getrusage' function exists, but is not declared
dnl in the respective header file `resource.h', as one would think when
dnl reading the man pages. Then we have to declare this function
dnl ourselves in those files that use this function. The question whether
dnl we have to do so is controlled by a preprocessor variable.
dnl
dnl If the function is not properly declared, then set the preprocessor
dnl variable NO_HAVE_GETRUSAGE
dnl
dnl Usage: DEAL_II_CHECK_GETRUSAGE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_GETRUSAGE, dnl
[
  AC_MSG_CHECKING(whether getrusage is properly declared)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <sys/resource.h>
    ],
    [
      rusage *ru;
      getrusage(RUSAGE_SELF,ru);
    ],
    [
        AC_MSG_RESULT(yes)
    ],
    [
        AC_MSG_RESULT(no)
        AC_DEFINE(NO_HAVE_GETRUSAGE, 1,
                  [On SunOS 4.x, the getrusage() function exists, but
                   is not declared in the respective header file
                   <resource.h>, as one would think when reading the
                   man pages. Then we have to declare this function
                   ourselves in those files that use this function.
                   The question whether we have to do so is controlled
                   by the preprocessor variable.])
    ])
]
)



dnl -------------------------------------------------------------
dnl On some systems (well, DEC Alphas are the only ones we know of),
dnl gcc2.95 throws the hands in the air if it sees one of the AssertThrow
dnl calls, and dies with an internal compiler error. If this is the case,
dnl we disable AssertThrow and simply replace it with an `abort' if the
dnl condition is not satisfied.
dnl
dnl Usage: DEAL_II_CHECK_ASSERT_THROW("description of options set",
dnl                                   "compiler options set",
dnl                                   action if compiler crashes)
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ASSERT_THROW, dnl
[
  AC_MSG_CHECKING(whether AssertThrow works with $1 flags)
  AC_LANG(C++)
  CXXFLAGS="$2"
  AC_TRY_COMPILE(
    [
#include <exception>
#include <iostream>
#include <cstdlib>

#ifndef __GNUC__
#  define __PRETTY_FUNCTION__ "(unknown)"
#endif
class ExceptionBase : public std::exception {
  public:
    ExceptionBase ();
    ExceptionBase (const char* f, const int l, const char *func,
                   const char* c, const char *e);
    virtual ~ExceptionBase () throw();
    void SetFields (const char *f, const int   l, const char *func,
                    const char *c, const char *e);
    void PrintExcData (std::ostream &out) const;
    virtual void PrintInfo (std::ostream &out) const;
    virtual const char * what () const throw ();
  protected:
    const char  *file;
    unsigned int line;
    const char  *function, *cond, *exc;
};

template <class exc>
void __IssueError_Assert (const char *file,
                          int         line,
                          const char *function,
                          const char *cond,
                          const char *exc_name,
                          exc         e){
  e.SetFields (file, line, function, cond, exc_name);
  std::cerr << "--------------------------------------------------------"
            << std::endl;
  e.PrintExcData (std::cerr);
  e.PrintInfo (std::cerr);
  std::cerr << "--------------------------------------------------------"
            << std::endl;
  std::abort ();
}

template <class exc>
void __IssueError_Throw (const char *file,
                         int         line,
                         const char *function,
                         const char *cond,
                         const char *exc_name,
                         exc         e) {
                                   // Fill the fields of the exception object
  e.SetFields (file, line, function, cond, exc_name);
  throw e;
}

#define AssertThrow(cond, exc)                                    \
  {                                                               \
    if (!(cond))                                                  \
      __IssueError_Throw (__FILE__,                               \
                          __LINE__,                               \
                          __PRETTY_FUNCTION__, #cond, #exc, exc); \
  }

#define DeclException0(Exception0)  \
class Exception0 :  public ExceptionBase {}

namespace StandardExceptions
{
  DeclException0 (ExcInternalError);
}
using namespace StandardExceptions;
    ],
    [
        AssertThrow (false, ExcInternalError());
    ],
    [
        AC_MSG_RESULT(yes)
    ],
    [
        AC_MSG_RESULT(no)
        $3
    ])
])



dnl -------------------------------------------------------------
dnl Sun's Forte compiler (at least up to the Version 7 Early Access
dnl release) has a problem with the following code, when compiling
dnl with debug output:
dnl
dnl /* ---------------------------------------------------------- */
dnl /* Internal compiler error in abi2_mangler::entity_expression */
dnl /* when compiled with -g.                                     */
dnl template < int dim > struct T {
dnl     typedef T<dim-1> SubT;
dnl     T (SubT);
dnl };
dnl
dnl template <int dim> T<dim>::T (SubT) {};
dnl
dnl template class T<3> ;
dnl /* ---------------------------------------------------------- */
dnl
dnl The compiler gets an internal compiler error, so we work around
dnl this problem by a really evil hack in the sources.
dnl
dnl Usage: DEAL_II_CHECK_LOCAL_TYPEDEF_COMP
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_LOCAL_TYPEDEF_COMP, dnl
[
  AC_MSG_CHECKING(for local computed template typedef bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -g"
  AC_TRY_COMPILE(
    [
        template < int dim > struct T {
            typedef T<dim-1> SubT;
            T (SubT);
        };

        template <int dim> T<dim>::T (SubT) {}

        template class T<3>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_LOCAL_TYPEDEF_COMP_WORKAROUND, 1,
                [Define if we have to work around a bug in Sun's Forte compiler.
                 See the aclocal.m4 file in the top-level directory for a
                 description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl Sun's Forte compiler (at least up to the Version 7 Early Access
dnl release) have a problem with the following code, when compiling
dnl with debug output:
dnl
dnl /* ----------------------------------------------- */
dnl /* Problem 14: Access control. Friendship is not   */
dnl /* granted although explicitly declared.           */
dnl template <int N, int M> class T      {    int bar ();  };
dnl
dnl template <int M>        class T<1,M> {
dnl   private:
dnl     static int i;
dnl     template <int N1, int N2> friend class T;
dnl };
dnl
dnl template <int N,int M> int T<N,M>::bar () {
dnl   return T<N-1,M>::i;
dnl };
dnl
dnl template class T<2,1> ;
dnl /* ---------------------------------------------------------- */
dnl
dnl The compiler does not allow access to T<1,1>::i, although the
dnl accessing class is explicitly marked friend.
dnl
dnl Usage: DEAL_II_CHECK_TEMPLATE_SPEC_ACCESS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPLATE_SPEC_ACCESS, dnl
[
  AC_MSG_CHECKING(for partially specialized template access control bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
        template <int N, int M> struct T      {    int bar ();  };

        template <int M>        struct T<1,M> {
            T ();
          private:
            static int i;
            template <int N1, int N2> friend class T;
        };

        template <int N,int M> int T<N,M>::bar () {
          return T<N-1,M>::i;
        }

        template class T<2,1>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPLATE_SPEC_ACCESS_WORKAROUND, 1,
                [Define if we have to work around a bug in Sun's Forte compiler.
                 See the aclocal.m4 file in the top-level directory for a
                 description of this bug.])
    ])
])


dnl -------------------------------------------------------------
dnl Versions of GCC before 3.0 had a problem with the following
dnl code:
dnl
dnl /* ----------------------------------------------- */
dnl namespace NS {
dnl   template <typename T>  class C  {
dnl       template <typename N> friend class C;
dnl   };
dnl };
dnl /* ----------------------------------------------- */
dnl
dnl This is fixed with gcc at least in snapshots before version 3.1,
dnl but the following bug remains:
dnl
dnl /* ----------------------------------------------- */
dnl namespace NS {  template <typename number> class C;  };
dnl
dnl template <typename T> class X {
dnl   template <typename N> friend class NS::C;
dnl };
dnl
dnl template class X<int>;
dnl /* ----------------------------------------------- */
dnl
dnl The compiler gets an internal error for these cases. Since we need this
dnl construct at various places, we check for it and if the compiler
dnl dies, we use a workaround that is non-ISO C++ but works for these
dnl compilers.
dnl
dnl Usage: DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for 1st template friend in namespace bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
        namespace NS {
          template <typename T>  class C  {
              C(const C<T>&);
              template <typename N> friend class C;
          };
        }

        namespace NS2 {  template <typename number> class C;  }

        template <typename T> class X {
            template <typename N> friend class NS2::C;
            template <typename N> friend class NS::C;
        };

        template class X<int>;

        namespace NS {
          template<typename T>
          inline C<T>::C(const C<T>&)
          {}
        }
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE_UNQUOTED(DEAL_II_NAMESP_TEMPL_FRIEND_BUG, 1,
                         [Define if we have to work around a bug in gcc with
                          marking all instances of a template class as friends
                          to this class if the class is inside a namespace.
                          See the aclocal.m4 file in the top-level directory
                          for a description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl Another bug in gcc with template and namespaces (fixed since 3.2,
dnl but present in 3.0):
dnl
dnl /* ----------------------------------------------- */
dnl namespace NS {
dnl   template <typename> struct Foo;
dnl }
dnl
dnl class Bar {
dnl     template <typename Y> friend struct NS::Foo;
dnl };
dnl
dnl namespace NS {
dnl   template <typename> struct Foo { Foo (); };
dnl }
dnl
dnl template struct NS::Foo<int>;
dnl /* ----------------------------------------------- */
dnl
dnl gcc2.95 provides a really unhelpful error message, 3.0 gets an
dnl internal compiler error.
dnl
dnl Usage: DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG2
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG2, dnl
[
  AC_MSG_CHECKING(for 2nd template friend in namespace bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -Werror"
  AC_TRY_COMPILE(
    [
        namespace NS {
          template <typename> struct Foo;
        }

        class Bar {
            template <typename Y> friend struct NS::Foo;
        };

        namespace NS {
          template <typename> struct Foo { Foo (); };
        }

        template struct NS::Foo<int>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE_UNQUOTED(DEAL_II_NAMESP_TEMPL_FRIEND_BUG2, 1,
                         [Define if we have to work around another bug in gcc with
                          marking all instances of a template class as friends
                          to this class if the class is inside a namespace.
                          See the aclocal.m4 file in the top-level directory
                          for a description of this bug.])
    ])
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
dnl DEC/Compaq's cxx compiler does not want us to implement
dnl virtual functions that were declared abstract before. We do
dnl this with the destructor of the Function class, since we want
dnl to avoid people making objects of that class, but all functions
dnl have default implementations, so the class is not abstract
dnl without that. Since every derived class has a destructor, it
dnl is sufficient to mark the destructor =0.
dnl
dnl Unfortunately, cxx refuses to grok that. It sees the respective
dnl function in the .cc file, but does not instantiate it, leading
dnl to linker errors. Thus, check this misfeature, and if present
dnl simply do not mark the function abstract for this particular
dnl compiler.
dnl
dnl Usage: DEAL_II_CHECK_IMPLEMENTED_PURE_FUNCTION_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_IMPLEMENTED_PURE_FUNCTION_BUG, dnl
[
  AC_MSG_CHECKING(for bug with implementing pure functions)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
        template <int dim>
        struct Function
        {
          public:
            virtual ~Function () = 0;
        };

        template <int dim>
        Function<dim>::~Function ()
        {}

        template class Function<1>;
        template Function<1>::~Function();
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_IMPLEMENTED_PURE_FUNCTION_BUG, 1,
                     [Defined if the compiler refuses to compile the definition
                      of a function that was previously declared abstract.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc 2.95 dies on this construct:
dnl -----------------------------
dnl template <int dim> struct TT { typedef int type; };
dnl
dnl template <template <int> class T> struct X {
dnl     typedef typename T<1>::type type;
dnl     void foo (type t);
dnl };
dnl
dnl template <template <int> class T>
dnl void X<T>::foo (type t) {};
dnl
dnl template struct X<TT>;
dnl -----------------------------
dnl
dnl We work around this problem, if we encounter it.
dnl
dnl Usage: DEAL_II_CHECK_TEMPLATE_TEMPLATE_TYPEDEF_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPLATE_TEMPLATE_TYPEDEF_BUG, dnl
[
  AC_MSG_CHECKING(for template template typedef bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
        template <int dim> struct TT { typedef int type; };

        template <template <int> class T> struct X {
            typedef typename T<1>::type type;
            void foo (type t);
        };

        template <template <int> class T>
        void X<T>::foo (type) {}

        template struct X<TT>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPLATE_TEMPLATE_TYPEDEF_BUG, 1,
                     [Defined if the compiler refuses to allow a typedef
                      to a template template class template parameter. For
                      the exact failure mode, look at aclocal.m4 in the
                      top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc 2.95 as well as some other compilers do not correctly implement
dnl the resolution of defect report #45 to the C++ standard (see
dnl http://anubis.dkuug.dk/jtc1/sc22/wg21/docs/cwg_active.html#45).
dnl try to detect this, and set a flag correspondingly. in short,
dnl the DR says that this is allowed, since member classes are
dnl implicitly also friends:
dnl -----------------------------
dnl struct X {
dnl     X ();
dnl   private:
dnl     static int f();
dnl
dnl     struct Y {
dnl         int g() { return f(); };
dnl     };
dnl };
dnl -----------------------------
dnl
dnl We work around this problem, if we encounter it.
dnl
dnl Usage: DEAL_II_CHECK_NESTED_CLASS_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NESTED_CLASS_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for nested classes are implicit friends bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
        struct X {
            X ();
          private:
            static int f();

            struct Y {
                int g() { return f(); };
            };
        };
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_NESTED_CLASS_FRIEND_BUG, 1,
                     [Defined if the compiler does not properly implement
                      the resolution of defect report #45 to the C++
                      standard, which makes nested types implicit friends
                      of the enclosing class.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc up to 3.3 won't accept the following code:
dnl -----------------------------
dnl template <typename> class X {
dnl     template <typename> class Y {};
dnl
dnl     template <typename T>
dnl     template <typename>
dnl     friend class X<T>::Y;
dnl };
dnl
dnl X<int> x;
dnl -----------------------------
dnl
dnl They don't accept the X<T>::Y here, probably because the class is
dnl not complete at this point. gcc3.4 gets it right, though. One
dnl can work around by simply saying
dnl   template <typename> friend class Y;
dnl but then icc doesn't understand this :-( So everyone's got a bug
dnl here. Also, note that the standard says that Y is an implicit friend
dnl of X, but again, many compiler don't implement this correctly, which
dnl is why we have to do something like the above in the first place...
dnl
dnl Usage: DEAL_II_CHECK_NESTED_CLASS_TEMPL_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NESTED_CLASS_TEMPL_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for nested template class friends bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <typename> class X {
          template <typename> class Y {};

          template <typename T>
          template <typename>
          friend class X<T>::Y;
      };

      X<int> x;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_NESTED_CLASS_TEMPL_FRIEND_BUG, 1,
                     [Defined if the compiler does not understand friend
                      declarations for nested member classes when giving
                      a full class specification.])
    ])
])



dnl -------------------------------------------------------------
dnl Many compilers get this wrong (see Section 14.7.3.1, number (4)):
dnl ---------------------------------
dnl   template <int dim> struct T {
dnl     static const int i;
dnl   };
dnl
dnl   template <> const int T<1>::i;
dnl   template <> const int T<1>::i = 1;
dnl ---------------------------------
dnl First, by Section 14.7.3.14 of the standard, the first template<>
dnl line must necessarily be the _declaration_ of a specialization,
dnl and the second is then its definition. There is therefore no
dnl reason to report a doubly defined variable (Intel ICC 6.0), or
dnl to choke on these lines at all (Sun Forte)
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_VAR_SPECIALIZATION_SPEC_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_VAR_SPECIALIZATION_BUG, dnl
[
  AC_MSG_CHECKING(for member variable specialization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
        template <int dim> struct T {
            static const int i;
        };

        template <> const int T<1>::i;
        template <> const int T<1>::i = 1;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG, 1,
                     [Defined if the compiler refuses to allow the
                      explicit specialization of static member
                      variables. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
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
dnl gcc 2.95 doesn't like it if we have a member template function
dnl and define it as a template while specializing the outer class
dnl template. This is a nasty bug that is hard to work around...
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_TEMPLATE_SPECIALIZATION_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_TEMPLATE_SPECIALIZATION_BUG, dnl
[
  AC_MSG_CHECKING(for template member function specialization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <int dim> struct X
      {
        template <typename T> void f(T);
      };

      template <>
      template <typename T>
      void X<1>::f (T)
      {}

      template void X<1>::f(int);
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_MEMBER_TEMPLATE_SPECIALIZATION_BUG, 1,
                     [Defined if the compiler refuses to specialize
                      an outer class template while keeping a member
                      as a template. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc3.1 (and maybe later compilers) has a bug with long double
dnl and optimization (see code below), when compiling on Sparc
dnl machines. Since it affects only one platform and one compiler,
dnl we take the liberty to disable the function in which the problem
dnl occurs (Polynomial::shift in base/source/polynomial.cc), since
dnl this is a function that is rarely used anyway.
dnl
dnl For more information: the bug just described is reported to
dnl the gcc project under number 7335.
dnl
dnl Usage: DEAL_II_CHECK_LONG_DOUBLE_LOOP_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_LONG_DOUBLE_LOOP_BUG, dnl
[
  AC_MSG_CHECKING(for long double optimization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSO"
  AC_TRY_COMPILE(
    [
        double* copy(long double* first, long double* last, double* result)
        {
          int n;
          for (n = last - first; n > 0; --n) {
            *result = *first;
            ++first;
            ++result;
          }
          return result;
        }

        void f()
        {
          long double *p1=0, *p2=0;
          double *p3=0;
          copy (p1, p2, p3);
          p3 = copy (p1, p2, p3);
        };
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. disabling respective functions)
      AC_DEFINE(DEAL_II_LONG_DOUBLE_LOOP_BUG, 1,
                     [Defined if the compiler gets an internal compiler
                      upon some code involving long doubles, and with
                      optimization. For the details, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc2.95 (but not later compilers) has a bug with taking the
dnl address of a function with template template parameters (or
dnl with calling this function by specifying explicitly the template
dnl arguments). This requires some working around that in turn does
dnl not work with later compilers.
dnl
dnl Usage: DEAL_II_CHECK_FUNPTR_TEMPLATE_TEMPLATE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_FUNPTR_TEMPLATE_TEMPLATE_BUG, dnl
[
  AC_MSG_CHECKING(for address of template template function bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <int> struct X {};

      template <int dim, template <int> class T>
      void f(T<dim>);

      template <int dim, template <int> class T>
      void* g()
      {
        void (*p) (T<dim>) = &f<dim,T>;
        return (void*)p;
      }

      template void* g<2,X> ();
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_FUNPTR_TEMPLATE_TEMPLATE_BUG, 1,
                     [Defined if the compiler needs a workaround for
                      certain problems with taking the address of
                      template template functions. For the details, look at
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
dnl A third test: MIPSpro gives a bogus warning about unused code
dnl on code that is clearly used, when inline member functions
dnl are in classes in anonymous namespaces. Check for this to allow
dnl us to work around this problem
dnl
dnl Usage: DEAL_II_CHECK_ANON_NAMESPACE_BUG3
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ANON_NAMESPACE_BUG3, dnl
[
  AC_MSG_CHECKING(for bogus warning with anonymous namespaces)

  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  case "$GXX_VERSION" in
    gcc*)
        CXXFLAGS="$CXXFLAGSG -Werror"
        ;;

    MIPSpro*)
        CXXFLAGS="$CXXFLAGSG -diag_error 1174"
        ;;
  esac

  AC_TRY_COMPILE(
   [
     template<int dim> class Point
     {
       public:
         Point() {};
     };

     class GridTools {
       public:
         template <int dim> static void shift (const Point<dim> &s);

         template <typename Predicate>
         static void transform (const Predicate    &predicate);
     };

     namespace {
       template <int dim> class ShiftPoint {
         public:
           ShiftPoint (const Point<dim>   &) {};
           void g() const;
       };

       template<int dim> void ShiftPoint<dim>::g() const {}
     }

     template <typename Predicate>
     void GridTools::transform (const Predicate &predicate)
     { predicate.g(); }

     template <int dim>
     void GridTools::shift (const Point<dim>   &s)
     { transform (ShiftPoint<dim>(s)); }

     template void GridTools::shift (const Point<1>   &);
   ],
   [
        ;
   ],
   [
        AC_MSG_RESULT(no)
   ],
   [
        AC_MSG_RESULT(yes)
        AC_DEFINE_UNQUOTED(DEAL_II_ANON_NAMESPACE_BOGUS_WARNING, 1,
                     [Flag indicating whether there is a bug in the
                      compiler that leads to bogus warnings for
                      inline class members in anonymous namespaces])
   ])
])



dnl -------------------------------------------------------------
dnl A test to identify Apples buggy gcc3.3 version. If the
dnl explicit instantiations are placed at the end of a source
dnl file, sometimes only weak symbols are generated, which
dnl lead to linker problems.
dnl For more information, see
dnl http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24331
dnl As this bug can be easily avoided by applying the
dnl "november2004gccupdate" Patch which can be found
dnl on http://www.apple.com/developer, this check will
dnl only produce a warning at the end of the configuration
dnl process.
dnl
dnl Usage: DEAL_II_CHECK_WEAK_LINKAGE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_WEAK_LINKAGE_BUG, dnl
[
  DARWIN_GCC_WEAK_LINKAGE_BUG="no"
  case "$target" in
    *-apple-darwin* )

    AC_MSG_CHECKING(for weak linkage bug (Apple gcc3.3) )

    dnl Create 1st testfile
    echo "template <int> void SYMBOL() {}" >  conftest.cc
    echo "int g() { SYMBOL<1> (); }" >>  conftest.cc
    echo "template void SYMBOL<1> ();"     >> conftest.cc

    dnl Compile it
    $CXX -c conftest.cc -o conftest.$ac_objext

    dnl and write symbols to file 1
    nm conftest.$ac_objext | grep SYMBOL > symb1

    dnl Create 2nd testfile
    echo "template <int> void SYMBOL() {}" >  conftest.cc
    echo "template void SYMBOL<1> ();"     >> conftest.cc
    echo "int g() { SYMBOL<1> (); }" >>  conftest.cc

    dnl Compile it
    $CXX -c conftest.cc -o conftest.$ac_objext

    dnl and write symbols to file 2
    nm conftest.$ac_objext | grep SYMBOL > symb2

    dnl Compare the relevant symbols. According to the C++
    dnl standard, both source codes should produce the
    dnl same symbols.
    check="`diff symb1 symb2`"

    dnl Then try to link everything
    if test "x$check" = "x" ;
    then
        AC_MSG_RESULT(no)
    else
        AC_MSG_RESULT(yes)
        DARWIN_GCC_WEAK_LINKAGE_BUG="yes"
        AC_DEFINE(DEAL_II_WEAK_LINKAGE_BUG, 1,
                       [This error appears in the Apple edition of the
                        gcc 3.3, which ships with Darwin7.9.0 and
                        probably previous version. It leads to
                        problems during linking.
                        For the details, look at aclocal.m4 in the
                        top-level directory.])
    fi
    rm -f conftest.$ac_objext
    rm -f symb1 symb2
    ;;
  esac
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
dnl Old versions of gcc had a problem with arrays inside ?:
dnl expressions: they decayed too quickly to pointers. This then
dnl leads to erroneous warnings :-(
dnl
dnl Usage: DEAL_II_CHECK_ARRAY_CONDITIONAL_DECAY_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ARRAY_CONDITIONAL_DECAY_BUG, dnl
[
  AC_MSG_CHECKING(for array assignment in conditional bug)
  AC_LANG(C++)
  CXXFLAGS="-W -Wall -Werror"
  AC_TRY_COMPILE(
    [
    ],
    [
  const int x[2][2] = {{1,1},{1,1}};
  const int (&y)[2] = (1 ? x[0] : x[1]);
  return y[0];
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_ARRAY_CONDITIONAL_DECAY_BUG, 1,
                     [Defined if the compiler has a problem with
                      assigning arrays in conditionals])
    ])
])



dnl -------------------------------------------------------------
dnl Old versions of gcc (in particular gcc 3.3.3) has a problem
dnl if an argument to a member function is an array with bounds
dnl that depend on a static const variable inside that class.
dnl
dnl Usage: DEAL_II_CHECK_ARRAY_CONDITIONAL_DECAY_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ARRAY_ARG_BUG, dnl
[
  AC_MSG_CHECKING(for array argument bug)
  AC_LANG(C++)
  CXXFLAGS="-W -Wall -Werror"
  AC_TRY_COMPILE(
    [
      template <int dim> struct X {
        static const unsigned int N = 1<<dim;
        void f(int (&)[N]);
      };

      template <int dim> void X<dim>::f(int (&)[N]) {}
    ],
    [
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_ARRAY_ARG_BUG, 1,
                     [Defined if the compiler has a problem with
                      using arrays as arguments in functions])
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
dnl This is gcc bug 18644:
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=18644
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
dnl Some older versions of GCC (3.x series) have trouble with a
dnl certain part of the BOOST GRAPH library, yielding errors of
dnl the kind "sorry, unimplemented: use of `enumeral_type' in
dnl template type unification".
dnl
dnl Usage: DEAL_II_CHECK_BOOST_GRAPH_COMPILER_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BOOST_GRAPH_COMPILER_BUG, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for boost::graph compiler internal error)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSO $BOOST_INCLUDE_DIR"
    AC_TRY_COMPILE(
      [
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

namespace types
{
  using namespace boost;

  typedef adjacency_list<vecS, vecS, undirectedS,
                         property<vertex_color_t, default_color_type,
                                  property<vertex_degree_t,int> > > Graph;
}

void create_graph (types::Graph)
{}
      ],
      [
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(DEAL_II_BOOST_GRAPH_COMPILER_BUG, 1,
                  [Defined if the compiler gets an internal error compiling
                   some code that involves boost::graph])
      ])
  fi
])



dnl -------------------------------------------------------------
dnl The boost::shared_ptr class has a templated assignment operator
dnl but no assignment operator matching the default operator
dnl signature (this if for boost 1.29 at least). So when using
dnl using a normal assignment between identical types, the
dnl compiler synthesizes teh default operator, rather than using
dnl the template. It doesn't matter here, but is probably an
dnl oversight on behalf of the operators. It should be fixed in newer
dnl versions of boost.
dnl
dnl With -Wsynth in gcc we then get a warning. So if we find that
dnl this is still the case, disable -Wsynth, i.e. remove it from
dnl the list of warning flags.
dnl
dnl Usage: DEAL_II_CHECK_BOOST_SHARED_PTR_ASSIGNMENT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BOOST_SHARED_PTR_ASSIGNMENT, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for boost::shared_ptr assignment operator= template buglet)
    AC_LANG(C++)
    CXXFLAGS="-Wsynth -Werror"
    AC_TRY_COMPILE(
      [
#       include <boost/shared_ptr.hpp>
      ],
      [
        boost::shared_ptr<int> a,b;
        a = b;
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
dnl Some test, but for the icc C compiler.
dnl
dnl Usage: DEAL_II_ICC_C_WD_1572
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_ICC_C_WD_1572, dnl
[
  AC_MSG_CHECKING(whether -wd1572 is allowed for the C compiler)
  AC_LANG(C)
  OLDCFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -wd1572"
  AC_TRY_COMPILE( [], [],
      [
        AC_MSG_RESULT(yes)
        dnl Keep -wd1572 in CFLAGS
      ],
      [
        AC_MSG_RESULT(no)

        dnl Remove -wd1572 again from flags
        CFLAGS="$OLDCFLAGS"
      ])
])




dnl -------------------------------------------------------------
dnl gcc versions up to 2.95.3 had a problem with the std::advance function,
dnl when the number of steps forward was given by an unsigned number, since
dnl a comparison >=0 was performed on this number which leads to a warning
dnl that this comparison is always true. This is, at the best, annoying
dnl since it crops up at several places where std::advance is called from
dnl inside the library. Check whether the present version of the compiler
dnl has this problem
dnl
dnl This test is only called if gcc is used.
dnl
dnl Usage: DEAL_II_CHECK_ADVANCE_WARNING
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ADVANCE_WARNING, dnl
[
  AC_MSG_CHECKING(for std::advance warning)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -Werror"
  AC_TRY_COMPILE(
        [
#include <map>
#include <vector>
        ],
        [
  std::map<unsigned int, double> m;
  std::vector<std::pair<unsigned int, double> > v;
  v.insert (v.end(), m.begin(), m.end());
        ],
        [
          dnl compilation succeeded, no warning occured for the above code
          AC_MSG_RESULT(no)
          DEAL_II_ADVANCE_WARNING=no
        ],
        [
          AC_MSG_RESULT(yes)
          DEAL_II_ADVANCE_WARNING=yes
        ]
  )
])




dnl -------------------------------------------------------------
dnl On some Cygwin systems, a system header file includes this
dnl preprocessor define:
dnl   #define quad quad_t
dnl This is of course silly, but beyond that it also hurts as
dnl since we have member functions and variables with that name
dnl and we get compile errors depending or not we have this
dnl particular header file included.
dnl
dnl Fortunately, the define is only active is _POSIX_SOURCE is
dnl not set, so check for this define, and if necessary set
dnl this flag. We check on all systems, since maybe there are
dnl other such systems elsewhere...
dnl
dnl Usage: DEAL_II_CHECK_QUAD_DEFINE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_QUAD_DEFINE, dnl
[
  AC_MSG_CHECKING(for quad vs. quad_t define)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <sys/types.h>
#if defined(quad)
    no good system;
#endif
    ],
    [
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. working around)
      CXXFLAGSG="$CXXFLAGSG -D_POSIX_SOURCE"
      CXXFLAGSO="$CXXFLAGSO -D_POSIX_SOURCE"
    ])
])



dnl -------------------------------------------------------------
dnl On DEC OSF1, when we specify "-std strict_ansi", we can include
dnl errno.h, and still not get the definition of the error codes
dnl such as EINTR, EPIPE, etc.
dnl
dnl In this case, use a workaround by explicitly including
dnl /usr/include/errno.h, instead of just errno.h, which the compiler
dnl maps to one of its own C/C++ compatibility headers, which only
dnl define 3 error codes (for reasons unknown)
dnl
dnl Usage: DEAL_II_CHECK_ERROR_CODES_DEFINITION
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ERROR_CODES_DEFINITION, dnl
[
  AC_MSG_CHECKING(for definitions of error codes in errno.h)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <iostream>
#include <errno.h>
using namespace std;
    ],
    [
    cout << EINTR << endl;
    ],
    [
      AC_MSG_RESULT(yes)
    ],
    [
      AC_MSG_RESULT(no. working around)
      AC_DEFINE(DEAL_II_USE_DIRECT_ERRNO_H, 1,
                [Define if the compiler provides a <errno.g> header file
                 which does not define all error codes such as EINTR. In
                 that case, use the system include file at /usr/include
                 instead. There is probably a better way to do this, but
                 it is not apparent by looking at the C/C++ compatibility
                 header provided by the compiler.])
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

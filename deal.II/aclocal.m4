dnl ----------------------------------------------------
dnl        Platform specific autoconf tests
dnl
dnl    This is part of the input for the ./configure script of
dnl    the deal.II libraries. All options and paths are stored in
dnl    the file common/Make.global_options.
dnl
dnl    In doc/Makefile some information on the kind of documentation
dnl    is stored.
dnl
dnl
dnl author: Wolfgang Bangerth, 2000
dnl
dnl $Id$




dnl In some cases, -threads (or whatever else command line option)
dnl switches on some preprocessor flags. If this is not the case,
dnl then define them explicitely.
dnl 
dnl Usage: DEAL_II_THREAD_CPPFLAGS
dnl
AC_DEFUN(DEAL_II_THREAD_CPPFLAGS, dnl
[
  AC_MSG_CHECKING(for platform specific multi-threading defines)
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE(
   [
#if !defined (_REENTRANT) && !defined (_THREAD_SAFE)
#error Neither _REENTRANT nor _THREAD_SAFE were defined.
nonsense
#endif
   ],
   [
	;
   ],
   [
	AC_MSG_RESULT("not necessary")
   ],
   [
	AC_MSG_RESULT("-D_REENTRANT -D_THREAD_SAFE")
	CXXFLAGSG="$CXXFLAGSG -D_REENTRANT -D_THREAD_SAFE"
	CXXFLAGSO="$CXXFLAGSO -D_REENTRANT -D_THREAD_SAFE"
   ])
])




dnl Versions of gcc on different platforms use a multitude of flags to
dnl denote thread safe libraries and the like. They are, for example
dnl -threads/-pthread/-mthread, sometimes *thread, sometimes *threads. 
dnl Find out which is the right one on the present platform
dnl
dnl Usage: DEAL_II_FIND_THREAD_FLAGS
dnl
AC_DEFUN(DEAL_II_GET_THREAD_FLAGS, dnl
[
  AC_MSG_CHECKING(for platform specific thread flags)
  AC_REQUIRE([AC_LANG_CPLUSPLUS])
  for i in threads mt pthread pthreads mthreads Kthread kthread invalid_last_entry; do
    CXXFLAGS="$CXXFLAGSG -$i"
    DEAL_II_TRY_COMPILER_FLAG(
	[
	  thread_flag=$i
	  CXXFLAGSG="$CXXFLAGSG -$i"
	  CXXFLAGSO="$CXXFLAGSO -$i"
          LDFLAGS="$LDFLAGS -$i"
	
	  dnl The right parameter was found, so exit
	  break
   	])
  done
  if test "$thread_flag" = invalid_last_entry ; then
	AC_MSG_RESULT("no flag found!")
	AC_MSG_ERROR("Could not determine multithreading flag for this platform. Aborting!")
  fi
  AC_MSG_RESULT("-$thread_flag")
])



dnl Try whether the set of compiler flags in CXXFLAGS is reasonable, i.e.
dnl does not result in compiler messages (which are then produced by
dnl unknown or unrecognized compiler flags. This macro is mostly copied
dnl from the definition of AC_TRY_COMPILE, but the first two arguments are
dnl omitted and some other things are also simplified.
dnl
dnl Usage:
dnl   DEAL_II_TRY_COMPILER_FLAG([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
AC_DEFUN(DEAL_II_TRY_COMPILER_FLAG, dnl
[
  echo 'int main() { return 0; }' > conftest.$ac_ext

  dnl lets see whether the compiler generates output to stdout or stderr
  deal_II_compiler_output=`eval $ac_compile 2>&1`
  if test ! "$deal_II_compiler_output"; then
    ifelse([$1], , :, 
	   [
		rm -rf conftest*
    		$1
	   ])
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat conftest.$ac_ext >&AC_FD_CC
    ifelse([$2], , ,  
	   [
		rm -rf conftest*
		$2
           ])
  fi
  rm -f conftest*
])



dnl On SunOS 4.x, the `getrusage' function exists, but is not declared
dnl in the respective header file `resource.h', as one would think when
dnl reading the man pages. Then we have to declare this function 
dnl ourselves...
dnl
dnl If the function is not properly declared, then we augment the
dnl CXXFLAGS[OG] by `-DNO_HAVE_GETRUSAGE'
dnl
dnl Usage: DEAL_II_CHECK_GETRUSAGE
dnl
AC_DEFUN(DEAL_II_CHECK_GETRUSAGE, dnl
  AC_MSG_CHECKING(whether getrusage is properly declared)
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE(
    [
#include <sys/resource.h>
    ],
    [
      rusage *ru;
      getrusage(RUSAGE_SELF,ru);
    ],
    [
	AC_MSG_RESULT("yes")
    ],
    [
	AC_MSG_RESULT("no")
	CXXFLAGSG="$CXXFLAGSG -DNO_HAVE_GETRUSAGE"
	CXXFLAGSO="$CXXFLAGSO -DNO_HAVE_GETRUSAGE"
    ])
)      




dnl We'd like to use the `isnan' function. On some systems, this is
dnl simply declared in <math.h> (or <cmath>, for what it's worth), but
dnl on Linux for example, it is only declared if we specifically require
dnl support for ISO C 99. This macro checks whether `isnan' is declared
dnl or whether we have to pass special compiler flags, namely 
dnl -D_ISOC99_SOURCE (Linux) or -D__EXTENSIONS__ (Sun). Note that when
dnl checking we have to use the two sets of compiler flags.
dnl
dnl Usage: DEAL_II_CHECK_ISNAN
dnl
AC_DEFUN(DEAL_II_CHECK_ISNAN, dnl
  DEAL_II_CHECK_ISNAN_FLAG(debug, $CXXFLAGSG,
			   CXXFLAGSG="$deal_II_isnan_flag $CXXFLAGSG")
  DEAL_II_CHECK_ISNAN_FLAG(optimized, $CXXFLAGSO,
		 	   CXXFLAGSO="$deal_II_isnan_flag $CXXFLAGSO")
)


dnl The following function actually performs the check for the right flag
dnl for `isnan'. If a flag is found, the third argument is executed, and 
dnl the right flag is available in $deal_II_isnan_flag.
dnl
dnl Usage: DEAL_II_CHECK_ISNAN_FLAG("description of options set",
dnl                                 "compiler options set",
dnl                                 action when flag is found)
dnl
AC_DEFUN(DEAL_II_CHECK_ISNAN_FLAG, dnl
  AC_MSG_CHECKING(whether isnan is declared with $1 flags)
  AC_LANG_CPLUSPLUS
  CXXFLAGS=$2
  deal_II_isnan_flag=""
  AC_TRY_COMPILE(
    [
#include <cmath>
    ],
    [
	double d=0;
	isnan (d);
    ],
    [
	AC_MSG_RESULT("yes")
	deal_II_isnan_flag="-DHAVE_ISNAN"
	$3
    ])


  if test "x$deal_II_isnan_flag" = "x" ; then
    dnl Simply using isnan doesn't work. On Microsoft Windows systems, the function
    dnl is called _isnan, so check that
    AC_TRY_COMPILE(
      [
#include <cmath>
      ],
      [
	  double d=0;
	  _isnan (d);
      ],
      [
	  AC_MSG_RESULT("yes")
	  deal_II_isnan_flag="-DHAVE_UNDERSCORE_ISNAN"
	  $3
      ])
  fi


  dnl Let's see whether we _now_ have found something
  if test "x$deal_II_isnan_flag" = "x" ; then
    dnl We need to define something. Unfortunately, this is system 
    dnl dependent (argh!)
    deal_II_isnan_flag=""
    for testflag in -D_ISOC99_SOURCE -D__EXTENSIONS__ ; do 
      CXXFLAGS="$2 $testflag"
      AC_TRY_COMPILE(
        [
#include <cmath>
	],
	[
	  double d=0;
	  isnan (d);
	],
	[
	  dnl We found a flag by which isnan is defined; store
	  dnl this flag and exit the loop
	  deal_II_isnan_flag="-DHAVE_ISNAN $testflag"
	  break;
	])

      dnl If that didn't work (and it didn't as we are still inside the
      dnl loop), then try the _isnan function (maybe we are on a
      dnl Microsoft Windows system)
      AC_TRY_COMPILE(
        [
#include <cmath>
	],
	[
	  double d=0;
	  _isnan (d);
	],
	[
	  dnl We found a flag by which isnan is defined; store
	  dnl this flag and exit the loop
	  deal_II_isnan_flag="-DHAVE_UNDERSCORE_ISNAN $testflag"
	  break;
	])
    done

    dnl if no such flag was found, then abort ./configure since
    dnl the library will not be compilable on this platform
    dnl without knowledge of the right flag
    if test "x$deal_II_isnan_flag" = "x" ; then
      AC_MSG_RESULT(no.)
    else
      dnl we found something, lets us it
      AC_MSG_RESULT(using $testflag)
      $3
    fi
  fi
)      


dnl rand_r is defined for some compiler flag combinations, but not for
dnl others. check that. note that since these are C++ flags, we can't
dnl just use AC_CHECK_FUNCS
AC_DEFUN(DEAL_II_CHECK_RAND_R, dnl
  AC_LANG_CPLUSPLUS
  CXXFLAGS=$CXXFLAGSG
  AC_MSG_CHECKING(for rand_r)
  AC_TRY_COMPILE(
	[
#include <cstdlib>
	],
	[
int seed = 0;
int i=rand_r(&i);
	],
	[
	  AC_MSG_RESULT(found)
	  AC_DEFINE(HAVE_RAND_R)
	],
	[
	  AC_MSG_RESULT(no)
	]
  )
)


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
AC_DEFUN(DEAL_II_CHECK_ASSERT_THROW, dnl
  AC_MSG_CHECKING(whether AssertThrow works with $1 flags)
  AC_LANG_CPLUSPLUS
  CXXFLAGS=$2
  AC_TRY_COMPILE(
    [
#include "base/include/base/exceptions.h"
    ],
    [
	AssertThrow (false, ExcInternalError());
    ],
    [
	AC_MSG_RESULT("yes")
    ],
    [
	AC_MSG_RESULT("no")
	$3
    ])
)



dnl IBM xlC 5.0 from the VisualAge C++ pack has a bug with the following 
dnl code. We can work around it if we insert code like "using namespace std;"
dnl in the right place, but we'd like to do so only if absolutely necessary.
dnl Check whether the compiler which we are using has this bug.
dnl
dnl Usage: DEAL_II_CHECK_IBM_XLC_ERROR
dnl
AC_DEFUN(DEAL_II_CHECK_IBM_XLC_ERROR, dnl
  AC_MSG_CHECKING(for std::vector bug)
  AC_LANG_CPLUSPLUS
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      namespace std {
        template <class _Ty>                             class allocator {};
        template<class _Ty, class _Ax = allocator<_Ty> > class vector{};
      }

      struct X {};
      template <int dim> void g (const std::vector<X> &x);

      void f ()  {
        std::vector<X> x;
        g<1> (x);
      };
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. trying to work around)
      AC_DEFINE(XLC_WORK_AROUND_STD_BUG)
    ])
)


dnl gcc2.95 doesn't have the std::iterator class, but the standard requires it, so
dnl check whether we have to work around it
dnl
dnl Usage: DEAL_II_HAVE_STD_ITERATOR
dnl
AC_DEFUN(DEAL_II_HAVE_STD_ITERATOR, dnl
  AC_MSG_CHECKING(for std::iterator class)
  AC_LANG_CPLUSPLUS
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <iterator>
      class MyIterator : public std::iterator<std::bidirectional_iterator_tag,int>
      {};
    ],
    [],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_STD_ITERATOR_CLASS)
    ],
    [
      AC_MSG_RESULT(no)
    ])
)


dnl For gcc 2.95, when using the random_shuffle function with -ansi compiler flag,
dnl there is an error saying that lrand48 is undeclared. Fix that by declaring
dnl that function ourselves if necessary.
dnl
dnl Usage: DEAL_II_HAVE_LRAND48_DECLARED
dnl
AC_DEFUN(DEAL_II_HAVE_LRAND48_DECLARED, dnl
  AC_MSG_CHECKING(whether lrand48 needs to be declared with -ansi)
  AC_LANG_CPLUSPLUS
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <vector>
#include <algorithm>

void f()
{
  std::vector<unsigned int> new_indices;
  std::random_shuffle (new_indices.begin(), new_indices.end());
};
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_DECLARE_LRAND48)
    ])
)



dnl When compiling with ACE thread support, there are many constructs
dnl that are not allowed in C++, or that yield warnings when compiling with
dnl -ansi -pedantic. Check this, and if that is the case, set the variables
dnl $deal_II_ace_remove_ansi and $deal_II_ace_remove_pedantic to "yes".
dnl
dnl Usage: DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS
AC_DEFUN(DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS, dnl
  AC_MSG_CHECKING(whether compilation with ACE disallows flags)
  AC_LANG_CPLUSPLUS
  CXXFLAGS="-ansi -I$withmultithreading"
  AC_TRY_COMPILE(
    [
#  include <ace/Thread_Manager.h>
#  include <ace/Synch.h>
    ],
    [],
    [
      deal_II_ace_remove_ansi="no"
    ],
    [
      deal_II_ace_remove_ansi="yes"
    ])
  CXXFLAGS="-pedantic -Werror -I$withmultithreading"
  AC_TRY_COMPILE(
    [
#  include <ace/Thread_Manager.h>
#  include <ace/Synch.h>
    ],
    [],
    [
      deal_II_ace_remove_pedantic="no"
    ],
    [
      deal_II_ace_remove_pedantic="yes"
    ])

  if test $deal_II_ace_remove_ansi = "no" ; then
    if test $deal_II_ace_remove_pedantic = "no" ; then 
      AC_MSG_RESULT(no)
    else
      AC_MSG_RESULT(-pedantic)
    fi
  else
    if test $deal_II_ace_remove_pedantic = "no" ; then 
      AC_MSG_RESULT(-ansi)
    else
      AC_MSG_RESULT(-ansi -pedantic)
    fi
  fi
)

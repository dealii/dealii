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
  AC_REQUIRE([AC_LANG_CPLUSPLUS])
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
  for i in threads mt pthread pthreads mthreads Kthread kthread; do
    CXXFLAGS="$CXXFLAGSG -$i"
    AC_TRY_COMPILE(
     	[],
	[;],
	[
	  thread_flag=$i
	  CXXFLAGSG="$CXXFLAGSG -$i"
	  CXXFLAGSO="$CXXFLAGSO -$i"
	
	  dnl The right parameter was found, so exit
	  break
   	])
  done
  AC_MSG_RESULT("$thread_flag")
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
  AC_REQUIRE([AC_LANG_CPLUSPLUS])
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
dnl -D_ISOC99_SOURCE. Note that when checking we have to use the two
dnl sets of compiler flags.
dnl
dnl Usage: DEAL_II_CHECK_ISNAN
dnl
AC_DEFUN(DEAL_II_CHECK_ISNAN, dnl
  AC_MSG_CHECKING(whether isnan is declared with debug flags)
  AC_REQUIRE([AC_LANG_CPLUSPLUS])
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <cmath>
    ],
    [
	double d;
	isnan (d);
    ],
    [
	AC_MSG_RESULT("yes")
    ],
    [
	AC_MSG_RESULT("no")
	CXXFLAGSG="$CXXFLAGSG -D_ISOC99_SOURCE"
    ])
  AC_MSG_CHECKING(whether isnan is declared with optimized flags)
  AC_REQUIRE([AC_LANG_CPLUSPLUS])
  CXXFLAGS="$CXXFLAGSO"
  AC_TRY_COMPILE(
    [
#include <cmath>
    ],
    [
	double d;
	isnan (d);
    ],
    [
	AC_MSG_RESULT("yes")
    ],
    [
	AC_MSG_RESULT("no")
	CXXFLAGSO="$CXXFLAGSO -D_ISOC99_SOURCE"
    ])
)      

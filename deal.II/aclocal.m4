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
  for i in threads mt pthread pthreads mthreads Kthread kthread invalid_last_entry; do
    CXXFLAGS="$CXXFLAGSG -$i"
    DEAL_II_TRY_COMPILER_FLAG(
	[
	  thread_flag=$i
	  CXXFLAGSG="$CXXFLAGSG -$i"
	  CXXFLAGSO="$CXXFLAGSO -$i"
	
	  dnl The right parameter was found, so exit
	  break
   	])
  done
  if test $thread_flag = invalid_last_entry ; then
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
  cat > conftest.$ac_ext <<EOF
       int main() { return 0; }
EOF

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
  AC_REQUIRE([AC_LANG_CPLUSPLUS])
  CXXFLAGS=$2
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
    ],
    [
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
		deal_II_isnan_flag="-DHAVE_ISNAN $testflag"
		break;
	    ],
	    [
	])
	done

	dnl if no such flag was found, then abort ./configure since
	dnl the library will not be compilable on this platform
	dnl without knowledge of the right flag
	if test "$deal_II_isnan_flag" = "" ; then
	  AC_MSG_RESULT(no.)
	else

  	  dnl we found something, lets us it
	  AC_MSG_RESULT(using $testflag)
	  $3
	fi
    ])
)      

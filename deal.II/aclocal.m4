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

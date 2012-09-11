dnl -------------------------------------------------------------
dnl In some cases, -threads (or whatever else command line option)
dnl switches on some preprocessor flags. If this is not the case,
dnl then define them explicitely.
dnl
dnl Usage: DEAL_II_THREAD_CPPFLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_THREAD_CPPFLAGS, dnl
[
  AC_MSG_CHECKING(for platform specific multi-threading defines)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
   [
#       if !defined (_REENTRANT) && !defined (_THREAD_SAFE)
#       error Neither _REENTRANT nor _THREAD_SAFE were defined.
        nonsense
#       endif
   ],
   [
        ;
   ],
   [
        AC_MSG_RESULT(not necessary)
   ],
   [
        AC_MSG_RESULT(-D_REENTRANT -D_THREAD_SAFE)
        CXXFLAGSG="$CXXFLAGSG -D_REENTRANT -D_THREAD_SAFE"
        CXXFLAGSO="$CXXFLAGSO -D_REENTRANT -D_THREAD_SAFE"
   ])
])






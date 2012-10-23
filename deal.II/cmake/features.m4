dnl -------------------------------------------------------------
dnl Check for the Tecplot API. If it is found we will be able to write
dnl Tecplot binary files directly.
dnl
dnl This is a little ugly since we aren't guaranteed that TECHOME
dnl will point to the installation directory.  It could just as
dnl easily be TEC80HOME, TEC90HOME, etc...  So, better check them all
dnl
dnl Usage: DEAL_II_CONFIGURE_TECPLOT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_TECPLOT, dnl
[
  for i [ in $TECHOME $TEC100HOME $TEC90HOME $TEC80HOME ] ; do
    AC_CHECK_FILE($i/lib/tecio.a,
                  TECPLOT_LIB="$i/lib/tecio.a")
    AC_CHECK_FILE($i/include/TECIO.h,
                  TECPLOT_INCLUDE_DIR=-I$i/include,
                  TECPLOT_LIB="")
    if test "x$TECPLOT_LIB" != "x" ; then
      break
    fi
  done

  if (test -r "$TECPLOT_LIB") ; then
    AC_DEFINE(DEAL_II_HAVE_TECPLOT, 1,
              [Flag indicating whether the library shall be compiled to use the Tecplot interface])
    DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT($TECPLOT_LIB)
  fi
])

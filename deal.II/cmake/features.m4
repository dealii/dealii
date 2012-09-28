dnl ------------------------------------------------------------
dnl Check whether some of the HSL functions have been dropped
dnl into their respective place in the contrib subdir.
dnl Check for the following functions to be there:
dnl     MA27: needs files ma27.f
dnl     MA47: needs files ma47.f ma47.dep
dnl
dnl Usage: DEAL_II_CONFIGURE_HSL
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_HSL, dnl
[
  AC_MSG_CHECKING(for HSL subroutines)
  hsl_subroutines=""
  if test -r contrib/hsl/source/ma27.f ; then
    hsl_subroutines="$hsl_subroutines MA27"
    AC_DEFINE(HAVE_HSL_MA27, 1,
              [Availability of the MA27 algorithm from HSL])
  fi

  if (test -r contrib/hsl/source/ma47.f && \
      test -r contrib/hsl/source/ma47dep.f) ; then
    hsl_subroutines="$hsl_subroutines MA47"
    AC_DEFINE(HAVE_HSL_MA47, 1,
              [Availability of the MA47 algorithm from HSL])
  fi

  if test "x$hsl_subroutines" != "x" ; then
    AC_MSG_RESULT($hsl_subroutines)
    USE_CONTRIB_HSL=yes
    AC_SUBST(HSL_INCLUDE_DIR,'-I$D/contrib/hsl/include')
    AC_SUBST(NEEDS_F77LIBS,"yes")
  else
    AC_MSG_RESULT(none found)
  fi
])



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



dnl ------------------------------------------------------------
dnl Check whether SLEPc is installed, and if so store the
dnl respective links.
dnl
dnl Usage: DEAL_II_CONFIGURE_SLEPC
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_SLEPC, dnl
[
  AC_MSG_CHECKING([for SLEPc include directory])
  AC_ARG_WITH(slepc,
              [AS_HELP_STRING([--with-slepc=path/to/slepc],
              [Specify the path to the SLEPc installation, for which the include directory is a subdir; use this if you want to override the SLEPC_DIR environment variable.])],
              [dnl Special case when someone does --with-slepc=no
               if test "x$withval" = "xno" ; then
                 AC_MSG_RESULT([explicitly disabled])
                 USE_CONTRIB_SLEPC=no
               else
                 USE_CONTRIB_SLEPC=yes
                 DEAL_II_SLEPC_DIR="$withval"
                 AC_MSG_RESULT($DEAL_II_SLEPC_DIR)

                 dnl Make sure that what was specified is actually correct
               if test ! -d $DEAL_II_SLEPC_DIR \
                    -o ! -d $DEAL_II_SLEPC_DIR/include \
                    ; then
                 AC_MSG_ERROR([Path to SLEPc specified with --with-slepc does not point to a complete SLEPc installation])
               fi

               if test ! -d $DEAL_II_SLEPC_DIR/$DEAL_II_PETSC_ARCH \
                    -o ! -d $DEAL_II_SLEPC_DIR/$DEAL_II_PETSC_ARCH/lib \
                    ; then
            AC_MSG_ERROR([SLEPc has not been compiled for the PETSc architecture])
               fi
               fi
              ],
              [dnl Take something from the environment variables, if it is there
               if test "x$SLEPC_DIR" != "x" ; then
                 USE_CONTRIB_SLEPC=yes
                 DEAL_II_SLEPC_DIR="$SLEPC_DIR"
                 AC_MSG_RESULT($DEAL_II_SLEPC_DIR)

               dnl Make sure that what this is actually correct
               if test ! -d $DEAL_II_SLEPC_DIR \
                    -o ! -d $DEAL_II_SLEPC_DIR/include \
                    ; then
                 AC_MSG_ERROR([The path to SLEPc specified in the SLEPC_DIR environment variable does not point to a complete SLEPc installation])
               fi
               else
                 USE_CONTRIB_SLEPC=no
                 DEAL_II_SLEPC_DIR=""
                 AC_MSG_RESULT(not found)
               fi
              ])

  if test "$USE_CONTRIB_SLEPC" = "yes" ; then
    AC_DEFINE([DEAL_II_USE_SLEPC], [1],
              [Defined if a SLEPc installation was found and is going to be used])

    dnl Set an additional variable (not via AC_DEFINE, since we don't want
    dnl to have it in config.h) which we can use in doc/doxygen/options.dox.in.
    dnl If we have SLEPc, then the value of this variable expands to
    dnl defining the string "DEAL_II_USE_SLEPC" for the preprocessor. If
    dnl we don't have no SLEPc, then it does not define this string.
    DEAL_II_DEFINE_DEAL_II_USE_SLEPC=DEAL_II_USE_SLEPC
  fi

  dnl If we have found SLEPc, determine additional pieces of data
  if test "$USE_CONTRIB_SLEPC" = "yes" \
       ; then
    DEAL_II_CONFIGURE_SLEPC_VERSION

    dnl Finally set with_slepc if this hasn't happened yet
    if test "x$with_slepc" = "x" ; then
      with_slepc="yes"
  fi fi
])

dnl ------------------------------------------------------------
dnl Figure out the version numbers of SLEPc and compare with
dnl version numbers of PETSc. This is not strictly necessary
dnl but highly recommended that major, minor, and subminor
dnl version match. We blissfully ignor patch versions and hope
dnl for the best. If you want to overide all this you can.
dnl
dnl Usage: DEAL_II_CONFIGURE_SLEPC_VERSION
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_SLEPC_VERSION, dnl
[
  AC_MSG_CHECKING([for SLEPc version])
  DEAL_II_SLEPC_VERSION_MAJOR=`cat $DEAL_II_SLEPC_DIR/include/slepcversion.h \
                               | grep "#define SLEPC_VERSION_MAJOR" \
                               | perl -pi -e 's/.*MAJOR\s+//g;'`
  DEAL_II_SLEPC_VERSION_MINOR=`cat $DEAL_II_SLEPC_DIR/include/slepcversion.h \
                               | grep "#define SLEPC_VERSION_MINOR" \
                               | perl -pi -e 's/.*MINOR\s+//g;'`
  DEAL_II_SLEPC_VERSION_SUBMINOR=`cat $DEAL_II_SLEPC_DIR/include/slepcversion.h \
                               | grep "#define SLEPC_VERSION_SUBMINOR" \
                               | perl -pi -e 's/.*MINOR\s+//g;'`
  SLEPC_VERSION="$DEAL_II_SLEPC_VERSION_MAJOR.$DEAL_II_SLEPC_VERSION_MINOR.$DEAL_II_SLEPC_VERSION_SUBMINOR"

  dnl Here is where we check if the SLEPc version we have is a
  dnl release but do nothing about it.
  SLEPC_RELEASE=`cat $DEAL_II_SLEPC_DIR/include/slepcversion.h \
               | grep "#define SLEPC_VERSION_RELEASE" \
               | perl -pi -e 's/.*RELEASE\s+//g;'`
  if test "$SLEPC_RELEASE" = "0" ; then
    SLEPC_VERSION+="-dev"
  else
    SLEPC_VERSION+=""
  fi
  AC_MSG_RESULT($SLEPC_VERSION)

  dnl Then check that PETSc and SLEPc versions are compatible ie. that
  dnl they are equivalent. Patch numbers don't count for anything anymore,
  dnl but, we do include whether PETSc and SLEPc are both release
  dnl versions in the check. If they are not, we vomit.
  if test "${PETSC_VERSION}" != "${SLEPC_VERSION}" \
       -o "${PETSC_RELEASE}" != "${SLEPC_RELEASE}" \
       ; then
          AC_MSG_ERROR([If SLEPc is used, you must use the same version number as your PETSc Installation])
  fi
])



dnl ------------------------------------------------------------
dnl Check whether MUMPS is installed; and, if so, then check for
dnl two known dependecies, namely, SCALAPACK and BLACS.
dnl
dnl Usage: DEAL_II_CONFIGURE_MUMPS
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_MUMPS, dnl
[
  AC_MSG_CHECKING([for MUMPS library directory])
  AC_ARG_WITH(mumps,
    [AS_HELP_STRING([--with-mumps=path/to/mumps],
    [Specify the path to the MUMPS installation, for which the include directory and lib directory are subdirs; use this if you want to override the MUMPS_DIR environment variable.])],
    [dnl action-if-given
     if test "x$withval" = "xno" ; then
       AC_MSG_RESULT([explicitly disabled])
       USE_CONTRIB_MUMPS=no
     else
       USE_CONTRIB_MUMPS=yes
       DEAL_II_MUMPS_DIR="$withval"
       AC_MSG_RESULT($DEAL_II_MUMPS_DIR)
       dnl Make sure that what was specified is actually correct
       if test ! -d $DEAL_II_MUMPS_DIR         \
            -o ! -d $DEAL_II_MUMPS_DIR/include \
            -o ! -d $DEAL_II_MUMPS_DIR/lib     \
          ; then
         AC_MSG_ERROR([Path to MUMPS specified with --with-mumps does not point to a complete MUMPS installation])
       fi
     fi
    ],
    [dnl action-if-not-given (do nothing)
     USE_CONTRIB_MUMPS=no
     AC_MSG_RESULT([no])
    ])


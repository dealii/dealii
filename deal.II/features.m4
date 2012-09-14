dnl -------------------------------------------------------------
dnl Check for Doxygen.
dnl
dnl Usage: DEAL_II_CHECK_DOXYGEN
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_DOXYGEN, dnl
[
  AC_ARG_WITH(doxygen,
              [AS_HELP_STRING([--with-doxygen=filename],
              [Use 'filename' for doxygen.])],
      DOXYGEN="$withval",
      DOXYGEN=)

  dnl lets see whether the file exists
  if test "x$DOXYGEN" != "x" ; then
    AC_MSG_CHECKING(for specified doxygen path)
    if test -r $DOXYGEN ; then
      AC_MSG_RESULT($DOXYGEN)
    else
      AC_MSG_RESULT(not found)
      AC_MSG_ERROR(Invalid doxygen path $DOXYGEN)
    fi
  else
    dnl Check doxygen from the regular path. If we can't find it, then
    dnl set a flag and come back to that at the end of the ./configure
    dnl call.
    AC_PATH_PROG(DOXYGEN,doxygen)
    if test "x$DOXYGEN" = "x" ; then
      doxygen_not_found=yes;
    fi
  fi

  if test "x$doxygen_not_found" != "xyes" ; then
    AC_MSG_CHECKING(doxygen version)

    DOXYGEN_VERSION_STRING=`($DOXYGEN -v 2>&1) | grep "oxygen version"`
    case "$DOXYGEN_VERSION_STRING" in
      *1.3.* | *1.4.*)
        DOXYGEN_OPTIONS="options.136"
        AC_MSG_RESULT(pre 1.5)
        ;;
      *)
        DOXYGEN_OPTIONS="options.dox"
        AC_MSG_RESULT(1.5.x or later)
        ;;
    esac
  fi

  dnl Doxygen needs 'dot' for inheritance graph generation
  DEAL_II_CHECK_DOT
])



dnl -------------------------------------------------------------
dnl Check for DOT.
dnl
dnl Usage: DEAL_II_CHECK_DOT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_DOT, dnl
[
  AC_CHECK_PROG(DOT,dot,dot)
  if test "x$DOT" = "x" ; then
    DEAL_II_HAVE_DOT=NO;
  else
    DEAL_II_HAVE_DOT=YES;
  fi
])



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
dnl Check whether PETSc is installed, and if so store the
dnl respective links
dnl
dnl Usage: DEAL_II_CONFIGURE_PETSC
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_PETSC, dnl
[
  dnl First check for the PETSc directory
  AC_MSG_CHECKING([for PETSc library directory])

  AC_ARG_WITH(petsc,
              [AS_HELP_STRING([--with-petsc=path/to/petsc],
              [Specify the path to the PETSc installation, of which the include and library directories are subdirs; use this if you want to override the PETSC_DIR environment variable.])],
     [
        dnl Special case when someone does --with-petsc=no
        if test "x$withval" = "xno" ; then
          AC_MSG_RESULT([explicitly disabled])
          USE_CONTRIB_PETSC=no
        else
          USE_CONTRIB_PETSC=yes
          DEAL_II_PETSC_DIR="$withval"
          AC_MSG_RESULT($DEAL_II_PETSC_DIR)

          dnl Make sure that what was specified is actually correct
          if test ! -d $DEAL_II_PETSC_DIR/include \
               ; then
            AC_MSG_ERROR([Path to PETSc specified with --with-petsc does not point to a complete PETSc installation])
          fi
        fi
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$PETSC_DIR" != "x" ; then
          USE_CONTRIB_PETSC=yes
          DEAL_II_PETSC_DIR="$PETSC_DIR"
          AC_MSG_RESULT($DEAL_II_PETSC_DIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_PETSC_DIR/include \
               ; then
            AC_MSG_ERROR([The path to PETSc specified in the PETSC_DIR environment variable does not point to a complete PETSc installation])
          fi
        else
          USE_CONTRIB_PETSC=no
          DEAL_II_PETSC_DIR=""
          AC_MSG_RESULT([not found])
        fi
     ])
  if test "$USE_CONTRIB_PETSC" = "yes" ; then
    AC_DEFINE([DEAL_II_USE_PETSC], [1],
              [Defined if a PETSc installation was found and is going to be used])

    dnl Set an additional variable (not via AC_DEFINE, since we don't want
    dnl to have it in config.h) which we can use in doc/doxygen/options.dox.in.
    dnl If we have PETSc, then the value of this variable expands to
    dnl defining the string "DEAL_II_USE_PETSC" for the preprocessor. If
    dnl we don't have no PETSc, then it does not define this string.
    DEAL_II_DEFINE_DEAL_II_USE_PETSC=DEAL_II_USE_PETSC

    dnl Also work around a stupidity in PETSc that makes sure it interferes in
    dnl a completely obnoxious way with boost.
    AC_DEFINE([PETSC_SKIP_UNDERSCORE_CHKERR], [1],
              [Make sure PETSc doesn't re-define the underscore through the preprocessor, since this interferes with boost. PETSc redefines the underscore to be "__gterr =", but then forgets to undef this thing. Boost simply wants to concatenate the underscore with another string to form a class name, which then of course isn't valid any more. See mails in early Feb 2006.])
  fi

  dnl If we have found PETSc, determine additional pieces of data
  if test "$USE_CONTRIB_PETSC" = "yes" ; then
    DEAL_II_CONFIGURE_PETSC_VERSION
    DEAL_II_CONFIGURE_PETSC_ARCH

    DEAL_II_CONFIGURE_PETSC_MPIUNI_LIB
    DEAL_II_CHECK_PETSC_MPI_CONSISTENCY
    DEAL_II_CONFIGURE_PETSC_COMPLEX

    DEAL_II_EXPAND_PETSC_VECTOR="PETScWrappers::Vector"
    DEAL_II_EXPAND_PETSC_MPI_VECTOR="PETScWrappers::MPI::Vector"
    DEAL_II_EXPAND_PETSC_BLOCKVECTOR="PETScWrappers::BlockVector"
    DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR="PETScWrappers::MPI::BlockVector"

    dnl Finally set with_petsc if this hasn't happened yet
    if test "x$with_petsc" = "x" ; then
      with_petsc="yes"
    fi
  fi

  dnl Make sure that the right values for PETSC vectors are written
  dnl into common/template-arguments.in
  AC_SUBST(DEAL_II_EXPAND_PETSC_VECTOR)
  AC_SUBST(DEAL_II_EXPAND_PETSC_MPI_VECTOR)
  AC_SUBST(DEAL_II_EXPAND_PETSC_BLOCKVECTOR)
  AC_SUBST(DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR)
])



dnl ------------------------------------------------------------
dnl Figure out the architecture used for PETSc, since that
dnl determines where object and configuration files will be found.
dnl
dnl Usage: DEAL_II_CONFIGURE_PETSC_ARCH
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_PETSC_ARCH, dnl
[
  AC_MSG_CHECKING([for PETSc library architecture])

  AC_ARG_WITH(petsc-arch,
              [AS_HELP_STRING([--with-petsc-arch=architecture],
              [Specify the architecture for your PETSc installation; use this if you want to override the PETSC_ARCH environment variable.])],
              [DEAL_II_PETSC_ARCH="$withval"
               AC_MSG_RESULT($DEAL_II_PETSC_ARCH)
              ],
              [dnl Take something from the environment variables
               if test "x$PETSC_ARCH" != "x" ; then
                 DEAL_II_PETSC_ARCH="$PETSC_ARCH"
                 AC_MSG_RESULT($DEAL_II_PETSC_ARCH)
               else
                 AC_MSG_ERROR([If PETSc is used, you must specify the architecture either through the PETSC_ARCH environment variable, or through the --with-petsc-arch flag])
               fi
              ])

  if test "x$PETSC_ARCH" != "x" ; then

    dnl PETSc change the locations where they store their libraries
    dnl from time-to-time; so make sure that what was specified is
    dnl actually correct.
    case "${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}" in
      3.*) dnl
        if test ! -d $DEAL_II_PETSC_DIR/$DEAL_II_PETSC_ARCH/lib \
           ; then
          AC_MSG_ERROR([PETSc has not been compiled for the architecture specified with --with-petsc-arch])
        fi
        ;;
      *)    dnl
        AC_MSG_ERROR([Unknown PETSc version ${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}])
        ;;
    esac
  fi
])



dnl ------------------------------------------------------------
dnl Figure out the version numbers of PETSc. This is unfortunately
dnl necessary since PETSc has a habit to change function signatures,
dnl library names, etc, in random ways between versions...
dnl
dnl Usage: DEAL_II_CONFIGURE_PETSC_VERSION
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_PETSC_VERSION, dnl
[
  AC_MSG_CHECKING([for PETSc version])
  DEAL_II_PETSC_VERSION_MAJOR=`cat $DEAL_II_PETSC_DIR/include/petscversion.h \
                               | grep "#define PETSC_VERSION_MAJOR" \
                               | perl -pi -e 's/.*MAJOR\s+//g;'`
  DEAL_II_PETSC_VERSION_MINOR=`cat $DEAL_II_PETSC_DIR/include/petscversion.h \
                               | grep "#define PETSC_VERSION_MINOR" \
                               | perl -pi -e 's/.*MINOR\s+//g;'`
  DEAL_II_PETSC_VERSION_SUBMINOR=`cat $DEAL_II_PETSC_DIR/include/petscversion.h \
                               | grep "#define PETSC_VERSION_SUBMINOR" \
                               | perl -pi -e 's/.*MINOR\s+//g;'`
  PETSC_VERSION="$DEAL_II_PETSC_VERSION_MAJOR.$DEAL_II_PETSC_VERSION_MINOR.$DEAL_II_PETSC_VERSION_SUBMINOR"

  dnl Here is where we check if the PETSc version we have is a
  dnl release but do nothing about it.
  PETSC_RELEASE=`cat $DEAL_II_PETSC_DIR/include/petscversion.h \
                     | grep "#define PETSC_VERSION_RELEASE" \
                     | perl -pi -e 's/.*RELEASE\s+//g;'`
  if test "$PETSC_RELEASE" = "0" ; then
    PETSC_VERSION+="-dev"
    DEAL_II_PETSC_VERSION_DEV=yes
    AC_DEFINE([DEAL_II_USE_PETSC_DEV], [1],
              [Defined if a PETSc installation was found and is not a release])
  else
    DEAL_II_PETSC_VERSION_DEV=no
  fi

  AC_MSG_RESULT($PETSC_VERSION)
])


dnl -------------------------------------------------------------
dnl Make sure that if PETSc and deal.II were built with the same
dnl MPI enabled (or disabled) functionality.
dnl
dnl Usage: DEAL_II_CHECK_PETSC_MPI_CONSISTENCY
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_PETSC_MPI_CONSISTENCY, dnl
[
  dnl Then check for MPI consistency.
  AC_MSG_CHECKING(for consistency of PETSc and deal.II MPI settings)

  if test "x$DEAL_II_USE_MPI" = "xyes" ; then

  dnl So we support MPI. Check that our PETSc installation does
  dnl too. PETSc sets the variable PETSC_HAVE_MPIUNI to 1 in case
  dnl he *does not* support MPI, so just read out that information.
  dnl
  dnl Like always, we have to cake care of version control!
    case "${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}" in
      3.*) dnl
        AC_TRY_COMPILE(
        [#include "$DEAL_II_PETSC_DIR/$DEAL_II_PETSC_ARCH/include/petscconf.h"
        ],
        [#ifdef PETSC_HAVE_MPIUNI
           compile error;
         #endif
        ],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR([PETSc was not built for MPI, but deal.II is!]
        )])
      ;;
      *) dnl
        AC_MSG_ERROR([Unknown PETSc version ${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}])
      ;;
    esac
  else
    case "${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}" in
      3.*) dnl
        AC_TRY_COMPILE(
        [#include "$DEAL_II_PETSC_DIR/$DEAL_II_PETSC_ARCH/include/petscconf.h"
        ],
        [#ifndef PETSC_HAVE_MPIUNI
           compile error;
         #endif],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR([PETSc was built for MPI, but deal.II is not!])
        ])
      ;;
      *) dnl
        AC_MSG_ERROR([Unknown PETSc version ${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}])
      ;;
    esac
  fi
])


dnl ------------------------------------------------------------
dnl See if there is a library libmpiuni.a/so available. We need
dnl to link with it on some systems where PETSc is built without
dnl a real MPI and we need to handle trivial (one process) MPI
dnl functionality.
dnl
dnl Usage: DEAL_II_CONFIGURE_PETSC_MPIUNI_LIB
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_PETSC_MPIUNI_LIB, dnl
[
  AC_MSG_CHECKING([for PETSc libmpiuni library])
  case "${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}" in
    3.*) dnl
      if test -f $DEAL_II_PETSC_DIR/$DEAL_II_PETSC_ARCH/lib/libmpiuni.a ; then
        DEAL_II_PETSC_MPIUNI_LIB="$DEAL_II_PETSC_DIR/$DEAL_II_PETSC_ARCH/lib/libmpiuni.a"
      fi
      ;;
    *)    dnl
      AC_MSG_ERROR([Unknown PETSc version ${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}])
    ;;
  esac

  if test "$DEAL_II_PETSC_MPIUNI_LIB" = "" ; then
     AC_MSG_RESULT([not found])
  else
     AC_MSG_RESULT($DEAL_II_PETSC_MPIUNI_LIB)
  fi
])


dnl ------------------------------------------------------------
dnl Figure out PETSc was compiled with --scalar-type=complex by
dnl scanning PETSc configuration file.
dnl
dnl Warning: Up tp now, PETSc>3.0.0 is being supported and
dnl deal.II will not safely compile if this option is "yes".
dnl
dnl Usage: DEAL_II_CONFIGURE_PETSC_COMPLEX
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_PETSC_COMPLEX, dnl
[
  case "${DEAL_II_PETSC_VERSION_MAJOR}.${DEAL_II_PETSC_VERSION_MINOR}.${DEAL_II_PETSC_VERSION_SUBMINOR}" in
    3.3*)
      AC_MSG_CHECKING([for PETSc scalar complex])
      DEAL_II_PETSC_COMPLEX=`cat $DEAL_II_PETSC_DIR/$DEAL_II_PETSC_ARCH/include/petscconf.h \
                               | grep "#define PETSC_USE_COMPLEX" \
                               | perl -pi -e 's/.*COMPLEX\s+//g;'`
      if test "$DEAL_II_PETSC_COMPLEX" = "1" ; then
         AC_MSG_RESULT(yes)
      else
         AC_MSG_RESULT(no)
      fi

      dnl If we have previously found PETSc and here with a complex
      dnl scalar type then set the DEAL_II_USE_COMPLEX macro
      if test "$USE_CONTRIB_PETSC" = "yes" ; then
        if test "$DEAL_II_PETSC_COMPLEX" = "1" ; then
        AC_DEFINE([DEAL_II_USE_PETSC_COMPLEX], [1],
                  [Defined if a PETSc installation was found with complex scalar type and is going to be used])
      fi fi
    ;;
  esac
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
dnl Check whether Trilinos is installed, and if so store the
dnl respective links
dnl
dnl Usage: DEAL_II_CONFIGURE_TRILINOS
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_TRILINOS, dnl
[
 AC_MSG_CHECKING(for Trilinos directory)

  AC_ARG_WITH(trilinos,
              [AS_HELP_STRING([--with-trilinos=/path/to/trilinos],
              [Specify the path to the Trilinos installation, of which the include and lib directories are subdirs; use this if you want to override the TRILINOS_DIR environment variable.])],
     [
        dnl Special case when someone does --with-trilinos=no
        if test "x$withval" = "xno" ; then
          AC_MSG_RESULT([explicitly disabled])
          USE_CONTRIB_TRILINOS=no
        else
          USE_CONTRIB_TRILINOS=yes
          DEAL_II_TRILINOS_DIR="$withval"
          AC_MSG_RESULT($DEAL_II_TRILINOS_DIR)

          dnl Make sure that what was specified is actually correct
          if test ! -d $DEAL_II_TRILINOS_DIR/include \
               -o ! -d $DEAL_II_TRILINOS_DIR/lib ; then
            AC_MSG_ERROR([Path to Trilinos specified with --with-trilinos does not point to a complete Trilinos installation])
          fi

          DEAL_II_TRILINOS_INCDIR="$DEAL_II_TRILINOS_DIR/include"
          DEAL_II_TRILINOS_LIBDIR="$DEAL_II_TRILINOS_DIR/lib"
        fi
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$TRILINOS_DIR" != "x" ; then
          dnl Special case when someone does --with-trilinos=no
          if test "x$withval" = "xno" ; then
            AC_MSG_RESULT([explicitly disabled])
            USE_CONTRIB_TRILINOS=no
          else
            USE_CONTRIB_TRILINOS=yes
            DEAL_II_TRILINOS_DIR="$TRILINOS_DIR"
            AC_MSG_RESULT($DEAL_II_TRILINOS_DIR)

            dnl Make sure that what this is actually correct
            if test ! -d $DEAL_II_TRILINOS_DIR/include \
                 -o ! -d $DEAL_II_TRILINOS_DIR/lib ; then
              AC_MSG_ERROR([The path to Trilinos specified in the TRILINOS_DIR environment variable does not point to a complete Trilinos installation])
            fi
            DEAL_II_TRILINOS_INCDIR="$DEAL_II_TRILINOS_DIR/include"
            DEAL_II_TRILINOS_LIBDIR="$DEAL_II_TRILINOS_DIR/lib"
          fi
        else
          USE_CONTRIB_TRILINOS=no
          DEAL_II_TRILINOS_DIR=""
          AC_MSG_RESULT(not found)
        fi
     ])

  AC_MSG_CHECKING(for Trilinos header directory)

  AC_ARG_WITH(trilinos-include,
              [AS_HELP_STRING([--with-trilinos-include=/path/to/trilinos],
              [Specify the path to the Trilinos include; use this if you want to override the TRILINOS_INCDIR environment variable.])],
     [
        dnl Special case when someone does --with-trilinos=no
        if test "x$withval" = "xno" ; then
          AC_MSG_RESULT([explicitly disabled])
          USE_CONTRIB_TRILINOS=no
        else
          USE_CONTRIB_TRILINOS=yes
          DEAL_II_TRILINOS_INCDIR="$withval"
          AC_MSG_RESULT($DEAL_II_TRILINOS_INCDIR)

          dnl Make sure that what was specified is actually correct
          if test ! -d $DEAL_II_TRILINOS_INCDIR ; then
            AC_MSG_ERROR([Path to Trilinos specified with --with-trilinos-include does not point to a complete Trilinos installation])
          fi
        fi
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$TRILINOS_INCDIR" != "x" ; then
          USE_CONTRIB_TRILINOS=yes
          DEAL_II_TRILINOS_INCDIR="$TRILINOS_INCDIR"
          AC_MSG_RESULT($DEAL_II_TRILINOS_INCDIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_TRILINOS_INCDIR ; then
            AC_MSG_ERROR([The path to Trilinos includes specified in the TRILINOS_INCDIR environment variable does not point to a valid directory])
          fi
        else
          dnl --with-trilinos-include not explicitly specified. do
          dnl nothing if --with-trilinos has previously been specified,
          dnl otherwise say no to trilinos
          if test "x${USE_CONTRIB_TRILINOS}" != "xyes" ; then
            USE_CONTRIB_TRILINOS=no
            DEAL_II_TRILINOS_INCDIR=""
            AC_MSG_RESULT(not found)
          else
            AC_MSG_RESULT(not explicitly specified)
          fi
        fi
     ])

  AC_MSG_CHECKING(for Trilinos library directory)

  AC_ARG_WITH(trilinos-libs,
              [AS_HELP_STRING([--with-trilinos-libs=/path/to/trilinos],
              [Specify the path to the Trilinos libraries; use this if you want to override the TRILINOS_LIBDIR environment variable.])],
     [
        dnl Special case when someone does --with-trilinos=no
        if test "x$withval" = "xno" ; then
          AC_MSG_RESULT([explicitly disabled])
          USE_CONTRIB_TRILINOS=no
        else
          USE_CONTRIB_TRILINOS=yes
          DEAL_II_TRILINOS_LIBDIR="$withval"
          AC_MSG_RESULT($DEAL_II_TRILINOS_LIBDIR)

          dnl Make sure that what was specified is actually correct
          if test ! -d $DEAL_II_TRILINOS_LIBDIR ; then
            AC_MSG_ERROR([Path to Trilinos libraries with --with-trilinos-libs does not point to a valid directory])
          fi
        fi
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$TRILINOS_LIBDIR" != "x" ; then
          USE_CONTRIB_TRILINOS=yes
          DEAL_II_TRILINOS_LIBDIR="$TRILINOS_LIBDIR"
          AC_MSG_RESULT($DEAL_II_TRILINOS_LIBDIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_TRILINOS_LIBDIR ; then
            AC_MSG_ERROR([The path to Trilinos specified in the TRILINOS_LIBDIR environment variable does not point to a complete Trilinos installation])
          fi
        else
          dnl --with-trilinos-libs not explicitly specified. do
          dnl nothing if --with-trilinos has previously been specified,
          dnl otherwise say no to trilinos
          if test "x${USE_CONTRIB_TRILINOS}" != "xyes" ; then
            USE_CONTRIB_TRILINOS=no
            DEAL_II_TRILINOS_LIBDIR=""
            AC_MSG_RESULT(not found)
          else
            AC_MSG_RESULT(not explicitly specified)
          fi
        fi
     ])

  dnl If we have found Trilinos, determine and set additional pieces of data
  if test "$USE_CONTRIB_TRILINOS" = "yes" ; then
    AC_DEFINE(DEAL_II_USE_TRILINOS, 1,
              [Defined if a Trilinos installation was found and is going
               to be used])

    dnl Set an additional variable (not via AC_DEFINE, since we don't want
    dnl to have it in config.h) which we can use in doc/doxygen/options.dox.in.
    dnl If we have Trilinos, then the value of this variable expands to
    dnl defining the string "DEAL_II_USE_TRILINOS" for the preprocessor. If
    dnl we don't have no Trilinos, then it does not define this string.
    DEAL_II_DEFINE_DEAL_II_USE_TRILINOS=DEAL_II_USE_TRILINOS

    DEAL_II_CONFIGURE_TRILINOS_VERSION
    DEAL_II_CHECK_TRILINOS_MPI_CONSISTENCY
    DEAL_II_CHECK_TRILINOS_SHARED_STATIC
    DEAL_II_CHECK_TRILINOS_LIBS
    DEAL_II_CHECK_TRILINOS_WARNINGS
    DEAL_II_CHECK_TRILINOS_HEADER_FILES

    DEAL_II_EXPAND_TRILINOS_VECTOR="TrilinosWrappers::Vector"
    DEAL_II_EXPAND_TRILINOS_MPI_VECTOR="TrilinosWrappers::MPI::Vector"
    DEAL_II_EXPAND_TRILINOS_BLOCKVECTOR="TrilinosWrappers::BlockVector"
    DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR="TrilinosWrappers::MPI::BlockVector"
    DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN="TrilinosWrappers::SparsityPattern"
    DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN="TrilinosWrappers::BlockSparsityPattern"

    dnl Finally set with_trilinos if this hasn't happened yet
    if test "x$with_trilinos" = "x" ; then
      with_trilinos="yes"
    fi
  fi

  dnl Make sure that the right values for Trilinos vectors are written into
  dnl common/template-arguments.in
  AC_SUBST(DEAL_II_EXPAND_TRILINOS_VECTOR)
  AC_SUBST(DEAL_II_EXPAND_TRILINOS_MPI_VECTOR)
  AC_SUBST(DEAL_II_EXPAND_TRILINOS_BLOCKVECTOR)
  AC_SUBST(DEAL_II_EXPAND_TRILINOS_MPI_BLOCKVECTOR)
  AC_SUBST(DEAL_II_EXPAND_TRILINOS_SPARSITY_PATTERN)
  AC_SUBST(DEAL_II_EXPAND_TRILINOS_BLOCK_SPARSITY_PATTERN)
])



dnl ------------------------------------------------------------
dnl Figure out the version numbers of Trilinos.
dnl
dnl Usage: DEAL_II_CONFIGURE_TRILINOS_VERSION
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_TRILINOS_VERSION, dnl
[
  AC_MSG_CHECKING([for Trilinos version])
  DEAL_II_TRILINOS_VERSION_MAJOR=`cat $DEAL_II_TRILINOS_INCDIR/Trilinos_version.h \
                               | grep "#define TRILINOS_MAJOR_VERSION" \
                               | perl -pi -e 's/.*VERSION\s+//g;'`
  DEAL_II_TRILINOS_VERSION_MINOR=`cat $DEAL_II_TRILINOS_INCDIR/Trilinos_version.h \
                               | grep "#define TRILINOS_MAJOR_MINOR_VERSION" \
                               | perl -pi -e 's/.*VERSION\s+\d?\d(\d\d)\d\d/\1/g;' \
                               | perl -pi -e 's/0(\d)/\1/g;'`
  DEAL_II_TRILINOS_VERSION_SUBMINOR=`cat $DEAL_II_TRILINOS_INCDIR/Trilinos_version.h \
                               | grep "#define TRILINOS_MAJOR_MINOR_VERSION" \
                               | perl -pi -e 's/.*VERSION\s+\d?\d\d\d(\d\d)/\1/g;' \
                               | perl -pi -e 's/0(\d)/\1/g;'`
  AC_MSG_RESULT([$DEAL_II_TRILINOS_VERSION_MAJOR.$DEAL_II_TRILINOS_VERSION_MINOR.$DEAL_II_TRILINOS_VERSION_SUBMINOR])

  dnl Verify that we have at least Trilinos 10. This is the
  dnl version where Trilinos started using cmake, which allow
  dnl us to figure out which libraries Trilinos has built
  dnl and in which order they need to be linked
  if test "$DEAL_II_TRILINOS_VERSION_MAJOR" -lt 10 ; then
    AC_MSG_ERROR([Trilinos versions prior to 10.0 are no longer supported with deal.II.])
  fi

  dnl Trilinos 10.6 had quite a number of bugs we ran into, see
  dnl for example
  dnl   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5062
  dnl   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5319
  dnl The same is unfortunately true for 10.8.[01]:
  dnl   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5370
  if test "$DEAL_II_TRILINOS_VERSION_MAJOR" = 10 -a "$DEAL_II_TRILINOS_VERSION_MINOR" = 6 ; then
    AC_MSG_ERROR([Trilinos versions 10.6.x have bugs that make it incompatible
                  with deal.II. Please use versions before 10.6 or after 10.8.1.])
  fi
  if test "$DEAL_II_TRILINOS_VERSION_MAJOR" = 10 \
       -a "$DEAL_II_TRILINOS_VERSION_MINOR" = 8  \
       -a "$DEAL_II_TRILINOS_VERSION_SUBMINOR" -lt 2 ; then
    AC_MSG_ERROR([Trilinos versions 10.8.0 and 10.8.1 have bugs that make it incompatible
                  with deal.II. Please use versions before 10.6 or after 10.8.1.])
  fi

  AC_SUBST(DEAL_II_TRILINOS_VERSION_MAJOR)
  AC_SUBST(DEAL_II_TRILINOS_VERSION_MINOR)
  AC_SUBST(DEAL_II_TRILINOS_VERSION_SUBMINOR)
  AC_SUBST(DEAL_II_TRILINOS_LIBPREFIX)
])



dnl -------------------------------------------------------------
dnl Make sure that if Trilinos was built with/without MPI, then
dnl deal.II was built with the same flags.
dnl
dnl Usage: DEAL_II_CHECK_TRILINOS_MPI_CONSISTENCY
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TRILINOS_MPI_CONSISTENCY, dnl
[
  dnl Check for presence of Epetra_config.h that we need to detect MPI
  dnl settings
  AC_MSG_CHECKING(Epetra_config.h presence)
  if test -f $DEAL_II_TRILINOS_INCDIR/Epetra_config.h ; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    exit 1;
  fi

  AC_MSG_CHECKING(for consistency of Trilinos and deal.II MPI settings)
  AC_LANG(C++)

  OLD_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS -I$DEAL_II_TRILINOS_INCDIR"

  dnl Trilinos Epetra's Epetra_config.h might provide
  dnl   #define PACKAGE_BUGREPORT
  dnl   #define PACKAGE_NAME
  dnl   #define PACKAGE_STRING
  dnl   #define PACKAGE_TARNAME
  dnl   #define PACKAGE_VERSION
  dnl which is already set for the deal.II package. So undefine them for
  dnl this test.
  cp confdefs.h confdefs.h.bak
  echo "#undef PACKAGE_BUGREPORT" >> confdefs.h
  echo "#undef PACKAGE_NAME" >> confdefs.h
  echo "#undef PACKAGE_STRING" >> confdefs.h
  echo "#undef PACKAGE_TARNAME" >> confdefs.h
  echo "#undef PACKAGE_VERSION" >> confdefs.h

  if test "x$DEAL_II_USE_MPI" = "xyes" ; then
    dnl So we support MPI. Check that our Trilinos installation
    dnl does too. Epetra sets the variable HAVE_MPI to 1 in case
    dnl supports MPI, and does not set it otherwise, so just read
    dnl out that information.
    AC_TRY_COMPILE(
      [
        #include <Epetra_config.h>
      ],
      [
        #ifndef HAVE_MPI
          compile error;
        #endif
      ],
      [
        AC_MSG_RESULT(yes)
      ],
      [
        AC_MSG_ERROR([Trilinos was not built for MPI, but deal.II is!])
        exit 1;
      ])
  else
    dnl So we don't support MPI. Check that our Trilinos installation
    dnl doesn't either.
    AC_TRY_COMPILE(
      [
        #include <Epetra_config.h>
      ],
      [
        #ifdef HAVE_MPI
          compile error;
        #endif
      ],
      [
        AC_MSG_RESULT(yes)
      ],
      [
        AC_MSG_ERROR([Trilinos was built for MPI, but deal.II is not!])
        exit 1;
      ])
  fi

  mv confdefs.h.bak confdefs.h
  CXXFLAGS="${OLD_CXXFLAGS}"
])



dnl ------------------------------------------------------------
dnl Check whether the installed version of Trilinos uses shared
dnl or static libs, or both. Produce an error if this doesn't
dnl match the kind of libraries we produce here
dnl
dnl Usage: DEAL_II_CHECK_TRILINOS_SHARED_STATIC
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TRILINOS_SHARED_STATIC, dnl
[
  dnl Check using the epetra library since that should always be there
  dnl The problem is that on Debian, the library isn't called
  dnl libepetra.so, but instead libtrilinos_epetra.so, so
  dnl record this prefix in a variable of its own
  AC_MSG_CHECKING(whether Trilinos uses shared libraries)
  if test -f $DEAL_II_TRILINOS_LIBDIR/libepetra${shared_lib_suffix} ; then
    AC_MSG_RESULT(yes)
    DEAL_II_TRILINOS_SHARED=yes
  elif test -f $DEAL_II_TRILINOS_LIBDIR/libtrilinos_epetra${shared_lib_suffix} ; then
    AC_MSG_RESULT(yes)
    DEAL_II_TRILINOS_SHARED=yes
    DEAL_II_TRILINOS_LIBPREFIX="trilinos_"
  else
    AC_MSG_RESULT(no)
  fi

  AC_MSG_CHECKING(whether Trilinos uses static libraries)
  if test -f $DEAL_II_TRILINOS_LIBDIR/libepetra${static_lib_suffix} ; then
    AC_MSG_RESULT(yes)
    DEAL_II_TRILINOS_STATIC=yes
  elif test -f $DEAL_II_TRILINOS_LIBDIR/libtrilinos_epetra${static_lib_suffix} ; then
    AC_MSG_RESULT(yes)
    DEAL_II_TRILINOS_STATIC=yes
    DEAL_II_TRILINOS_LIBPREFIX="libtrilinos_"
  else
    AC_MSG_RESULT(no)
  fi

  dnl Make sure something is set at least
  if test "x${DEAL_II_TRILINOS_SHARED}${DEAL_II_TRILINOS_STATIC}" = "x" ; then
    AC_MSG_ERROR([Unable to determine whether Trilinos uses shared or static libraries.])
  fi


  dnl Now make sure the Trilinos libs are of the same kind as the ones we
  dnl produce here
  if test "x$enableshared" = "xyes" -a "x$DEAL_II_TRILINOS_SHARED" != "xyes" ; then
    AC_MSG_ERROR([When building deal.II with shared libraries, Trilinos also needs to be built with shared libraries])
  fi

  if test "x$enableshared" = "xno" -a "x$DEAL_II_TRILINOS_STATIC" != "xyes" ; then
    AC_MSG_ERROR([When building deal.II with shared libraries, Trilinos also needs to be built with shared libraries])
  fi

  dnl If we use shared libs (and we've made sure above that Trilinos provides
  dnl these as well), then set some of the LD_FLAGS and similar
  if test "x$enableshared" = "xyes" ; then
    LDFLAGS="$LDFLAGS -L$DEAL_II_TRILINOS_LIBDIR"
    if test "x$LD_PATH_OPTION" != "xno" ; then
      LDFLAGS="$LDFLAGS $LD_PATH_OPTION$DEAL_II_TRILINOS_LIBDIR"
    fi
  fi
])



dnl ------------------------------------------------------------
dnl Figure out which libraries Trilinos has built and that we
dnl need to link against. Also make sure we know their order
dnl
dnl Usage: DEAL_II_CHECK_TRILINOS_LIBS
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TRILINOS_LIBS, dnl
[
  AC_MSG_CHECKING(for the set of Trilinos libraries)

  dnl Trilinos' cmake invokation stores the set of libraries
  dnl in a special file for consumption of cmake at a later
  dnl time. We'll simply grep through it. Unfortunately, it has
  dnl changed place in Trilinos starting from version 10.8.0, and
  dnl at least 10.0.4 did not have the file at all. So test
  dnl whether the file is available at one location or another
  dnl and if it is at neither then fall back to a hardcoded
  dnl list of libraries
  if test -f $DEAL_II_TRILINOS_LIBDIR/cmake/Trilinos/TrilinosConfig.cmake ; then
    dnl This is the location for 10.8 and following
    DEAL_II_TRILINOS_LIBS="`grep Trilinos_LIBRARIES $DEAL_II_TRILINOS_LIBDIR/cmake/Trilinos/TrilinosConfig.cmake \
      | perl -pi -e 's/.*\"(.*)\".*/\1/g; s/;/ /g;'`"

    dnl Above, we may have set DEAL_II_TRILINOS_LIBPREFIX="trilinos_" to deal
    dnl with alternate naming conventions on Debian. However, if
    dnl we can get the list of libraries from Trilinos directly we
    dnl don't really need to do this any more
    DEAL_II_TRILINOS_LIBPREFIX=""
  else
    if test -f $DEAL_II_TRILINOS_INCDIR/TrilinosConfig.cmake ; then
      dnl The location prior to 10.8
      DEAL_II_TRILINOS_LIBS="`grep Trilinos_LIBRARIES $DEAL_II_TRILINOS_INCDIR/TrilinosConfig.cmake \
        | perl -pi -e 's/.*\"(.*)\".*/\1/g; s/;/ /g;'`"

      dnl Same as above
      DEAL_II_TRILINOS_LIBPREFIX=""
    else
      dnl Fall back to the fixed list. This should only be necessary
      dnl for Trilinos versions before 10.4. If, for a later version,
      dnl this happens, we don't want to fall back to this list, so
      dnl add an assertion.
      if test $DEAL_II_TRILINOS_VERSION_MAJOR = 10 -a \
              $DEAL_II_TRILINOS_VERSION_MINOR -lt 4 ; then
        :
      else
        AC_MSG_ERROR([package file TrilinosConfig.cmake not found])
      fi
      DEAL_II_TRILINOS_LIBS="stratimikosamesos stratimikosaztecoo stratimikosifpack stratimikosml stratimikos ml amesos belos ifpack aztecoo rtop sacado thyra thyraepetra thyraepetraext epetraext epetra teuchos triutils"
    fi
  fi
  AC_MSG_RESULT([$DEAL_II_TRILINOS_LIBS])
  AC_SUBST(DEAL_II_TRILINOS_LIBS)
])



dnl ------------------------------------------------------------
dnl Trilinos has headers that produce tons of warnings when
dnl used with -W -Wall (which includes -Wunused). Regrettable
dnl though it may be, these warnings pretty much drown everything
dnl else and we better disable some of the warnings to enable us
dnl to see through the clutter.
dnl
dnl Usage: DEAL_II_CHECK_TRILINOS_WARNINGS
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TRILINOS_WARNINGS, dnl
[
  OLD_CXXFLAGS="$CXXFLAGS"

  if test "$GXX" = yes ; then
    AC_MSG_CHECKING(whether we can use -Wno-unused to suppress Trilinos warnings)
    CXXFLAGS=-Wno-unused
    AC_TRY_COMPILE([], [;],
      [
        AC_MSG_RESULT(yes)
        CXXFLAGSG="$CXXFLAGSG -Wno-unused"
      ],
      [
        AC_MSG_RESULT(no)
      ])

    AC_MSG_CHECKING(whether we can use -Wno-overloaded-virtual to suppress Trilinos warnings)
    CXXFLAGS=-Wno-overloaded-virtual
    AC_TRY_COMPILE([], [;],
      [
        AC_MSG_RESULT(yes)
        CXXFLAGSG="$CXXFLAGSG -Wno-overloaded-virtual"
      ],
      [
        AC_MSG_RESULT(no)
      ])

    AC_MSG_CHECKING(whether we can use -Wno-extra to suppress Trilinos warnings)
    CXXFLAGS=-Wno-extra
    AC_TRY_COMPILE([], [;],
      [
        AC_MSG_RESULT(yes)
        CXXFLAGSG="$CXXFLAGSG -Wno-extra"
      ],
      [
        AC_MSG_RESULT(no)
      ])
  fi

  CXXFLAGS="${OLD_CXXFLAGS}"
])



dnl ------------------------------------------------------------
dnl Trilinos consists of a number of individual packages. We
dnl need several of those so we should make sure at configure
dnl that all necessary Trilinos packages were compiled and installed
dnl when Trilinos was built.
dnl
dnl Usage: DEAL_II_CHECK_TRILINOS_HEADER_FILES
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TRILINOS_HEADER_FILES, dnl
[
  OLD_CXXFLAGS="$CXXFLAGS"
  OLD_CPPFLAGS="$CPPFLAGS"

  CPPFLAGS="-I$DEAL_II_TRILINOS_INCDIR $CPPFLAGS"
  CXXFLAGS="-I$DEAL_II_TRILINOS_INCDIR $CXXFLAGS"

  dnl This is a subtle problem: Trilinos ML's ml_config.h has a
  dnl   #define HAVE_INTTYPES_H
  dnl without giving it a value. On the other hand, the previous
  dnl autoconf test for this header file will have put a
  dnl   #define HAVE_INTTYPES_H 1
  dnl into confdefs.h, which will lead to an error. Avoid this
  dnl problem by #undefining HAVE_INTTYPES_H for now and undoing
  dnl this later on again.
  dnl
  dnl Note that we have to do essentially the same trick as
  dnl well during compile time; see the block in AH_BOTTOM in
  dnl configure.in that goes into base/include/base/config.h
  cp confdefs.h confdefs.h.bak
  echo "#ifdef HAVE_INTTYPES_H" >> confdefs.h
  echo "#undef HAVE_INTTYPES_H" >> confdefs.h
  echo "#endif" >> confdefs.h

  dnl Sacado_cmath.hpp does things that aren't compatible
  dnl with the -std=c++0x flag of GCC, see deal.II FAQ.
  dnl Test whether that is indeed the case
  if test -f $DEAL_II_TRILINOS_INCDIR/Sacado_cmath.hpp ; then
    CXX_FLAGS_SAVED="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGSG -I$DEAL_II_TRILINOS_INCDIR"
    AC_MSG_CHECKING([whether Sacado_cmath.hpp is C++11 compatible])
    AC_TRY_COMPILE(
      [
        #include <Sacado_cmath.hpp>
      ],
      [;],
      [
        AC_MSG_RESULT([yes])
      ],
      [
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([*** Your Trilinos installation is not compatible with the C++ standard selected for this compiler. See the deal.II FAQ page for a solution. ***])
      ])
  else
    AC_MSG_ERROR([File $DEAL_II_TRILINOS_INCDIR/Sacado_cmath.hpp not found.])
  fi

  dnl Now just check that all headers we need are in fact there
  AC_CHECK_HEADERS([Amesos.h \
                    Epetra_CrsGraph.h \
                    Epetra_CrsMatrix.h \
                    Epetra_Import.h \
                    Epetra_LinearProblem.h \
                    Epetra_Map.h \
                    Epetra_MultiVector.h \
                    Epetra_Operator.h \
                    Epetra_SerialComm.h \
                    Epetra_Vector.h \
                    Ifpack.h \
                    ml_MultiLevelPreconditioner.h \
                    AztecOO.h \
                    AztecOO_Operator.h \
                    Sacado.hpp \
                    Teuchos_ParameterList.hpp \
                    Teuchos_RCP.hpp \
                    Teuchos_RefCountPtr.hpp
                   ],
                   [],
                   [
                     AC_MSG_ERROR([The Trilinos installation is missing one or more header files necessary for the deal.II Trilinos interfaces. Please re-install Trilinos with the missing Trilinos sub-packages enabled.])
                   ],
                   [])
  mv confdefs.h.bak confdefs.h

  CPPFLAGS="${OLD_CPPFLAGS}"
  CXXFLAGS="${OLD_CXXFLAGS}"
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

  dnl ------------------------------------------------------------
  dnl If MUMPS was requested and found, we had better check for
  dnl dependencies right here. First, SCALAPACK and then BLACS.
  if test "$USE_CONTRIB_MUMPS" = "yes" ; then
  dnl So here it goes...

    dnl ------------------------------------------------------------
    dnl Check whether SCALAPACK is installed
    AC_MSG_CHECKING([for SCALAPACK library directory])
    AC_ARG_WITH(scalapack,
      [AS_HELP_STRING([--with-scalapack=path/to/scalapack],
      [Specify the path to the scalapack installation; use this if you want to override the SCALAPACK_DIR environment variable.])],
      [dnl action-if-given (test)
       DEAL_II_SCALAPACK_DIR="$withval"
       AC_MSG_RESULT($DEAL_II_SCALAPACK_DIR)
       dnl Make sure that what was specified is actually correct
       if test ! -d $DEAL_II_SCALAPCK_DIR ; then
         AC_MSG_ERROR([The path to SCALAPACK specified with --with-scalapack does t point to a complete SCALAPACK installation])
       fi
      ],
      [dnl action-if-not-given (bail out)
       AC_MSG_ERROR([If MUMPS is used, the path to SCALAPACK must be specified with --with-scalapack])
      ])
    dnl ------------------------------------------------------------

    dnl ------------------------------------------------------------
    dnl Check whether BLACS is installed and BLAS architecture type
    AC_MSG_CHECKING([for BLACS library directory])
    AC_ARG_WITH(blacs,
      [AS_HELP_STRING([--with-blacs=path/to/blacs],
      [Specify the path to the BLACS installation; use this if you want to override the BLACS_DIR environment variable.])],
      [dnl action-if-given
       DEAL_II_BLACS_DIR="$withval"
       AC_MSG_RESULT($DEAL_II_BLACS_DIR)
       dnl Make sure that what was specified is actually correct
       if test ! -d $DEAL_II_BLACS_DIR     \
            -o ! -d $DEAL_II_BLACS_DIR/LIB ; then
       AC_MSG_ERROR([The path to BLACS specified with --with-blacs does not point to a complete BLACS installation])
       fi
      ],
      [dnl action-if-not-given (bail out)
       AC_MSG_ERROR([If MUMPS is used, the path to BLACS must be specified with --with-blacs])
      ])

    dnl BLACS labels libraries with "communications library",
    dnl "platform type" and "debug level" (see BLACS Bmake.inc for
    dnl details of the meaning of these terms). Finally, determine what
    dnl these are:
    AC_MSG_CHECKING([for BLACS library architecture])
    BLACS_COMM=`cat $DEAL_II_BLACS_DIR/Bmake.inc \
              | grep "COMMLIB = " \
              | perl -pi -e 's/.*LIB =\s+//g;'`
    BLACS_PLAT=`cat $DEAL_II_BLACS_DIR/Bmake.inc \
              | grep "PLAT = " \
              | perl -pi -e 's/.*PLAT =\s+//g;'`
    BLACS_DEBUG=`cat $DEAL_II_BLACS_DIR/Bmake.inc \
              | grep "BLACSDBGLVL = " \
              | perl -pi -e 's/.*DBGLVL =\s+//g;'`
    dnl and patch that together to make the BLACS architecture type:
    DEAL_II_BLACS_ARCH="$BLACS_COMM-$BLACS_PLAT-$BLACS_DEBUG"
    AC_MSG_RESULT($DEAL_II_BLACS_ARCH)
    dnl ------------------------------------------------------------

  fi
  dnl End check for MUMPS dependencies
  dnl ------------------------------------------------------------

  dnl If we do get this far then define a macro that says so:
  if test "$USE_CONTRIB_MUMPS" = "yes" ; then
    AC_DEFINE([DEAL_II_USE_MUMPS], [1],
              [Defined if a MUMPS installation was found and is going to be used              ])
    dnl and set an additional variable:
    DEAL_II_DEFINE_DEAL_II_USE_MUMPS=DEAL_II_USE_MUMPS
    dnl and finally set with_mumps if this hasn't happened yet:
    if test "x$with_mumps" = "x" ; then
      with_mumps="yes"
    fi
  fi
])

dnl ------------------------------------------------------------
dnl Check whether ARPACK is installed, and if so store the
dnl respective links
dnl
dnl Usage: DEAL_II_CONFIGURE_ARPACK
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_ARPACK, dnl
[
  AC_MSG_CHECKING([for ARPACK library directory])
  AC_ARG_WITH(arpack,
    [AS_HELP_STRING([--with-arpack=path/to/arpack],
    [Specify the path to the ARPACK installation, for which the include directory and lib directory are subdirs; use this if you want to override the ARPACK_DIR environment variable.])],
    [dnl action-if-given
     if test "x$withval" = "xno" ; then
       AC_MSG_RESULT([explicitly disabled])
       USE_CONTRIB_ARPACK=no
     else
       DEAL_II_ARPACK_DIR="$withval"
       AC_MSG_RESULT($DEAL_II_ARPACK_DIR)
       dnl Make sure that what was specified is actually correct
       if test ! -d $DEAL_II_ARPACK_DIR \
          ; then
         AC_MSG_ERROR([Path to ARPACK specified with --with-arpack does not point to a complete ARPACK installation])
       fi

       dnl ------------------------------------------------------------
       dnl Grab a meaningful name of the architeture ARPACK was
       dnl compiled for
       AC_MSG_CHECKING([for ARPACK library architecture])
       ARPACK_ARCH=`cat $DEAL_II_ARPACK_DIR/ARmake.inc \
                 | grep "PLAT = " \
                 | perl -pi -e 's/.*PLAT =\s+//g;'`
       DEAL_II_ARPACK_ARCH="$ARPACK_ARCH"
       AC_MSG_RESULT($DEAL_II_ARPACK_ARCH)

       USE_CONTRIB_ARPACK=yes
       AC_DEFINE(DEAL_II_USE_ARPACK, 1,
                 [Defined if an ARPACK installation was found and is
                  going to be used])
     fi
    ],
    [dnl action-if-not-given (do nothing)
     USE_CONTRIB_ARPACK=no
     AC_MSG_RESULT([no])
    ])
])


dnl ------------------------------------------------------------
dnl Check whether Metis is installed, and if so store the
dnl respective links
dnl
dnl Usage: DEAL_II_CONFIGURE_METIS
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_METIS, dnl
[
  dnl First check for the Metis directory

  AC_ARG_WITH(metis,
              [AS_HELP_STRING([--with-metis=/path/to/metis],
              [Specify the path to the Metis installation, of which the include and library directories are subdirs; use this if you want to override the METIS_DIR environment variable.])],
     [
        AC_MSG_CHECKING([for METIS library directory])
        USE_CONTRIB_METIS=yes
        DEAL_II_METIS_DIR="$withval"
        AC_MSG_RESULT($DEAL_II_METIS_DIR)

        dnl Make sure that what was specified is actually correct. The
        dnl libraries could be in either $DEAL_II_METIS_DIR/lib (metis was
        dnl make installed) or $DEAL_II_METIS_DIR/libmetis (metis was make
        dnl only and PETSc).
        if test ! -d $DEAL_II_METIS_DIR/lib \
           && test ! -d $DEAL_II_METIS_DIR/libmetis ; then
          AC_MSG_ERROR([Path to Metis specified with --with-metis does not point to a complete Metis installation])
        fi

        dnl If lib is not found, we must have libraries in libmetis
        dnl (which was found above).
        if test -d $DEAL_II_METIS_DIR/lib ; then
          DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR/lib"
        else
          DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR/libmetis"
        fi

        if test ! -d $DEAL_II_METIS_DIR/include ; then
          AC_MSG_ERROR([Path to Metis specified with --with-metis does not point to a complete Metis installation])
        fi

        DEAL_II_METIS_INCDIR="$DEAL_II_METIS_DIR/include"
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$METIS_DIR" != "x" ; then
          AC_MSG_CHECKING([for METIS from the environment])
          USE_CONTRIB_METIS=yes
          DEAL_II_METIS_DIR="$METIS_DIR"
          AC_MSG_RESULT($DEAL_II_METIS_DIR)

          dnl Make sure that what this is actually correct (see notes above).
          if test ! -d $DEAL_II_METIS_DIR/lib \
             && test ! -d $DEAL_II_METIS_DIR/libmetis ; then
            AC_MSG_ERROR([The path to Metis specified in the METIS_DIR environment variable does not point to a complete Metis installation])
          fi

          if test -d $DEAL_II_METIS_DIR/lib ; then
            DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR/lib"
          else
            DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR/libmetis"
          fi

          DEAL_II_METIS_INCDIR="$DEAL_II_METIS_DIR/include"

        else
          USE_CONTRIB_METIS=no
          DEAL_II_METIS_DIR=""
        fi
     ])

  if test "x$USE_CONTRIB_METIS" = "xyes" ; then
    AC_DEFINE(DEAL_II_USE_METIS, 1,
              [Defined if a Metis installation was found and is going
               to be used])
    LDFLAGS="$LDFLAGS -L$DEAL_II_METIS_LIBDIR -lmetis"

    if test "x$DEAL_II_LD_UNDERSTANDS_RPATH" = "xyes" ; then
      LDFLAGS="$LDFLAGS $LD_PATH_OPTION$DEAL_II_METIS_LIBDIR"
    fi

    dnl AC_MSG_CHECKING(for Metis version)
    dnl DEAL_II_METIS_VERSION=`cat $DEAL_II_METIS_DIR/VERSION`
    dnl AC_MSG_RESULT($DEAL_II_METIS_VERSION)
  fi
])


dnl ------------------------------------------------------------
dnl Check whether P4EST is to be used to parallelize meshes
dnl
dnl Usage: DEAL_II_CONFIGURE_P4EST
dnl
dnl ------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_P4EST, dnl
[
  AC_MSG_CHECKING(whether p4est will be used)

  AC_ARG_WITH(p4est,
             [AS_HELP_STRING([--with-p4est=/path/to/p4est],
             [Specify the path to the p4est installation; use this to distribute meshes on a cluster computer.])],
              use_p4est=$withval,
              use_p4est=no)

  if test "x$use_p4est" != "xno" ; then
    AC_MSG_RESULT(yes)

    dnl Verify that the p4est files are actually there
    if test ! -d "${use_p4est}/DEBUG" -o ! -d "${use_p4est}/FAST" ; then
      AC_MSG_ERROR([p4est directories $use_p4est/DEBUG or $use_p4est/FAST not found])
    fi

    dnl Make sure that we have also enabled MPI
    if test "x${DEAL_II_COMPILER_SUPPORTS_MPI}" != "x1" ; then
      AC_MSG_ERROR([When using p4est you also need to enable MPI.])
    fi

    dnl Right now, we always build p4est as shared lib, so make sure we
    dnl have built deal.II as a shared lib as well
    if test "x$enableshared" != "xyes" ; then
      AC_MSG_ERROR([When using p4est with shared libraries, you need to build
                    deal.II with shared libraries as well.])
    fi

    AC_DEFINE(DEAL_II_USE_P4EST, 1,
              [Defined if we are to use the p4est library to distribute
               meshes on a cluster computer.])
    USE_CONTRIB_P4EST=yes
    AC_SUBST(USE_CONTRIB_P4EST)

    DEAL_II_P4EST_DIR=${use_p4est}
    AC_SUBST(DEAL_II_P4EST_DIR)

    CXXFLAGSG="$CXXFLAGSG -I$use_p4est/DEBUG/include"
    CXXFLAGSO="$CXXFLAGSO -I$use_p4est/FAST/include"

    AC_MSG_CHECKING(for p4est library directory)
    if test -d "$use_p4est/DEBUG/lib64" -a -d "$use_p4est/FAST/lib64" ; then
      AC_MSG_RESULT(lib64)
      DEAL_II_P4EST_LIBDIR_NAME=lib64
    else
      if test -d "$use_p4est/DEBUG/lib" -a -d "$use_p4est/FAST/lib" ; then
        AC_MSG_RESULT(lib)
        DEAL_II_P4EST_LIBDIR_NAME=lib
      else
        AC_MSG_ERROR(could not determine whether p4est stores its library in lib/ or lib64/ directories)
      fi
    fi
    AC_SUBST(DEAL_II_P4EST_LIBDIR_NAME)

  else
    AC_MSG_RESULT(no)
  fi
])

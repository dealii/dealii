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

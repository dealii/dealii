
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


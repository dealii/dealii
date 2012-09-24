
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



    DEAL_II_EXPAND_PETSC_VECTOR="PETScWrappers::Vector"
    DEAL_II_EXPAND_PETSC_MPI_VECTOR="PETScWrappers::MPI::Vector"
    DEAL_II_EXPAND_PETSC_BLOCKVECTOR="PETScWrappers::BlockVector"
    DEAL_II_EXPAND_PETSC_MPI_BLOCKVECTOR="PETScWrappers::MPI::BlockVector"

    AC_DEFINE([DEAL_II_USE_PETSC_DEV], [1],
              [Defined if a PETSc installation was found and is not a release])



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


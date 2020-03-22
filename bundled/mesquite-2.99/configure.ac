dnl Process this file with autoconf 2.60 to produce a configure script
dnl for the CCA tutorial software installation.

AC_INIT(Mesquite, 2.99)
AC_CONFIG_AUX_DIR(./config_aux)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR(src/include/Mesquite.hpp)
AC_CANONICAL_SYSTEM
m4_ifdef([AM_SILENT_RULES],
[AM_INIT_AUTOMAKE([color-tests])],
[AM_INIT_AUTOMAKE])

if test "x$CXXFLAGS" = "x"; then
  have_user_CXXFLAGS=no
else
  have_user_CXXFLAGS=yes
fi

AC_DEFINE(MSQ_VERSION,["AC_PACKAGE_VERSION"],"Mesquite Version")
VERSION_MAJOR=`expr AC_PACKAGE_VERSION : '\([[0-9]]*\)'`
VERSION_MINOR=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.\([[0-9]]*\)'`
VERSION_PATCH=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.[[0-9]]*\.\(.*\)'`
test "x" != "x$VERSION_MAJOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
test "x" != "x$VERSION_MINOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
AC_DEFINE_UNQUOTED(MSQ_VERSION_MAJOR,$VERSION_MAJOR,"Mesquite Major Version")
AC_DEFINE_UNQUOTED(MSQ_VERSION_MINOR,$VERSION_MINOR,"Mesquite Minor Version")
if test "x" != "x$VERSION_PATCH"; then
  AC_DEFINE_UNQUOTED(MSQ_VERSION_PATCH,$VERSION_PATCH,"Mesquite Patch Level")
fi
AC_DEFINE(MSQ_VERSION_STRING,["AC_PACKAGE_STRING"],"Mesquite Version String")

AC_ARG_ENABLE( [namespace],
               [AC_HELP_STRING([--enable-namespace=NAME],[Specify alternate namespace])],
               [namespace=$enableval],
               [namespace=yes] )
case "x$namespace" in
  xyes)
    namespace="Mesquite$VERSION_MAJOR"
    ;;
  xno|x)
    AC_MSG_ERROR([Invalid value for --enable-namespace.  Mesquite cannot be configured with no namespace])
    ;;
esac
AC_DEFINE_UNQUOTED(MESQUITE_NS,[$namespace],"Mesquite namespace")
if test "x$namespace" != "xMesquite"; then
AC_DEFINE(MESQUITE_NS_ALIAS,[1],"Mesquite namespace alias")
fi


#------------------------------------------------------------------------------
# Standard stuff
#------------------------------------------------------------------------------
AC_PROG_MAKE_SET
AM_DISABLE_SHARED
LIBS="-lm"

#------------------------------------------------------------------------------
# Compile flag options -- need to do this *before* detecting the compiler
# otherwise cannot tell if user has CXXFLAGS defined already or if config
# set to a default value.
#------------------------------------------------------------------------------
AC_ARG_ENABLE( [release], 
  [AC_HELP_STRING([--enable-release],[(default), implies:
                                  --enable-compile-optimized
                                  --disable-debug-assertions])],
  [MSQ_DO_RELEASE="$enableval"], [MSQ_DO_RELEASE=""] )

AC_ARG_ENABLE( [debug], 
  [AC_HELP_STRING([--enable-debug],[Developer build, implies:
                                 --enable-debug-symbols
                                 --enable-debug-output=1,2
                                 --enable-trap-fpe
                                 --enable-silent-rules])], 
  [MSQ_DO_DEBUG="$enableval"], [MSQ_DO_DEBUG=""] )

m4_ifdef([AM_SILENT_RULES],[
  if test "xyes" = "x$MSQ_DO_DEBUG"; then
    AM_SILENT_RULES(yes)
  else
    AM_SILENT_RULES(no)
  fi
])

AC_ARG_ENABLE( [compile-optimized],
  [AC_HELP_STRING([--enable-compile-optimized],[Compile optimized (-O2)])],
  [MSQ_DO_OPTIMIZE="$enableval"], [MSQ_DO_OPTIMIZE=""] )

AC_ARG_ENABLE( [fast-native],
  [AC_HELP_STRING([--enable-fast-native],[Compile with fast math options, target build CPU, etc.])],
  [MSQ_DO_FAST="$enableval"], [MSQ_DO_FAST=""] )

AC_ARG_ENABLE( [debug-symbols], 
  [AC_HELP_STRING([--enable-debug-symbols],[Compile with debug symbols (-g)])],
  [MSQ_DEBUG_SYMBOLS="$enableval"], [MSQ_DEBUG_SYMBOLS=""] )

AC_ARG_ENABLE( [debug-assertions], 
  [AC_HELP_STRING([--disable-debug-assertions],[Compile without asserts (-DNDEBUG)])],
  [MSQ_DEBUG_ASSERTS="$enableval"], [MSQ_DEBUG_ASSERTS=""] )

AC_ARG_ENABLE( [debug-output], 
  [AC_HELP_STRING([--enable-debug-output(=n,m,...)],
                  [Enable debug output.  Optionally specify comma-separated 
                   list of flags to enable.  --enable-debug implies 
                   --enable-debug-output=1,2])],
  [MSQ_DEBUG_OUT="$enableval"], [MSQ_DEBUG_OUT=""] )

AC_ARG_ENABLE( [function-timers], 
  [AC_HELP_STRING([--enable-function-timers],[Enable function timers])], 
  [MSQ_DO_TIMERS="$enableval"], [MSQ_DO_TIMERS="no"] )

AC_ARG_ENABLE( [trap-fpe], 
  [AC_HELP_STRING([--enable-trap-fpe],
                  [Enable trapping of floating point exceptions.
                   Default is no unless --enable-debug.])], 
  [MSQ_TRAP_FPE="$enableval"], [MSQ_TRAP_FPE=""] )

AC_ARG_ENABLE( [quiet-make], 
  [AC_HELP_STRING([--enable-quiet-make],
                  [Much less verbose output during build.
                   Default is no unless --enable-debug.])], 
  [MSQ_QUITE_MAKE="$enableval"], [MSQ_QUITE_MAKE=""] )

AC_ARG_ENABLE( [32bit],
  [AC_HELP_STRING([--enable-32bit], [Force 32bit object code])],
  [MSQ_DO_32BIT=yes],[MSQ_DO_32BIT=no] )

AC_ARG_ENABLE( [64bit],
  [AC_HELP_STRING([--enable-64bit], [Force 64bit object code])],
  [MSQ_DO_64BIT=yes],[MSQ_DO_64BIT=no] )

# If neather debug or release is specified, enable release
if test -z "$MSQ_DO_RELEASE" -a -z "$MSQ_DO_DEBUG"; then
  if test "x$have_user_CXXFLAGS" = "xno"; then
    MSQ_DO_RELEASE=yes
  fi
fi

# if release, then enable appropriate sub-options
if test "$MSQ_DO_RELEASE" = yes; then
  if test -z "$MSQ_DO_OPTIMIZE"; then
    MSQ_DO_OPTIMIZE=yes
  fi
  if test -z "$MSQ_DEBUG_ASSERTS"; then
    MSQ_DEBUG_ASSERTS=no
  fi
fi

# if debug, then enable appropriate sub-options
if test "$MSQ_DO_DEBUG" = "yes"; then
  if test "x$enable_dependency_tracking" = "x"; then
    enable_dependency_tracking=yes
  fi
  if test -z "$MSQ_DEBUG_SYMBOLS"; then
    MSQ_DEBUG_SYMBOLS=yes
  fi
  if test -z "$MSQ_DEBUG_OUT"; then
    MSQ_DEBUG_OUT=yes
  fi
  if test -z "$MSQ_TRAP_FPE"; then
    MSQ_TRAP_FPE=yes
  fi
  if test -z "$MSQ_QUITE_MAKE"; then
    MSQ_QUITE_MAKE=yes
  fi
fi

CC_LIST="cc gcc cl egcs pgcc"
CXX_LIST="aCC cxx xlC_r xlC c++ g++ CC pgCC gpp cc++ cl FCC KCC RCC"

#------------------------------------------------------------------------------
# MPI part 1 (before compiler checks)
#------------------------------------------------------------------------------

AC_ARG_WITH( [mpi], [AC_HELP_STRING([--with-mpi(=DIR)],[Enable parallel support])],
             [WITH_MPI=$withval],[WITH_MPI=no] )

# Modify compiler list if MPI is requested
if test "xyes" = "x$WITH_MPI"; then  
  CC_LIST="mpicc $CC_LIST"
  CXX_LIST="mpicxx $CXX_LIST"
elif test "xno" != "x$WITH_MPI"; then
  if test "x" = "x$CC"; then
    CC="$WITH_MPI/bin/mpicc"
  fi
  if test "x" = "x$CXX"; then
    CXX="$WITH_MPI/bin/mpicxx"
  fi
fi


#------------------------------------------------------------------------------
# Construct compiler flags from options
#------------------------------------------------------------------------------

# Construct compiler options from above configure options
AC_PROG_CC( [$CC_LIST] )
AC_PROG_CXX( [$CXX_LIST] )
# get rid of annoying autoconf default CXXFLAGS
if test "x$have_user_CXXFLAGS" = "xno"; then
  CXXFLAGS=
fi
SNL_DETECT_CXX
SNL_REQUIRED_CXX_FLAGS
if test "$MSQ_DO_FAST" = "yes"; then
  SNL_CXX_COMPILE_FAST
elif test "$MSQ_DO_OPTIMIZE" = "yes"; then
  SNL_CXX_COMPILE_OPTIMIZED
fi
if test "$MSQ_DEBUG_SYMBOLS" = "yes"; then
  SNL_CXX_DEBUG_SYMBOLS
fi
if test "$MSQ_DEBUG_ASSERTS" = "no"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DNDEBUG"
fi
if test -n "$MSQ_DEBUG_OUT"; then
  if test "$MSQ_DEBUG_OUT" = "yes"; then 
    AM_CPPFLAGS="$AM_CPPFLAGS -DMSQ_ENABLE_DEBUG"
  elif test "$MSQ_DEBUG_OUT" != "no"; then
    AM_CPPFLAGS="$AM_CPPFLAGS -DMSQ_ENABLE_DEBUG=$MSQ_DEBUG_OUT"
  fi
fi
if test "$MSQ_TRAP_FPE" = "yes"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DMSQ_TRAP_FPE"
fi
if test "$MSQ_DO_TIMERS" = "yes"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DMSQ_USE_FUNCTION_TIMERS"
fi
if test "x$MSQ_QUITE_MAKE" = "xyes"; then
  MSQ_LIBTOOL_PREFIX='@echo "building $@...";'
  MSQ_LIBTOOL_FLAGS='--silent'
else
  MSQ_LIBTOOL_PREFIX=
  MSQ_LIBTOOL_FLAGS=
fi

if test "x$MSQ_DO_32BIT" = "xyes"; then
  SNL_CXX_COMPILE_32BIT
fi
if test "x$MSQ_DO_64BIT" = "xyes"; then
  if test "x$MSQ_DO_32BIT" = "xyes"; then
    AC_MSG_ERROR([--enable-32bit and --enable-64bit are mutually exclusive.])
  fi
  SNL_CXX_COMPILE_64BIT
fi
AC_LANG_CPLUSPLUS

#------------------------------------------------------------------------------
# Fortran stuff for ITAPS
#------------------------------------------------------------------------------
dnl AC_PROG_FC
dnl AC_ARG_VAR([FC], [FORTRAN compiler command])
dnl if test "x$FC" != "x"; then
dnl   AC_FC_WRAPPERS
dnl fi

#------------------------------------------------------------------------------
# LIBTOOL
#------------------------------------------------------------------------------
# Don't set up libtool until we're done messing with the compiler flags
# so we're sure libtool is set up correctly for the flags we want.
# Mostly need to make sure ar and nm are correct on IBM AIX, and that's
# handled inside our compiler detection code.

# We never use a fortran compiler, so skip testing for it
m4_undefine([AC_PROG_F77])
m4_defun([AC_PROG_F77],[])

# C++ compiler should already be set up, don't reconfigure
# with different settings
m4_undefine([AC_PROG_CXX])
m4_defun([AC_PROG_CXX],[])

AC_LIBTOOL_WIN32_DLL
AM_PROG_LIBTOOL


#------------------------------------------------------------------------------
# MPI part 2 (after compiler checks)
#------------------------------------------------------------------------------
if test "xno" != "x$WITH_MPI"; then
  AC_MSG_CHECKING([mpi.h])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_PREPROC_IFELSE([AC_LANG_PROGRAM( [[#include <stdio.h>
     #include <mpi.h>]],[])],
    [AC_MSG_RESULT([yes])],
    [save_CPPFLAGS="$CPPFLAGS"
     CPPFLAGS="$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"
     AC_PREPROC_IFELSE([AC_LANG_PROGRAM( [[#include <stdio.h>
         #include <mpi.h>]],[])],
         [AC_MSG_RESULT([yes (with -DMPICH_IGNORE_CXX_SEEK)]); AM_CPPFLAGS="$AM_CPPLFLAGS -DMPICH_IGNORE_CXX_SEEK" ],
         [AC_MSG_ERROR([mpi.h is not usable])]),
     CPPFLAGS="$save_CPPFLAGS"])
  AC_LANG_RESTORE
fi
AM_CONDITIONAL([WITH_MPI],[test "xno" != "x$WITH_MPI"])


#------------------------------------------------------------------------------
# Documentation
#------------------------------------------------------------------------------

# If old version of autoconf (automake?), these won't be defined for us
if test "x$htmldir" = "x"; then
  htmldir="$docdir"
  AC_SUBST(htmldir)
fi
if test "x$pdfdir" = "x"; then
  pdfdir="$docdir"
  AC_SUBST(pdfdir)
fi
if test "x$psdir" = "x"; then
  psdir="$docdir"
  AC_SUBST(psdir)
fi

AC_ARG_VAR([DOXYGEN],[command to generate Doxygen docs])
AC_PATH_PROG([DOXYGEN], [doxygen], [no])
AC_SUBST([DOXYGEN])

AC_ARG_VAR([LATEX],[command to generate latex documentation])
AC_PATH_PROG([LATEX], [latex], [no])  
AC_SUBST([LATEX])

AC_ARG_VAR([BIBTEX],[command to generate latex bibliography])
AC_PATH_PROG([BIBTEX], [bibtex], [no])  
AC_SUBST([BIBTEX])

AC_ARG_VAR([PS2PDF],[command to generate PDF from PostScript])
AC_PATH_PROG([PS2PDF], [ps2pdf], [no])  
AC_SUBST([PS2PDF])

AC_ARG_VAR([DVIPS],[command to generate PostScript from latex output])
AC_PATH_PROG([DVIPS], [dvips], [no])  
AC_SUBST([DVIPS])

AC_ARG_VAR([DOT],[Graphiz Utility (for Doxygen class graphs)])
AC_PATH_PROG([DOT], [dot], [no])  
HAVE_DOT=no
test "x$DOT" = "xno" || HAVE_DOT=yes
AC_SUBST([HAVE_DOT])
# Doxygen requires the directory, not the complete path
DOT_DIR=`expr "x$DOT" : 'x\(.*\)/[^/]*'`
AC_SUBST([DOT_DIR])

AC_ARG_ENABLE([api-doc],[AC_HELP_STRING([--enable-api-doc=FORMATLIST],
  [Generate Doxygen docs for API.
   Available formats are HTML,PDF,PS,ALL. Default is no docs.
   Multiple values must be separated by commas.])],
  [ARG_ENABLE_DOX=${enableval}],[ARG_ENABLE_DOX=no])
   
AC_ARG_ENABLE([user-guide],[AC_HELP_STRING([--enable-user-guide=FORMATLIST],
  [Generate Mesquite User's Guide
   Available formats are PDF,PS,ALL. Default is no docs.
   Multiple values must be separated by commas.])],
  [ARG_ENABLE_GUIDE=${enableval}],[ARG_ENABLE_GUIDE=no])

DOXYGEN_HTML=no
DOXYGEN_PDF=no
DOXYGEN_PS=no
USER_GUIDE_PDF=no
USER_GUIDE_PS=no
# No configuration options for these 
# - need to handle installation correctly before enabling
# - need to check for extra tools/support in installed doxygen before enabling?
DOXYGEN_HTMLHELP=no
DOXYGEN_CHI=no
DOXYGEN_RTF=no
DOXYGEN_MAN=no
DOXYGEN_XML=no
DOXYGEN_LATEX=no

# If user didn't explicitly request API dox
if test "x$ARG_ENABLE_DOX" = "x"; then
  # check which formats we can do
  if test "x$DOXYGEN" != "xno"; then
    if test "x$LATEX" != "xno"; then
      DOXYGEN_HTML=yes
    fi
  fi
# If user did reqest API docs explicitly
elif test "x$ARG_ENABLE_DOX" != "xno"; then
  # make sure we can do the requested formats
  if test "x$DOXYGEN" = "xno"; then
    AC_MSG_ERROR("doxygen not found.  Cannot generate API docs.")
  fi
  if test "x$LATEX" = "xno"; then
    AC_MSG_ERROR("latex not found.  Cannot generate API docs.")
  fi
  if test "x$ARG_ENABLE_DOX" == "xyes"; then
    if test "x$LATEX" != "xno"; then
      DOXYGEN_HTML=yes
    fi
  else
    list=`echo $ARG_ENABLE_DOX | tr :, ' '`
    for type in $list; do
      case $type in
        HTML|Html|html)
          DOXYGEN_HTML=yes
          ;;
        PDF|Pdf|pdf)
          DOXYGEN_PDF=yes
          DOXYGEN_LATEX=yes
          if test "x$DVIS" = "xno"; then AC_MSG_ERROR("Doxygen PDF requires dvips"); fi
          if test "x$PS2PDF" = "xno"; then AC_MSG_ERROR("Doxygen PDF requires ps2pdf"); fi
          ;;
        PS|Ps|ps)
          DOXYGEN_PS=yes
          DOXYGEN_LATEX=yes
          if test "x$DVIPS" = "xno"; then AC_MSG_ERROR("Doxygen PostScript requires dvips"); fi
          ;;
        ALL|All|all)
          DOXYGEN_HTML=yes
          if test "x$DVIPS" != "xno"; then 
            DOXYGEN_LATEX=yes
            DOXYGEN_PS=yes
            if test "x$PS2PDF" != "xno"; then DOXYGEN_PDF=yes; fi
          fi
          ;;
        *)
          AC_MSG_ERROR("Invalid API doc format: $type")
          ;;
      esac
    done
  fi
fi


# If user didn't explicitly request user's guide
if test "x$ARG_ENABLE_GUIDE" = "x"; then
  # check which formats we can do
  if test "x$LATEX" != "xno" -o "x$BIBTEX" != "xno" -o "x$DVIPS" != "xno"; then
    if test "x$PS2PDF" != "xno"; then
      USER_GUIDE_PDF=yes
    else
      USER_GUIDE_PS=yes
    fi
  fi
# If user did reqest user's guide explicitly
elif test "x$ARG_ENABLE_GUIDE" != "xno"; then
  # make sure we can do the requested formats
  if test "x$LATEX" = "xno"; then
    AC_MSG_ERROR("latex not found.  Cannot generate User's Guide.")
  fi
  if test "x$BIBTEX" = "xno"; then
    AC_MSG_ERROR("bibtex not found.  Cannot generate User's Guide.")
  fi
  if test "x$DVIPS" = "xno"; then
    AC_MSG_ERROR("dvips not found.  Cannot generate User's Guide.")
  fi
  if test "x$ARG_ENABLE_GUIDE" == "xyes"; then
    if test "x$PS2PDF" != "xno"; then
      USER_GUIDE_PDF=yes
    else
      USER_GUIDE_PS=yes
    fi
  else
    list=`echo $ARG_ENABLE_GUIDE | tr :, ' '`
    for type in $list; do
      case $type in
        PDF|Pdf|pdf)
          USER_GUIDE_PDF=yes
          if test "x$PS2PDF" = "xno"; then AC_MSG_ERROR("PDF User's Guide requires ps2pdf"); fi
          ;;
        PS|Ps|ps)
          USER_GUIDE_PS=yes
          ;;
        ALL|All|all)
          USER_GUIDE_PS=yes
          if test "x$PS2PDF" != "xno"; then USER_GUIDE_PDF=yes; fi
          ;;
        *)
          AC_MSG_ERROR("Invalid User's Guide format: $type")
          ;;
      esac
    done
  fi
fi

AM_CONDITIONAL([ENABLE_DOXYGEN],     [echo "$DOXYGEN_HTML $DOXYGEN_LATEX $DOXYGEN_HTMLHELP $DOXYGEN_CHI $DOXYGEN_RTF $DOXYGEN_MAN $DOXYGEN_XML" | grep yes > /dev/null] )
AM_CONDITIONAL([WITH_DOXYGEN_HTML],  [test "x$DOXYGEN_HTML"   = "xyes"])
AM_CONDITIONAL([WITH_DOXYGEN_PDF],   [test "x$DOXYGEN_PDF"    = "xyes"])
AM_CONDITIONAL([WITH_DOXYGEN_PS],    [test "x$DOXYGEN_PS"     = "xyes"])
AM_CONDITIONAL([WITH_USER_GUIDE_PDF],[test "x$USER_GUIDE_PDF" = "xyes"])
AM_CONDITIONAL([WITH_USER_GUIDE_PS], [test "x$USER_GUIDE_PS"  = "xyes"])
AC_SUBST(ENABLE_DOXYGEN)
AC_SUBST(DOXYGEN_HTML)
AC_SUBST(DOXYGEN_LATEX)
AC_SUBST(DOXYGEN_PDF)
AC_SUBST(DOXYGEN_PS)
AC_SUBST(DOXYGEN_HTMLHELP)
AC_SUBST(DOXYGEN_CHI)
AC_SUBST(DOXYGEN_RTF)
AC_SUBST(DOXYGEN_MAN)
AC_SUBST(DOXYGEN_XML)
AC_SUBST(USER_GUIDE_PDF)
AC_SUBST(USER_GUIDE_PS)

#-----------------------------------------------------------------------------
# Check for required headers
#-----------------------------------------------------------------------------
AC_HEADER_STDC
AC_CHECK_HEADERS(fcntl.h strings.h sys/file.h sys/time.h unistd.h)
AC_CHECK_HEADER(   ieeefp.h, AC_DEFINE(MSQ_HAVE_IEEEFP_H,  [1], [System has ieeefp.h]  ))
AC_CHECK_HEADER(   stdint.h, AC_DEFINE(MSQ_HAVE_STDINT_H,  [1], [System has stdint.h]  ))
AC_CHECK_HEADER( inttypes.h, AC_DEFINE(MSQ_HAVE_INTTYPES_H,[1], [System has inttypes.h]))
AC_CHECK_HEADER(  termios.h, [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_TERMIOS_H"])
AC_CHECK_HEADER( [sys/ioctl.h], [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_SYS_IOCTL_H"])

#-----------------------------------------------------------------------------
# Check for required functions
#-----------------------------------------------------------------------------
AC_CHECK_FUNCS( vsnprintf, [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_VSNPRINTF"] )
AC_CHECK_FUNCS( vsprintf, [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_VSPRINTF"] )
AC_CHECK_FUNCS( clock, [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_CLOCK"] )
AC_CHECK_FUNCS( times, [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_TIMES"] )
AC_CHECK_FUNCS( cbrt,  [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_CBRT"] )

# Use C++ compiler because C allows undefined functions, so these
# checks don't achive much if compiled as C.  Also, if C and C++
# compilers are mis-matched, what works for one may not work for the
# other and the C++ one is what is actually used for this in Mesquite.
AC_LANG_PUSH(C++)
AC_MSG_CHECKING( for fpsetmask );
AC_TRY_COMPILE( [#include <ieeefp.h>],
                [fpsetmask(FP_X_INV|FP_X_OFL|FP_X_DZ);],
                [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_FPSETMASK"
                 AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)] )
AC_MSG_CHECKING( for feenableexcept );
AC_TRY_COMPILE( [#define _GNU_SOURCE
                 #include <fenv.h>  ],
                [feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);],
                [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_FEENABLEEXCEPT"
                 AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)] )
AC_LANG_POP(C++)

#-----------------------------------------------------------------------------
# Check C++ environment
#-----------------------------------------------------------------------------
MSQ_CPLUSPLUS_FUNC

#-----------------------------------------------------------------------------
# Check for C++ features
#-----------------------------------------------------------------------------
AC_MSG_CHECKING( [for fd() in std::basic_file] )
AC_TRY_COMPILE([#include <fstream>],
               [using namespace std; ofstream f; int fd = f.rdbuf()->fd();],
               [CPPFLAGS="$CPPFLAGS -DFSTREAM_HAS_FD"; AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no)])

#------------------------------------------------------------------------------
# Other build options
#------------------------------------------------------------------------------

AC_ARG_WITH(exodus,
  [AC_HELP_STRING([--with-exodus(=DIR)],
                   [Enable exodusII support and specifiy directory])
AC_HELP_STRING([--without-exodus],[Disable exodusII support (default)])],
  [EXODUS_ARG=${withval}],[EXODUS_ARG=no])

AC_ARG_WITH(netcdf,
  [AC_HELP_STRING([--with-netcdf(=DIR)],
                  [ExodusII requires NetCDF - defaults to values for --with-exodus])
AC_HELP_STRING([--without-netcdf],[Skip NetCDF check])],
  [NETCDF_ARG=$withval], [NETCDF_ARG=])

AC_ARG_WITH(cppunit,
  [AC_HELP_STRING([--with-cppunit(=DIR)],[Specify directory where CppUnit is installed.])
AC_HELP_STRING([--without-cppunit],   [Disable CppUnit tests])],
  [CPPUNIT_ARG=${withval}], [CPPUNIT_ARG=])


#-------------------------------------------------------------------------------
# Configure different options
#-------------------------------------------------------------------------------

# Configure Exodus and NetCDF
AM_CONDITIONAL([WITH_EXODUS],[test "x$EXODUS_ARG" != "xno"])
if test "x$EXODUS_ARG" != "xno"; then
  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  AM_CPPFLAGS="$AM_CPPFLAGS -DMSQ_USING_EXODUS"
  
    # If user specified path for Exodus, add to appropriate vars
  if test "x$EXODUS_ARG" != "xyes"; then
    EXODUS_INC="-I${EXODUS_ARG}/include"
    EXODUS_LNK="-L${EXODUS_ARG}/lib"
    CPPFLAGS="$CPPFLAGS $EXODUS_INC"
    LDFLAGS="$LDFLAGS $EXODUS_LNK"
    AM_CPPFLAGS="$AM_CPPFLAGS $EXODUS_INC"
    AM_LDFLAGS="$AM_LDFLAGS $EXODUS_LNK"
  fi
  AM_LDFLAGS="$AM_LDFLAGS -lexoIIv2c"
    
    # Check for netcdf unless user explicitly said --without-netcdf
  if test "x$NETCDF_ARG" != "xno"; then
      # If user specified path for NetCDF, add to appropriate vars
    if test "x$NETCDF_ARG" != "xyes" -a "x$NETCDF_ARG" != "x" ; then
      NETCDF_INC="-I${NETCDF_ARG}/inc -I${NETCDF_ARG}/include"
      NETCDF_LNK="-L${NETCDF_ARG}/lib"
      CPPFLAGS="$CPPFLAGS $NETCDF_INC"
      LDFLAGS="$LDFLAGS $NETCDF_LNK"
      AM_CPPFLAGS="$AM_CPPFLAGS $NETCDF_INC"
      AM_LDFLAGS="$AM_LDFLAGS $NETCDF_LNK"
    fi
      # Check for NetCDF, but don't fail configure unless user
      # explicitly specified that netcdf should be included
    AC_CHECK_LIB( [netcdf], [nc_open], 
      [LIBS="$LIBS -lnetcdf"; AM_LDFLAGS="$AM_LDFLAGS -lnetcdf"],
      [if test "x$NETCDF_ARG" != "x"; then AC_MSG_ERROR([NetCDF library not found]); fi])
    AC_CHECK_HEADER( [netcdf.h], [], 
      [if test "x$NETCDF_ARG" != "x"; then AC_MSG_ERROR([netcdf.h not found]); fi ])
  fi
  
    # Check for ExodusII
  AC_CHECK_LIB( [exoIIv2c], [ex_close], [], [AC_MSG_ERROR([ExodusII library not found])] )
  AC_CHECK_HEADER( [exodusII.h], [], [AC_MSG_ERROR([exodusII.h not found])] )
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
fi



# CPPUnit
CPPUNIT_CPPFLAGS=
CPPUNIT_LDFLAGS=
if test "x$CPPUNIT_ARG" = "xno"; then
  HAVE_CPPUNIT=no
else
    # Incorporate user-specified search paths
  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  if test "x$CPPUNIT_ARG" != "x"; then
    if test "x$CPPUNIT_ARG" != "xyes"; then
      CPPUNIT_CPPFLAGS="-I$CPPUNIT_ARG/include"
      CPPUNIT_LDFLAGS="-L$CPPUNIT_ARG/lib"
      CPPFLAGS="$CPPFLAGS $CPPUNIT_CPPFLAGS"
      LDFLAGS="$CPPFLAGS $CPPUNIT_LDFLAGS"
    fi
  fi
    # Set some variables 
  HAVE_CPPUNIT=yes
  CPPUNIT_LDFLAGS="$CPPUNIT_LDFLAGS -lcppunit"
    # Check that CppUnit exists
  AC_CHECK_LIB( [dl], [dlopen], 
    [LDFLAGS="$LDFLAGS -ldl"; CPPUNIT_LDFLAGS="$CPPUNIT_LDFLAGS -ldl"] )
  AC_LANG_PUSH(C++)
  AC_CHECK_HEADER( [cppunit/Test.h], [], [HAVE_CPPUNIT=no] )
  AC_CHECK_LIB( [cppunit], [main], [], [HAVE_CPPUNIT=no] )
  AC_LANG_POP(C++)
    # If user explicitly specified --with-cppunit, fail if it was not found
  if test "x$CPPUNIT_ARG" != "x"; then
    if test "x$HAVE_CPPUNIT" == "xno"; then
      AC_MSG_ERROR([CppUnit not found])
    fi
  fi
    # restore some modified state
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
fi
AM_CONDITIONAL([WITH_CPPUNIT],[test "x$HAVE_CPPUNIT" = "xyes"])

#------------------------------------------------------------------------------
# ITAPS
#------------------------------------------------------------------------------

ITAPS_API([imesh], [IMESH], [iMesh])
ITAPS_API([igeom], [IGEOM], [iGeom])
ITAPS_API([irel],  [IREL],  [iRel],  [IMESH], [IGEOM])
ITAPS_API([imeshp],[IMESHP],[iMeshP],[IMESH])

if test "x$ENABLE_ITAPS" = "xyes"; then
  if test "x$IBASE_INCL" = "x"; then
    AC_MSG_ERROR([Cannot enable support for any ITAPS API w/out iBase.h])
  fi
fi
AM_CPPFLAGS="${AM_CPPFLAGS}"
AM_CONDITIONAL([ENABLE_ITAPS],[test "x$ENABLE_ITAPS" = "xyes"])
AC_SUBST(IBASE_INCL)

#------------------------------------------------------------------------------
# The End
#------------------------------------------------------------------------------
if test "x" = "x$docdir"; then
  docdir='${datadir}/doc/mesquite'
  AC_SUBST(docdir)
fi

MESQUITE_LIBS="$MESQUITE_LIBS -lm"
AC_SUBST(SED)
AC_SUBST(DOXYGEN)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)
AC_SUBST(MESQUITE_LIBS)
AC_SUBST(MSQ_LIBTOOL_PREFIX)
AC_SUBST(MSQ_LIBTOOL_FLAGS)
AC_SUBST(CPPUNIT_CPPFLAGS)
AC_SUBST(CPPUNIT_LDFLAGS)
AC_SUBST(MESQUITE_IMESH_MAKE_INCLUDE)

# We don't want to automatically generate mesquite_config.h.in
# automatically using autoheader.  It puts #defines in the file
# that are not appropriate for a library (will conflict with app's
# #defines).  Try to disable it
#AUTOHEADER="touch \$@"
AUTOHEADER=":"
AC_SUBST(AUTOHEADER)

TEST_MAKE_INCLUDE="include ../../src/cppflags.make"
UTIL_MAKE_INCLUDE="include ../src/cppflags.make"
AC_SUBST(TEST_MAKE_INCLUDE)
AC_SUBST(UTIL_MAKE_INCLUDE)

# The first header is the header that autoheader creates.
# It should not be used because it contains stuff that is
# not restricted to the Mesquite namespace and therefore
# may conflict with the config.h of the application using
# mesquite.
AC_CONFIG_HEADERS([src/include/junk_config.h])
AC_CONFIG_HEADERS([src/include/mesquite_config.h
                   src/include/mesquite_version.h])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 itaps/Makefile
                 meshFiles/Makefile
                 testSuite/Makefile
                 testSuite/meshfiles.h
                 testSuite/ActiveSetTest/Makefile
                 testSuite/algorithm_test/Makefile
                 testSuite/analytical_grad_3D/Makefile
                 testSuite/convert/Makefile
                 testSuite/benchmark_tests/Makefile
                 testSuite/paraboloid_domain_test/Makefile
                 testSuite/headers/Makefile
                 testSuite/idft_time/Makefile
                 testSuite/igeom/Makefile
                 testSuite/imesh/Makefile
                 testSuite/jacobi/Makefile
                 testSuite/laplacian_test/Makefile
                 testSuite/laplacian_polygon_test/Makefile
                 testSuite/nongradient_test/Makefile
                 testSuite/pyramid/Makefile
                 testSuite/simple_hybrid_test/Makefile
                 testSuite/test_1/Makefile
                 testSuite/transform/Makefile
                 testSuite/tutorial/Makefile
                 testSuite/tutorial/tutorial.make
                 testSuite/unit/Makefile
                 testSuite/untangle_test/Makefile
                 testSuite/parallel_untangle_shape/Makefile 
                 testSuite/parallel_smooth_laplace/Makefile
                 testSuite/wedge/Makefile
                 testSuite/wrapper_tests/Makefile
                 testSuite/2d_target/Makefile
                 testSuite/2d_metrics/Makefile
                 testSuite/2d_formulation/Makefile
                 testSuite/synchronous/Makefile
                 testSuite/high_aspect_ratio/Makefile
                 testSuite/higher_order/Makefile
                 testSuite/slaved/Makefile
                 utils/Makefile
                 doc/Makefile
                 doc/user/Makefile
                 doc/user/config.tex
                 doc/user/doxygen/Makefile
                 doc/user/doxygen/Mesquite-user.dox
                 ])
                 
# Remove any source depenency data so that it is regenerated.
# This isn't necessary for a correctly working configure script,
# but from a support POV, it is a easier to explain and remember
# that the configure script should be re-run without having to
# also remember to do "rm .deps/*" when a source file is removed
# from the build.
AC_CONFIG_COMMANDS_POST([rm -f .deps/*])                 
                 
# This header is the header that autoheader creates.
# It should not be used because it contains stuff that is
# not restricted to the Mesquite namespace and therefore
# may conflict with the config.h of the application using
# mesquite.
AC_CONFIG_COMMANDS([junk_config.h],
                   [echo '#error "Dont use this header"' > src/include/junk_config.h])

# Generate output files
AC_OUTPUT

echo "CXXFLAGS = $CXXFLAGS"

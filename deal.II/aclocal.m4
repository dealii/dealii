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



dnl -------------------------------------------------------------
dnl Determine the C++ compiler in use. Return the name and possibly
dnl version of this compiler in GXX_VERSION.
dnl
dnl Usage: DEAL_II_DETERMINE_CXX_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_DETERMINE_CXX_BRAND, dnl
[
  if test "$GXX" = yes ; then
    dnl find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    case "$GXX_VERSION_STRING" in
      *"egcs-1.1"*)
  	AC_MSG_RESULT(C++ compiler is egcs-1.1)
  	GXX_VERSION=egcs1.1
  	;;
      *2.95*)
  	AC_MSG_RESULT(C++ compiler is gcc-2.95)
  	GXX_VERSION=gcc2.95
  	;;
      *2.96*)
  	AC_MSG_RESULT(C++ compiler is gcc-2.96)
  	GXX_VERSION=gcc2.96
  	;;
      *2.97*)
  	AC_MSG_RESULT(C++ compiler is gcc-2.97)
  	GXX_VERSION=gcc2.97
  	;;
      *3.0*)
  	AC_MSG_RESULT(C++ compiler is gcc-3.0)
  	GXX_VERSION=gcc3.0
  	;;
      *3.1*)
  	AC_MSG_RESULT(C++ compiler is gcc-3.1)
  	GXX_VERSION=gcc3.1
  	;;
      *3.2*)
  	AC_MSG_RESULT(C++ compiler is gcc-3.2)
  	GXX_VERSION=gcc3.2
  	;;
      *2.4* | *2.5* | *2.6* | *2.7* | *2.8*)
  	dnl These compilers are too old to support a useful subset
  	dnl of modern C++, so we don't support them
  	AC_MSG_RESULT(C++ compiler is $GXX_VERSION_STRING)
  	AC_MSG_ERROR(C++ compiler is not supported)
  	;;
      *)
  	AC_MSG_RESULT(C++ compiler is unknown but accepted gcc version)
  	GXX_VERSION=gcc-other
  	;;
    esac
  else
    dnl Check other (non-gcc) compilers
  
    dnl Check for IBM xlC. For some reasons, depending on some environment
    dnl variables, moon position, and other reasons unknown to me, the
    dnl compiler displays different names in the first line of output, so
    dnl check various possibilities
    is_ibm_xlc="`($CXX 2>&1) | egrep 'VisualAge C++|C Set ++|C for AIX Compiler'`"
    if test "x$is_ibm_xlc" != "x"  ; then
      dnl Ah, this is IBM's C++ compiler. Unfortunately, we don't presently
      dnl know how to check the version number, so assume that is sufficiently
      dnl high...
      AC_MSG_RESULT(C++ compiler is IBM xlC)
      GXX_VERSION=ibm_xlc
    else
  
      dnl Check whether we are dealing with the MIPSpro C++ compiler
      is_mips_pro="`($CXX -version 2>&1) | grep MIPSpro`"
      if test "x$is_mips_pro" != "x" ; then
        AC_MSG_RESULT(C++ compiler is MIPSpro C++ compiler)
        GXX_VERSION=MIPSpro
      else
  
        dnl Intel's ICC C++ compiler?
        is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R) C++ Compiler'`"
        if test "x$is_intel_icc" != "x" ; then
          AC_MSG_RESULT(C++ compiler is Intel ICC)
          GXX_VERSION=intel_icc
        else
  
          dnl Or DEC's cxx compiler?
          is_dec_cxx="`($CXX -V 2>&1) | grep 'Compaq C++'`"
          if test "x$is_dec_cxx" != "x" ; then
            AC_MSG_RESULT(C++ compiler is Compaq cxx)
            GXX_VERSION=compaq_cxx
          else
  
  	    dnl Sun Workshop?
            is_sun_cc="`($CXX -V 2>&1) | grep 'Sun WorkShop'`"
            if test "x$is_sun_cc" != "x" ; then
              AC_MSG_RESULT(C++ compiler is Sun Workshop compiler)
              GXX_VERSION=sun_workshop
            else
  
  	      dnl Sun Forte?
              is_sun_forte_cc="`($CXX -V 2>&1) | grep 'Forte'`"
              if test "x$is_sun_forte_cc" != "x" ; then
                AC_MSG_RESULT(C++ compiler is Sun Forte compiler)
                GXX_VERSION=sun_forte
              else
  
  	      dnl KAI C++?
  	      is_kai_cc="`($CXX -V 2>&1) | grep 'KAI C++'`"
  	      if test "x$is_kai_cc" != "x" ; then
  	        AC_MSG_RESULT(compile is KAI C++)
  	        GXX_VERSSION=kai_cc
  	      else
  
                  dnl  Aw, nothing suitable found...
                  AC_MSG_ERROR(Unrecognized compiler, sorry)
                  exit 1
                fi
              fi
  	    fi
          fi
        fi
      fi
    fi
  fi
])





dnl -------------------------------------------------------------
dnl Set C++ compiler flags to their default values. They will be 
dnl modified according to other options in later steps of
dnl configuration
dnl
dnl CXXFLAGSO  : flags for optimized mode
dnl CXXFLAGSG  : flags for debug mode
dnl CXXFLAGSPIC: flags for generation of object files that are suitable
dnl              for shared libs
dnl LDFLAGSPIC : flags needed for linking of object files to shared
dnl              libraries
dnl
dnl Usage: DEAL_II_SET_CXX_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_SET_CXX_FLAGS, dnl
[
  dnl First the flags for gcc compilers
  if test "$GXX" = yes ; then
    CXXFLAGSO="$CXXFLAGS -O2 -Wuninitialized -felide-constructors -ftemplate-depth-32"
    CXXFLAGSG="$CXXFLAGS -DDEBUG -ansi -pedantic -Wall -W -Wpointer-arith -Wwrite-strings -Wmissing-prototypes -Winline -Woverloaded-virtual -Wstrict-prototypes -Wsynth -Wsign-compare -Wconversion -Wswitch -ftemplate-depth-32"
    CXXFLAGSPIC="-fPIC"  
    LDFLAGSPIC="-fPIC"

    dnl set some flags that are specific to some versions of the
    dnl compiler:
    dnl - egcs1.1 yielded incorrect code with some loop unrolling
    dnl - after egcs1.1, the optimization flag -fstrict-aliasing was
    dnl   introduced, which enables better optimizations for
    dnl   well-written C++ code. we believe that deal.II falls into that
    dnl   category and thus enable the flag 
    dnl - egcs1.1 yielded incorrect code with vtable-thunks. thus disable
    dnl   them for egcs1.1. however, if on Linux, disabling them
    dnl   prevents programs from being linked, so take the risk of broken
    dnl   thunks on this platform
  
    case "$GXX_VERSION" in
      egcs1.1)
          case "$target" in
            *linux*)
                ;;
  
            *)
                CXXFLAGSG = "$CXXFLAGSG -fno-vtable-thunks"
                CXXFLAGSO = "$CXXFLAGSO -fno-vtable-thunks"
                ;;
          esac
          ;;
  
      dnl All other gcc versions
      *)
          CXXFLAGSO="$CXXFLAGSO -funroll-loops -funroll-all-loops -fstrict-aliasing"
          ;;
    esac
  
    dnl - after gcc2.95, some flags were deemed obsolete for C++
    dnl   (and are only supported for C any more), so only define them for
    dnl   previous compilers
  
    case "$GXX_VERSION" in
      egcs1.1 | gcc2.95)
          CXXFLAGSG="$CXXFLAGSG -Wmissing-declarations -Wbad-function-cast -Wtraditional -Wnested-externs"
          CXXFLAGSO="$CXXFLAGSO -fnonnull-objects"
          ;;
  
      *)
          ;;
    esac
  
    dnl Some gcc compiler versions have a problem when using an unsigned count
    dnl in the std::advance function. Unfortunately, this also happens 
    dnl occasionally from within the standard library, so we can't prevent the
    dnl warning messages. Since this is annoying, switch of the flag -W which
    dnl causes this.
    DEAL_II_CHECK_ADVANCE_WARNING
    if test "x$DEAL_II_ADVANCE_WARNING" = "xyes" ; then
      CXXFLAGSG="`echo $CXXFLAGSG | perl -pi -e 's/-W //g;'`"
    fi
  
  else
    dnl Non-gcc compilers
  
    case "$GXX_VERSION" in
      ibm_xlc)
          CXXFLAGSG="$CXXFLAGS -DDEBUG -check=bounds -info=all -qrtti=all"
          CXXFLAGSO="$CXXFLAGS -O2 -w -qansialias -qrtti=all"
          CXXFLAGSPIC="-fPIC"
          LDFLAGSPIC="-fPIC"
          ;;
  
      MIPSpro)
          CXXFLAGSG="$CXXFLAGS -DDEBUG -LANG:std"
          CXXFLAGSO="$CXXFLAGS -LANG:alias_const=ON -LANG:std -w"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-KPIC"
          ;;
  
      intel_icc)
          dnl Disable some compiler warnings, as they often are wrong on
          dnl our code:
          dnl #175: `subscript out of range' (doesn't take into account that
          dnl       some code is only reachable for some dimensions)
          dnl #327: `NULL reference is not allowed' (this happens when we
          dnl       write "*static_cast<double*>(0)" or some such thing,
          dnl       which we do to create invalid references)
          dnl #525: `type "DataOutBase::DataOutBase" is an inaccessible type
          dnl       (allowed for compatibility)' (I don't understand what the
          dnl       compiler means)
          CXXFLAGSG="$CXXFLAGS -Kc++eh -Krtti -w1 -wd175 -wd525 -wd327 -DDEBUG -inline_debug_info"
          CXXFLAGSO="$CXXFLAGS -Kc++eh -Krtti -O2 -tpp6 -axiMK -unroll -w0"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-KPIC"
          ;;
  
      compaq_cxx)
          dnl Disable some warning messages:
          dnl #175: `subscript out of range' (detected when instantiating a
          dnl       template and looking at the indices of an array of
          dnl       template dependent size, this error is triggered in a
          dnl       branch that is not taken for the present space dimension)
          dnl #236 and
          dnl #237: `controlling expression is constant' (in while(true), or
          dnl       switch(dim))
          dnl #381: `extra ";" ignored' (at function or namespace closing
          dnl       brace)
          dnl #487: `Inline function ... cannot be explicitly instantiated'
          dnl       (also reported when we instantiate the entire class)
          dnl #1136:`conversion to integral type of smaller size could lose data'
          dnl       (occurs rather often in addition of int and x.size(),
          dnl       because the latter is size_t=long unsigned int on Alpha)
          dnl #1156:`meaningless qualifiers not compatible with "..." and "..."'
          dnl       (cause unknown, happens when taking the address of a
          dnl       template member function)
          dnl #111 and
          dnl #1182:`statement either is unreachable or causes unreachable code'
          dnl       (happens in switch(dim) clauses for other dimensions than
          dnl       the present one)
          dnl
          dnl Also disable the following error:
          dnl #265: `class "..." is inaccessible' (happens when we try to
          dnl       initialize a static member variable in terms of another
          dnl       static member variable of the same class if the latter is
          dnl       not public and therefore not accessible at global scope in
          dnl       general. I nevertheless think that this is valid.)
          dnl
          dnl Besides this, choose the most standard conforming mode of the
          dnl compiler, i.e. -model ansi and -std strict_ansi. Unfortunately,
          dnl we have to also add the flag -implicit_local (generating implicit
          dnl instantiations of template with the `weak' link flag) since
          dnl otherwise not all templates are instantiated (also some from the
          dnl standards library are missing).
  
          CXXFLAGSG="$CXXFLAGS -model ansi -std strict_ansi -w1 -msg_display_number -timplicit_local -DDEBUG"
          CXXFLAGSO="$CXXFLAGS -model ansi -std strict_ansi -w2 -msg_display_number -timplicit_local -fast"
  
          for i in 175 236 237 381 487 1136 1156 111 1182 265 ; do
            CXXFLAGSG="$CXXFLAGSG -msg_disable $i"
            CXXFLAGSO="$CXXFLAGSO -msg_disable $i"
          done
  
          dnl If we use -model ansi to compile the files, we also have to
          dnl specify it for linking
          LDFLAGS="$LDFLAGS -model ansi"
  
          dnl For some reason, cxx also forgets to add the math lib to the
          dnl linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"
  
          dnl if necessary: -shared, -pthread 
          dnl Should one use -compress? -distinguish_nested_enums?
          dnl                -nousing_std? -pch? -noimplicit_include?

          CXXFLAGSPIC="-shared"
          LDFLAGSPIC="-shared"
          ;;
  
      sun_workshop | sun_forte)
          CXXFLAGSG="$CXXFLAGS -DDEBUG"
          CXXFLAGSO="$CXXFLAGS -fast"
          CXXFLAGSPIC="-KPIC"
          ;;
  
      *)
          AC_MSG_ERROR(No compiler options for this C++ compiler
                       specified at present)
          ;;
    esac
  fi
])



dnl -------------------------------------------------------------
dnl Depending on what compiler we use, set the flag necessary to
dnl obtain debug output. We will usually use -ggdb, but on DEC Alpha 
dnl OSF1, this leads to stabs symbols that are too long for the system 
dnl assembler. We will therefore check whether the assembler can handle
dnl these symbols by a rather perverse function with many templates, 
dnl and if the assembler can't handle them, then use -gstabs 
dnl instead. This reduces debugging possibilities, but no other
dnl way is known at present.
dnl
dnl For all compilers other than gcc, use -g instead and don't check.
dnl 
dnl Usage: DEAL_II_SET_CXX_DEBUG_FLAG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_SET_CXX_DEBUG_FLAG, dnl
[
  if test "$GXX" = yes ; then
    AC_MSG_CHECKING(whether -ggdb works for long symbols)
    case "$target" in
       dnl On Alpha, use the special treatment
       alpha*-osf*)
            CXXFLAGS="-ggdb $CXXFLAGSG"
            AC_TRY_COMPILE(
              [
#               include <string>
#               include <map>
                using namespace std;

                typedef map<string,map<string,pair<string,string> > > T;

                bool f(T& t1, const T* t2) {
                  t1["s"] = map<string,pair<string,string> >();
                  map<string,map<string,pair<string,string> > >
                        ::const_iterator i2=t1.begin();
                  map<string,map<string,pair<string,string> > >
                        ::const_iterator i1=t2->begin();
                  return (i1==i2);
                }
              ],
            [
                  ;
            ],
            [
                CXXFLAGSG="-ggdb $CXXFLAGSG"
              AC_MSG_RESULT(yes)
            ],
            [
                CXXFLAGSG="-gstabs $CXXFLAGSG"
              AC_MSG_RESULT(no -- using -gstabs)
            ])
          ;;
  
       dnl For all other systems assume that -ggdb works (we can't make the
       dnl test  above the default, as stabs are not the default debugging
       dnl format on many systems, and we only want to use it where necessary
       *)
          AC_MSG_RESULT(yes)
            CXXFLAGSG="-ggdb $CXXFLAGSG"
          ;;
    esac
  
  else
    dnl Non-gcc compilers use -g instead of -ggdb
    CXXFLAGSG="-g $CXXFLAGSG"
  fi
])



dnl -------------------------------------------------------------
dnl Determine the F77 compiler in use. Return the name and possibly
dnl version of this compiler in F77_VERSION.
dnl
dnl Usage: DEAL_II_DETERMINE_F77_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_DETERMINE_F77_BRAND, dnl
[
  if test "x$F77" != "x" ; then
    F77_VERSION_STRING="`($F77 -v 2>&1)`"
    if test -n "`echo $F77_VERSION_STRING | grep \"GNU F77\"`" ; then
      dnl Yes, this is a GNU g77 version. find out the right version
      G77_VERSION_STRING="`($F77 -v 2>&1) | grep \"gcc version\"`"
      case "$G77_VERSION_STRING" in
        *"egcs-1.1"*)
            AC_MSG_RESULT(F77 compiler is egcs-1.1)
            F77_VERSION=egcs1.1
            ;;
        *2.95*)
            AC_MSG_RESULT(F77 compiler is gcc-2.95)
  	  F77_VERSION=gcc2.95
  	  ;;
        *2.96*)
  	  AC_MSG_RESULT(F77 compiler is gcc-2.96)
  	  F77_VERSION=gcc2.96
  	  ;;
        *2.97*)
  	  AC_MSG_RESULT(F77 compiler is gcc-2.97)
  	  F77_VERSION=gcc2.97
  	  ;;
        *3.0*)
  	  AC_MSG_RESULT(F77 compiler is gcc-3.0)
  	  F77_VERSION=gcc3.0
  	  ;;
        *3.1*)
  	  AC_MSG_RESULT(F77 compiler is gcc-3.1)
  	  F77_VERSION=gcc3.1
  	  ;;
        *3.2*)
  	  AC_MSG_RESULT(F77 compiler is gcc-3.2)
  	  F77_VERSION=gcc3.2
  	  ;;
        *2.4* | *2.5* | *2.6* | *2.7* | *2.8*)
  	  dnl These compilers are too old to support a useful subset
  	  dnl of modern C++, so we don't support them
  	  AC_MSG_RESULT(F77 compiler is $G77_VERSION_STRING)
  	  AC_MSG_ERROR(F77 compiler is not supported)
  	  ;;
        *)
  	  AC_MSG_RESULT(F77 compiler is unknown but accepted gcc version)
  	  F77_VERSION=gcc-other
  	  ;;
      esac
  
    else
  
      dnl No GNU g77 version, something else. Try to find out what it is:
      if test -n "`echo $F77_VERSION_STRING | grep \"XL Fortran for AIX\"`" ; then
  
        dnl This is the "XL Fortran for AIX" compiler
        F77_VERSION=AIXF77
        AC_MSG_RESULT(F77 compiler is AIX Fortran77)
  
      else
  
        dnl Umh, still something else unknown. Try to find it out by a
        dnl different method (-V instead of -v):
        F77_VERSION_STRING=`($F77 -V 2>&1)`
        if test -n "`echo $F77_VERSION_STRING | grep \"WorkShop Compilers\"`" ; then
  
          dnl OK, this is the Sun Fortran77 compiler
  	  AC_MSG_RESULT(F77 compiler is Sun WorkShop f77)
  	  F77_VERSION="SunF77"
  
        else
  
          dnl If we can detect IRIX's f77 somehow, then the following flags
          dnl might be appropriate:
          F77_VERSION_STRING=`($F77 -version 2>&1)`
  	  if test -n "`echo $F77_VERSION_STRING | grep MIPSpro`" ; then
  	    AC_MSG_RESULT(F77 compiler is MIPSpro f77)
  	    F77_VERSION="MIPSproF77"

          else
  
           dnl Now, this is a hard case, we have no more clues...
            F77_VERSION=
  	    AC_MSG_RESULT(F77 compiler is unkown. no flags set!)
          fi
        fi
      fi
    fi
  fi
])



dnl -------------------------------------------------------------
dnl Set F77 compiler flags to their default values.
dnl
dnl F77FLAGSO  : flags for optimized mode
dnl F77FLAGSG  : flags for debug mode
dnl F77LIBS    : libraries to link with when using F77
dnl F77FLAGSPIC: since we don't presently use libtool to link F77 libs,
dnl              this is the flag to use for f77 compilers when using
dnl              shared libraries
dnl
dnl Usage: DEAL_II_SET_F77_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_SET_F77_FLAGS, dnl
[
  case "$F77_VERSION" in
    egcs-1.1 | gcc2.95 | gcc2.96 | gcc2.97 | gcc3.0 | gcc3.1 | gcc3.2)
        F77FLAGSG="$FFLAGS -ggdb -DDEBUG -pedantic -W -Wall"
        F77FLAGSO="$FFLAGS -O2"
  
        dnl Some flags can only be set for some compilers as others either
        dnl did not accept them or were buggy on them (see the explanation
        dnl for CXXFLAGS for an explanation of some of these cases)
        if test "x$F77_VERSION" != "xegcs1.1" ; then
          F77FLAGSO="$F77FLAGSO -funroll-loops -funroll-all-loops -fstrict-aliasing"
        fi
  
	F77FLAGSPIC="-fPIC"
        F77LIBS="$F77LIBS -lg2c"
  
	;;

    AIXF77)
        F77FLAGSG="$FFLAGS -g"
        F77FLAGSO="$FFLAGS -O3 -w"
        F77LIBS="$F77LIBS -lxlf90"  

        F77FLAGSPIC="unknown!"
        ;;
  
    SunF77)
  	F77FLAGSG="$FFLAGS -silent -g"
  	F77FLAGSO="$FFLAGS -silent -O3 -w"
  	F77LIBS="$F77LIBS -lF77 -lsunmath -lM77"
        F77FLAGSPIC="-PIC"
	;;

    MIPSproF77)
  	F77FLAGSG="$FFLAGS -ansi -g"
  	F77FLAGSO="$FFLAGS -O3 -woffall"
  	F77LIBS="$F77LIBS -lftn"
  
  	F77FLAGSPIC="shared -KPIC"
	;;

    *)
	AC_MSG_ERROR(No compiler options for F77 compiler
                     "$F77_VERSION" specified at present)
        ;;
  esac
])



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




dnl -------------------------------------------------------------
dnl Versions of gcc on different platforms use a multitude of flags to
dnl denote thread safe libraries and the like. They are, for example
dnl -threads/-pthread/-mthread, sometimes *thread, sometimes *threads. 
dnl Find out which is the right one on the present platform
dnl
dnl Usage: DEAL_II_FIND_THREAD_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_GET_THREAD_FLAGS, dnl
[
  AC_MSG_CHECKING(for platform specific thread flags)
  AC_LANG(C++)
  for i in threads mt pthread pthreads mthreads Kthread kthread invalid_last_entry; do
    CXXFLAGS="$CXXFLAGSG -$i"
    DEAL_II_TRY_COMPILER_FLAG(
	[
	  thread_flag=$i
	  CXXFLAGSG="$CXXFLAGSG -$i"
	  CXXFLAGSO="$CXXFLAGSO -$i"
          LDFLAGS="$LDFLAGS -$i"
	
	  dnl The right parameter was found, so exit
	  break
   	])
  done
  if test "$thread_flag" = invalid_last_entry ; then
	AC_MSG_RESULT(no flag found!)
	AC_MSG_ERROR(Could not determine multithreading flag for this platform. Aborting!)
  fi
  AC_MSG_RESULT(-$thread_flag)
])



dnl -------------------------------------------------------------
dnl Test whether multithreading support is requested. This
dnl does not tell deal.II to actually use it, but the
dnl compiler flags can be set to allow for it. If the user specified
dnl --enable-multithreading, then set $enablemultithreading=yes,
dnl otherwise $enablemultithreading=no.
dnl
dnl Usage:
dnl   DEAL_II_CHECK_MULTITHREADING
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MULTITHREADING, dnl
[
  AC_ARG_ENABLE(multithreading,
  [  --enable-multithreading Set compiler flags to allow for
                             multithreaded programs],
    enablemultithreading=$enableval,
    enablemultithreading=no)
])



dnl -------------------------------------------------------------
dnl If multithreading support is requested, figure out the right
dnl compiler flags to use:
dnl - Use the right -threads/-pthread/-mthread option
dnl - Set preprocessor directives if necessary
dnl - __USE_MALLOC tells gcc to use thread safe STL allocators
dnl - _REENTRANT is a flag that is used in the standard UNIX
dnl   headers to make reentrant functions (with suffix _r) declared
dnl
dnl Do not use AC_DEFINE when using these two flags, since that would
dnl put them into config.h, instead of the compiler flags. Then, however,
dnl it would be necessary to include config.h as _first_ file in all
dnl files, since the flags change the things we include. It is therefore
dnl better to put the flags into the command line, since then we have them
dnl defined for all include files.
dnl
dnl Usage:
dnl   DEAL_II_SET_MULTITHREADING_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_SET_MULTITHREADING_FLAGS, dnl
[
  if test "$enablemultithreading" = yes ; then
    if test "$GXX" = yes ; then
      DEAL_II_GET_THREAD_FLAGS
      DEAL_II_THREAD_CPPFLAGS
  
      CXXFLAGSG="$CXXFLAGSG -D__USE_MALLOC -D_REENTRANT"
      CXXFLAGSO="$CXXFLAGSO -D__USE_MALLOC -D_REENTRANT"
    else
      if test "x$GXX_VERSION" = "xibm_xlc" ; then
        CXXFLAGSG = "$CXXFLAGSG -threaded"  
        CXXFLAGSO = "$CXXFLAGSO -threaded"
      else
        dnl Other compiler
        AC_MSG_ERROR(No threading compiler options for this C++ compiler
                     specified at present)
        exit 1
      fi
    fi
  fi
])



dnl -------------------------------------------------------------
dnl Test whether MT code shall be used through ACE. If so, set
dnl $withmultithreading to the path to ACE
dnl
dnl Usage:
dnl   DEAL_II_CHECK_USE_ACE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_USE_ACE, dnl
[
  AC_ARG_WITH(multithreading,
  [  --with-multithreading=DIR use DIR as path to the ACE library],
      withmultithreading=$withval,
      withmultithreading=no)
  if test "$withmultithreading" != no ; then
    AC_MSG_CHECKING(for ACE)
    if test -d "$withmultithreading" ; then
      AC_MSG_RESULT(found)
    else
      AC_MSG_RESULT(not found)
      AC_MSG_ERROR(Invalid ACE path)
    fi
  
    AC_DEFINE(DEAL_II_USE_MT, 1, 
              [Flag indicating whether the library shall be compiler for
               multithreaded applications])
  fi
])



dnl -------------------------------------------------------------
dnl If the MT code shall be used through ACE, we might have to 
dnl adjust the compiler flags: Possibly remove -ansi -pedantic from
dnl compiler flags again, since ACE yields hundreds of error messages
dnl with these flags.
dnl
dnl Usage:
dnl   DEAL_II_SET_ACE_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_SET_ACE_FLAGS, dnl
[
  if test "$enablemultithreading" = yes ; then
    if test "$GXX" = yes ; then
      DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS
      if test "x$deal_II_ace_remove_ansi" = "xyes" ; then
        CXXFLAGSG="`echo $CXXFLAGSG | perl -pi -e 's/-ansi//g;'`"
        CXXFLAGSO="`echo $CXXFLAGSO | perl -pi -e 's/-ansi//g;'`"
      fi
      if test "x$deal_II_ace_remove_pedantic" = "xyes" ; then
        CXXFLAGSG="`echo $CXXFLAGSG | perl -pi -e 's/-pedantic//g;'`"
        CXXFLAGSO="`echo $CXXFLAGSO | perl -pi -e 's/-pedantic//g;'`"
      fi
    fi
  fi
])




dnl -------------------------------------------------------------
dnl Check whether some backward compatibility features are disabled
dnl Usage:
dnl   DEAL_II_CHECK_COMPAT_BLOCKER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_COMPAT_BLOCKER, dnl
[
  AC_ARG_ENABLE(compat-blocker,
  [  --enable-compat-blocker=mapping Block functions that implicitely
                                     assume a Q1 mapping],
      enable_compat_blocker=$enableval,
      enable_compat_blocker="")

  dnl Replace the comma-separated list by a space-separated one
  disable_compat=`echo $enable_compat_blocker | perl -pi -e 's/,/ /g;'`

  dnl Check that each entry is an allowed one
  for i in $disable_compat ; do 
    case $i in
      mapping)
        AC_MSG_RESULT(Disabling backward compatibility feature: "$i")
        ;;
      *) 
        AC_MSG_ERROR(Backward compatibility feature "$i" unknown)
        ;;
    esac
  done

  dnl Now for each known feature, either disable it or enable it.
  dnl Default is to enable. In order to have these flags in config.h,
  dnl it is necessary to AC_DEFINE them actually by name, rather than
  dnl by some loop variable, since otherwise autoheader can't generate
  dnl an entry for config.h for this variable
  for i in mapping ; do
    uppercase=`echo $i | perl -pi -e 'tr/a-z/A-Z/;'`
    if test -n "`echo $disable_compat | grep $i`"  ; then
      compat_value=false
    else
      compat_value=true
    fi
  
    case $i in
      mapping)
        AC_DEFINE_UNQUOTED(DEAL_II_COMPAT_MAPPING,$compat_value,
                           [Backward compatibility support for functions
                              and classes that do not take an explicit
                              mapping variable, but rather use a default
                              Q1 mapping instead])
  	;;

      *)
  	AC_MSG_ERROR(Backward compatibility feature "$i" unknown)
  	;;
    esac
  done
])



dnl -------------------------------------------------------------
dnl Try whether the set of compiler flags in CXXFLAGS is reasonable, i.e.
dnl does not result in compiler messages (which are then produced by
dnl unknown or unrecognized compiler flags. This macro is mostly copied
dnl from the definition of AC_TRY_COMPILE, but the first two arguments are
dnl omitted and some other things are also simplified.
dnl
dnl Usage:
dnl   DEAL_II_TRY_COMPILER_FLAG([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_TRY_COMPILER_FLAG, dnl
[
  echo 'int main() { return 0; }' > conftest.$ac_ext

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



dnl -------------------------------------------------------------
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
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_GETRUSAGE, dnl
[
  AC_MSG_CHECKING(whether getrusage is properly declared)
  AC_LANG(C++)
  AC_TRY_COMPILE(
    [
#include <sys/resource.h>
    ],
    [
      rusage *ru;
      getrusage(RUSAGE_SELF,ru);
    ],
    [
	AC_MSG_RESULT(yes)
    ],
    [
	AC_MSG_RESULT(no)
	CXXFLAGSG="$CXXFLAGSG -DNO_HAVE_GETRUSAGE"
	CXXFLAGSO="$CXXFLAGSO -DNO_HAVE_GETRUSAGE"
    ])
]
)      




dnl -------------------------------------------------------------
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
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ISNAN, dnl
[
  DEAL_II_CHECK_ISNAN_FLAG(debug, $CXXFLAGSG,
			   CXXFLAGSG="$deal_II_isnan_flag $CXXFLAGSG")
  DEAL_II_CHECK_ISNAN_FLAG(optimized, $CXXFLAGSO,
		 	   CXXFLAGSO="$deal_II_isnan_flag $CXXFLAGSO")
]
)



dnl -------------------------------------------------------------
dnl The following function actually performs the check for the right flag
dnl for `isnan'. If a flag is found, the third argument is executed, and 
dnl the right flag is available in $deal_II_isnan_flag.
dnl
dnl Usage: DEAL_II_CHECK_ISNAN_FLAG("description of options set",
dnl                                 "compiler options set",
dnl                                 action when flag is found)
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ISNAN_FLAG, dnl
[
  AC_MSG_CHECKING(whether isnan is declared with $1 flags)
  AC_LANG(C++)
  CXXFLAGS=$2
  deal_II_isnan_flag=""
  AC_TRY_COMPILE(
    [
#include <cmath>
    ],
    [
	double d=0;
	isnan (d);
    ],
    [
	AC_MSG_RESULT(yes)
	deal_II_isnan_flag="-DHAVE_ISNAN"
	$3
    ])


  if test "x$deal_II_isnan_flag" = "x" ; then
    dnl Simply using isnan doesn't work. On Microsoft Windows systems, the
    dnl function is called _isnan, so check that
    AC_TRY_COMPILE(
      [
#include <cmath>
      ],
      [
	  double d=0;
	  _isnan (d);
      ],
      [
	  AC_MSG_RESULT(yes)
	  deal_II_isnan_flag="-DHAVE_UNDERSCORE_ISNAN"
	  $3
      ])
  fi


  dnl Let's see whether we _now_ have found something
  if test "x$deal_II_isnan_flag" = "x" ; then
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
	  dnl We found a flag by which isnan is defined; store
	  dnl this flag and exit the loop
	  deal_II_isnan_flag="-DHAVE_ISNAN $testflag"
	  break;
	])

      dnl If that didn't work (and it didn't as we are still inside the
      dnl loop), then try the _isnan function (maybe we are on a
      dnl Microsoft Windows system)
      AC_TRY_COMPILE(
        [
#include <cmath>
	],
	[
	  double d=0;
	  _isnan (d);
	],
	[
	  dnl We found a flag by which isnan is defined; store
	  dnl this flag and exit the loop
	  deal_II_isnan_flag="-DHAVE_UNDERSCORE_ISNAN $testflag"
	  break;
	])
    done

    dnl if no such flag was found, then abort ./configure since
    dnl the library will not be compilable on this platform
    dnl without knowledge of the right flag
    if test "x$deal_II_isnan_flag" = "x" ; then
      AC_MSG_RESULT(no.)
    else
      dnl we found something, lets us it
      AC_MSG_RESULT(using $testflag)
      $3
    fi
  fi
])      



dnl -------------------------------------------------------------
dnl rand_r is defined for some compiler flag combinations, but not for
dnl others. check that. note that since these are C++ flags, we can't
dnl just use AC_CHECK_FUNCS
dnl
dnl Usage: DEAL_II_CHECK_RAND_R
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_RAND_R, dnl
[
  AC_LANG(C++)
  CXXFLAGS=$CXXFLAGSG
  AC_MSG_CHECKING(for rand_r)
  AC_TRY_COMPILE(
	[
#include <cstdlib>
	],
	[
int seed = 0;
int i=rand_r(&i);
	],
	[
	  AC_MSG_RESULT(found)
	  AC_DEFINE(HAVE_RAND_R, 1, 
                    [Define if you have the rand_r function])
	],
	[
	  AC_MSG_RESULT(no)
	]
  )
])



dnl -------------------------------------------------------------
dnl On some systems (well, DEC Alphas are the only ones we know of),
dnl gcc2.95 throws the hands in the air if it sees one of the AssertThrow
dnl calls, and dies with an internal compiler error. If this is the case,
dnl we disable AssertThrow and simply replace it with an `abort' if the
dnl condition is not satisfied.
dnl
dnl Usage: DEAL_II_CHECK_ASSERT_THROW("description of options set",
dnl                                   "compiler options set",
dnl                                   action if compiler crashes)
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ASSERT_THROW, dnl
[
  AC_MSG_CHECKING(whether AssertThrow works with $1 flags)
  AC_LANG(C++)
  CXXFLAGS=$2
  AC_TRY_COMPILE(
    [
#include <exception>
#include <iostream>
#include <cstdlib>

#ifndef __GNUC__
#  define __PRETTY_FUNCTION__ "(unknown)"
#endif
class ExceptionBase : public std::exception {
  public:
    ExceptionBase ();
    ExceptionBase (const char* f, const int l, const char *func,
		   const char* c, const char *e);
    virtual ~ExceptionBase () throw();
    void SetFields (const char *f, const int   l, const char *func,
		    const char *c, const char *e);
    void PrintExcData (std::ostream &out) const;
    virtual void PrintInfo (std::ostream &out) const;
    virtual const char * what () const throw ();
  protected:
    const char  *file;
    unsigned int line;
    const char  *function, *cond, *exc;
};

template <class exc>
void __IssueError_Assert (const char *file,
			  int         line,
			  const char *function,
			  const char *cond,
			  const char *exc_name,
			  exc         e){
  e.SetFields (file, line, function, cond, exc_name);
  std::cerr << "--------------------------------------------------------"
	    << std::endl;
  e.PrintExcData (std::cerr);
  e.PrintInfo (std::cerr);
  std::cerr << "--------------------------------------------------------"
	    << std::endl;  
  std::abort ();
};

template <class exc>
void __IssueError_Throw (const char *file,
			 int         line,
			 const char *function,
			 const char *cond,
			 const char *exc_name,
			 exc         e) {
				   // Fill the fields of the exception object
  e.SetFields (file, line, function, cond, exc_name);
  throw e;
};

#define AssertThrow(cond, exc)                                    \
  {                                                               \
    if (!(cond))                                                  \
      __IssueError_Throw (__FILE__,                               \
			  __LINE__,                               \
			  __PRETTY_FUNCTION__, #cond, #exc, exc); \
  }

#define DeclException0(Exception0)  \
class Exception0 :  public ExceptionBase {}

namespace StandardExceptions 
{
  DeclException0 (ExcInternalError);
};
using namespace StandardExceptions;
    ],
    [
	AssertThrow (false, ExcInternalError());
    ],
    [
	AC_MSG_RESULT(yes)
    ],
    [
	AC_MSG_RESULT(no)
	$3
    ])
])



dnl -------------------------------------------------------------
dnl IBM xlC 5.0 from the VisualAge C++ pack has a bug with the following 
dnl code. We can work around it if we insert code like "using namespace std;"
dnl in the right place, but we'd like to do so only if absolutely necessary.
dnl Check whether the compiler which we are using has this bug.
dnl
dnl Usage: DEAL_II_CHECK_IBM_XLC_ERROR
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_IBM_XLC_ERROR, dnl
[
  AC_MSG_CHECKING(for std::vector bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      namespace std {
        template <class _Ty>                             class allocator {};
        template<class _Ty, class _Ax = allocator<_Ty> > class vector{};
      }

      struct X {};
      template <int dim> void g (const std::vector<X> &x);

      void f ()  {
        std::vector<X> x;
        g<1> (x);
      };
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. trying to work around)
      AC_DEFINE(XLC_WORK_AROUND_STD_BUG, 1, 
                [Define if we have to work around a bug in IBM's xlC compiler])
    ])
])



dnl -------------------------------------------------------------
dnl gcc2.95 doesn't have the std::iterator class, but the standard
dnl requires it, so check whether we have to work around it
dnl
dnl Usage: DEAL_II_HAVE_STD_ITERATOR
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_STD_ITERATOR, dnl
[
  AC_MSG_CHECKING(for std::iterator class)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <iterator>
      class MyIterator : public std::iterator<std::bidirectional_iterator_tag,int>
      {};
    ],
    [],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_STD_ITERATOR_CLASS, 1,
	        [Define if the compiler's library in use provides a
                 std::iterator class (early gcc versions did not)])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])




dnl -------------------------------------------------------------
dnl Up to early gcc2.95 releases, the i/ostringstream classes were not
dnl available. check their availability, or whether we have to fall back
dnl to the old strstream classes.
dnl
dnl Usage: DEAL_II_HAVE_STD_STRINGSTREAM
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_STD_STRINGSTREAM, dnl
[
  AC_MSG_CHECKING(for std::i/ostringstream classes)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <sstream>
    ],
    [
	std::istringstream i;
	std::ostringstream o;
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_STD_STRINGSTREAM, 1, 
                [Define if the compiler's library in use provides
                 std::i/ostringstream classes (early gcc versions did not)])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Check for existence of the __builtin_expect facility of newer
dnl gcc compilers. This can be used to hint the compiler's branch
dnl prediction unit in some cases. We use it in the AssertThrow
dnl macros.
dnl
dnl Usage: DEAL_II_HAVE_BUILTIN_EXPECT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_BUILTIN_EXPECT, dnl
[
  AC_MSG_CHECKING(for __builtin_expect)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	bool f();
    ],
    [
	if (__builtin_expect(f(),false));
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_BUILTIN_EXPECT, 1, 
                [Define if the compiler provides __builtin_expect])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl For gcc 2.95, when using the random_shuffle function with -ansi
dnl compiler flag, there is an error saying that lrand48 is
dnl undeclared. Fix that by declaring that function ourselves if
dnl necessary.
dnl
dnl Usage: DEAL_II_HAVE_LRAND48_DECLARED
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_LRAND48_DECLARED, dnl
[
  AC_MSG_CHECKING(whether lrand48 needs to be declared with -ansi)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <vector>
#include <algorithm>

void f()
{
  std::vector<unsigned int> new_indices;
  std::random_shuffle (new_indices.begin(), new_indices.end());
};
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_DECLARE_LRAND48, 1, 
                [Define if you have the rand_r function])
    ])
])



dnl -------------------------------------------------------------
dnl When compiling with ACE thread support, there are many constructs
dnl that are not allowed in C++, or that yield warnings when compiling with
dnl -ansi -pedantic. Check this, and if that is the case, set the variables
dnl $deal_II_ace_remove_ansi and $deal_II_ace_remove_pedantic to "yes".
dnl
dnl Usage: DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS, dnl
[
  AC_MSG_CHECKING(whether compilation with ACE disallows flags)
  AC_LANG(C++)
  CXXFLAGS="-ansi -I$withmultithreading"
  AC_TRY_COMPILE(
    [
#  include <ace/Thread_Manager.h>
#  include <ace/Synch.h>
    ],
    [],
    [
      deal_II_ace_remove_ansi="no"
    ],
    [
      deal_II_ace_remove_ansi="yes"
    ])
  CXXFLAGS="-pedantic -Werror -I$withmultithreading"
  AC_TRY_COMPILE(
    [
#  include <ace/Thread_Manager.h>
#  include <ace/Synch.h>
    ],
    [],
    [
      deal_II_ace_remove_pedantic="no"
    ],
    [
      deal_II_ace_remove_pedantic="yes"
    ])

  if test $deal_II_ace_remove_ansi = "no" ; then
    if test $deal_II_ace_remove_pedantic = "no" ; then 
      AC_MSG_RESULT(no)
    else
      AC_MSG_RESULT(-pedantic)
    fi
  else
    if test $deal_II_ace_remove_pedantic = "no" ; then 
      AC_MSG_RESULT(-ansi)
    else
      AC_MSG_RESULT(-ansi -pedantic)
    fi
  fi
])



dnl -------------------------------------------------------------
dnl gcc versions up to 2.95.3 had a problem with the std::advance function,
dnl when the number of steps forward was given by an unsigned number, since
dnl a comparison >=0 was performed on this number which leads to a warning
dnl that this comparison is always true. This is, at the best, annoying
dnl since it crops up at several places where std::advance is called from
dnl inside the library. Check whether the present version of the compiler
dnl has this problem
dnl
dnl This test is only called if gcc is used.
dnl
dnl Usage: DEAL_II_CHECK_ADVANCE_WARNING
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ADVANCE_WARNING, dnl
[
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -Werror"
  AC_MSG_CHECKING(for std::advance warning)
  AC_TRY_COMPILE(
	[
#include <map>
#include <vector>
	],
	[
  std::map<unsigned int, double> m;
  std::vector<std::pair<unsigned int, double> > v;
  v.insert (v.end(), m.begin(), m.end());
	],
	[
	  dnl compilation succeeded, no warning occured for the above code
	  AC_MSG_RESULT(no)
          DEAL_II_ADVANCE_WARNING=no
	],
	[
	  AC_MSG_RESULT(yes)
          DEAL_II_ADVANCE_WARNING=yes
	]
  )
])



dnl -------------------------------------------------------------
dnl Check whether the numeric_limits classes are available
dnl
dnl Usage: DEAL_II_HAVE_STD_NUMERIC_LIMITS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_STD_NUMERIC_LIMITS, dnl
[
  AC_MSG_CHECKING(for std::numeric_limits classes)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <limits>
    ],
    [
	unsigned int i = std::numeric_limits<unsigned int>::min();
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_STD_NUMERIC_LIMITS, 1, 
                [Define if the compiler's library in use provides
                 std::numeric_limits classes in the appropriate header file])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Check whether the numeric_limits classes are available
dnl
dnl Usage: DEAL_II_HAVE_STD_OSTREAM_HEADER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_STD_OSTREAM_HEADER, dnl
[
  AC_MSG_CHECKING(for <ostream> header)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <ostream>
void f (const std::ostream &out);
    ],
    [
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_STD_OSTREAM_HEADER, 1, 
                [Define if the compiler provides an <ostream> header file])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Check whether the numeric_limits classes are available
dnl
dnl Usage: DEAL_II_HAVE_STD_OSTREAM_HEADER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_STD_IOSFWD_HEADER, dnl
[
  AC_MSG_CHECKING(for <iosfwd> header)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <iosfwd>
void f (const std::ostream &out);
    ],
    [
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_STD_IOSFWD_HEADER, 1, 
                [Define if the compiler provides an <iosfwd> header file])
    ],
    [
      AC_MSG_RESULT(no)
    ])
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
  else
    AC_MSG_RESULT(none found)
    USE_CONTRIB_HSL=no
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
  AC_CHECK_FILE($TECHOME/lib/tecio.a,
		TECPLOT_LIBRARY_PATH=$TECHOME/lib/tecio.a)
  AC_CHECK_FILE($TEC80HOME/lib/tecio.a,
		TECPLOT_LIBRARY_PATH=$TEC80HOME/lib/tecio.a)
  AC_CHECK_FILE($TEC90HOME/lib/tecio.a,
		TECPLOT_LIBRARY_PATH=$TEC90HOME/lib/tecio.a)
  AC_CHECK_FILE($TECHOME/include/TECIO.h,
		TECPLOT_INCLUDE_PATH=$TECHOME/include)
  AC_CHECK_FILE($TEC80HOME/include/TECIO.h,
	        TECPLOT_INCLUDE_PATH=$TEC80HOME/include)
  AC_CHECK_FILE($TEC90HOME/include/TECIO.h,
	        TECPLOT_INCLUDE_PATH=$TEC90HOME/include)

  if (test -r $TECPLOT_LIBRARY_PATH && \
      test -r $TECPLOT_INCLUDE_PATH/TECIO.h) ; then
    AC_DEFINE(DEAL_II_HAVE_TECPLOT, 1,
	      [Flag indicating whether the library shall be compiled to use the Tecplot interface])
  fi
])



dnl -------------------------------------------------------------
dnl Check whether CXXFLAGSG and CXXFLAGSO are a valid combination
dnl of flags or if there are contradicting flags in them.
dnl
dnl Usage: DEAL_II_CHECK_CXXFLAGS_CONSISTENCY
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CXXFLAGS_CONSISTENCY, dnl
[
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_MSG_CHECKING(for consistency of CXXFLAGSG flags)
  AC_TRY_COMPILE(
    [],
    [;],
    [
      AC_MSG_RESULT(yes)
    ],
    [
      AC_MSG_ERROR(invalid combination of flags!)
      exit 1;
    ])
  
  CXXFLAGS="$CXXFLAGSO"
  AC_MSG_CHECKING(for consistency of CXXFLAGSO flags)
  AC_TRY_COMPILE(
    [],
    [;],
    [
      AC_MSG_RESULT(yes)
    ],
    [
      AC_MSG_ERROR(invalid combination of flags!)
      exit 1;
    ])
])



dnl -------------------------------------------------------------
dnl Check whether F77FLAGS are a valid combination of flags or if
dnl there are contradicting flags in them.
dnl
dnl Usage: DEAL_II_CHECK_F77FLAGS_CONSISTENCY
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_F77FLAGS_CONSISTENCY, dnl
[
dnl     Err, well -- we'd like to have checked these flags, but autoconf 2.13
dnl     has a problem here: when writing the compile file, a newline is
dnl     missing, leading to an error when executing ./configure. So:
dnl     while this is not fixed in autoconf, disable the respective tests
dnl
dnl     Postface: with autoconf 2.50, AC_LANG_PROGRAM(Fortran 77) is called
dnl     internally, but this unfortunately results in a warning from
dnl     autoconf. So still do not enable the code
dnl  AC_LANG(Fortran 77)
dnl  FFLAGS="$F77FLAGSG"
dnl  AC_MSG_CHECKING(for consistency of F77FLAGSG flags)
dnl  AC_TRY_COMPILE(
dnl    [],
dnl    [],
dnl    [
dnl      AC_MSG_RESULT(yes)
dnl    ],
dnl    [
dnl      AC_MSG_ERROR(invalid combination of flags!)
dnl      exit 1;
dnl    ])
dnl  
dnl  FFLAGS="$F77FLAGSO"
dnl  AC_MSG_CHECKING(for consistency of F77FLAGSO flags)
dnl  AC_TRY_COMPILE(
dnl    [],
dnl    [],
dnl    [
dnl      AC_MSG_RESULT(yes)
dnl    ],
dnl    [
dnl      AC_MSG_ERROR(invalid combination of flags!)
dnl      exit 1;
dnl    ])
])




dnl -------------------------------------------------------------
dnl Check for KDOC.
dnl
dnl Usage: DEAL_II_CHECK_KDOC
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_KDOC, dnl
[
  dnl    Find the kdoc directory for documentation. kdoc2 is in
  dnl    the contrib directory, but you might want another one
  AC_ARG_WITH(kdoc,
  [  --with-kdoc=DIR use kdoc installed in DIR],
      kdocdir=$withval,
      kdocdir=${DEAL2_DIR}/contrib/kdoc/bin)
  AC_MSG_CHECKING(for kdoc)

  dnl lets see whether the file exists if not the default was taken
  if test "$kdocdir" != ${DEAL2_DIR}/contrib/kdoc/bin ; then
    if test -r $kdocdir/kdoc ; then
      AC_MSG_RESULT(found)
    else
      AC_MSG_RESULT(not found)
      AC_MSG_ERROR(Invalid kdoc path $kdocdir/kdoc)
    fi
  
    if test -r "$kdocdir/Version" ; then
      kdocversion=`cat $kdocdir/Version` ;
    else
      kdocversion=1;
    fi
  else
    kdocversion=`cat ${DEAL2_DIR}/contrib/kdoc/src/Version`
    AC_MSG_RESULT(using default version $kdocversion)
  fi
])



dnl -------------------------------------------------------------
dnl Check for Doc++.
dnl
dnl Usage: DEAL_II_CHECK_DOCXX
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_DOCXX, dnl
[
  AC_ARG_WITH(docxx,
  [  --with-docxx=PATH use the doc++ executable pointed to by PATH],
      docxx=$withval,
      docxx=to-be-determined)
  if test "$docxx" = to-be-determined ; then
    AC_PATH_PROG(docxx,"doc++")
  else
    AC_MSG_CHECKING(for doc++)
    if test -x "$docxx" ; then
      AC_MSG_RESULT(yes)
    else
      AC_MSG_RESULT(no)
      docxx=
    fi
  fi
])

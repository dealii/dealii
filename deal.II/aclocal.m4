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
      *3.3*)
  	AC_MSG_RESULT(C++ compiler is gcc-3.3)
  	GXX_VERSION=gcc3.3
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
  
          dnl Intel's ECC C++ compiler for the Itanium?
          is_intel_ecc="`($CXX -V 2>&1) | grep 'Intel(R) C++ Itanium(TM) Compiler'`"
          if test "x$is_intel_ecc" != "x" ; then
            AC_MSG_RESULT(C++ compiler is Intel Itanium ECC)
            GXX_VERSION=intel_ecc
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
  
  	          dnl Portland Group C++?
  	          is_pgcc="`($CXX -V 2>&1) | grep 'Portland Group'`"
  	          if test "x$is_pgcc" != "x" ; then
  	            AC_MSG_RESULT(C++ compiler is Portland Group C++)
  	            GXX_VERSION=portland_group
  	          else
  
  	            dnl HP aCC?
  	            is_aCC="`($CXX -V 2>&1) | grep 'aCC'`"
  	            if test "x$is_aCC" != "x" ; then
  	              AC_MSG_RESULT(C++ compiler is HP aCC)
  	              GXX_VERSION=hp_aCC
  	            else
  
  	              dnl KAI C++? It seems as if the documented options
		      dnl -V and --version are not always supported, so give
	              dnl the whole thing a second try by looking for /KCC/
	 	      dnl somewhere in the paths that are output by -v. This
	              dnl is risky business, since this combination of
	              dnl characters might appear on other installations
                      dnl of other compilers as well, so put this test to
                      dnl the very end of the detection chain for the
                      dnl various compilers
  	              is_kai_cc="`($CXX --version 2>&1) | grep 'KAI C++'`"
  	              is_kai_cc="$is_kai_cc`($CXX -v 2>&1) | grep /KCC/`"
  	              if test "x$is_kai_cc" != "x" ; then
  	                AC_MSG_RESULT(C++ compiler is KAI C++)
  	                GXX_VERSION=kai_cc
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
    CXXFLAGSG="$CXXFLAGS -DDEBUG -pedantic -Wall -W -Wpointer-arith -Wwrite-strings -Wmissing-prototypes -Winline -Woverloaded-virtual -Wstrict-prototypes -Wsynth -Wsign-compare -Wconversion -Wswitch -ftemplate-depth-32"
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
                CXXFLAGSG="$CXXFLAGSG -fno-vtable-thunks"
                CXXFLAGSO="$CXXFLAGSO -fno-vtable-thunks"
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

    dnl In order to link shared libraries, almost all versions of gcc can
    dnl use CXX, i.e. the C++ compiler. The exception is gcc2.95, for which
    dnl we have to use the C compiler, unless we want to get linker errors
    SHLIBLD="$CXX"
    if test "$GXX_VERSION" = "gcc2.95"; then
      SHLIBLD="$CC"
    fi
  
  else
    dnl Non-gcc compilers. By default, use the C++ compiler also for linking
    dnl shared libraries. If some compiler cannot do that and needs something
    dnl different, then this must be specified in the respective section
    dnl below, overriding this define:
    SHLIBLD="$CXX"
  
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
          dnl #424: `extra ";" ignored'
          dnl #525: `type "DataOutBase::DataOutBase" is an inaccessible type
          dnl       (allowed for compatibility)' (I don't understand what the
          dnl       compiler means)
	  dnl
          dnl Note: we would really like to use  -ansi -Xc, since that
	  dnl is _very_ picky about standard C++, and is thus very efficient
          dnl in detecting slight standard violations, but these flags are
          dnl also very efficient in crashing the compiler (it generates a
          dnl segfault)
          CXXFLAGSG="$CXXFLAGS -Kc++eh -Krtti -w1 -wd175 -wd525 -wd327 -wd424 -DDEBUG -inline_debug_info"
          CXXFLAGSO="$CXXFLAGS -Kc++eh -Krtti -O2 -tpp6 -axiMK -ip -unroll -w0 -wd424"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-KPIC -shared"
          ;;
  
      intel_ecc)
          dnl Disable some compiler warnings, as they often are wrong on
          dnl our code:
          dnl #175: `subscript out of range' (doesn't take into account that
          dnl       some code is only reachable for some dimensions)
          dnl #327: `NULL reference is not allowed' (this happens when we
          dnl       write "*static_cast<double*>(0)" or some such thing,
          dnl       which we do to create invalid references)
          dnl #424: `extra ";" ignored'
          dnl #525: `type "DataOutBase::DataOutBase" is an inaccessible type
          dnl       (allowed for compatibility)' (I don't understand what the
          dnl       compiler means)
	  dnl
          dnl Note: we would really like to use  -ansi -Xc, since that
	  dnl is _very_ picky about standard C++, and is thus very efficient
          dnl in detecting slight standard violations, but these flags are
          dnl also very efficient in crashing the compiler (it generates a
          dnl segfault)
          CXXFLAGSG="$CXXFLAGS -Kc++eh -Krtti -w1 -wd175 -wd525 -wd327 -wd424 -DDEBUG -inline_debug_info"
          CXXFLAGSO="$CXXFLAGS -Kc++eh -Krtti -O2 -ip -unroll -w0 -wd424"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-KPIC -shared"
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
  
          CXXFLAGSPIC="-shared"
          LDFLAGSPIC="-shared"
          ;;
  
      sun_workshop | sun_forte)
          CXXFLAGSG="$CXXFLAGS -DDEBUG -w"
          CXXFLAGSO="$CXXFLAGS -w"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-G"
          ;;
  
      portland_group)
	  dnl Suppress warnings:
          dnl #111: "Statement unreachable": we use return statements
	  dnl       occasionally after case-switches where you cannot
 	  dnl       fall though, but other compilers cometimes complain
	  dnl       that the function might not return with a value, if
          dnl       it can't figure out that the function always uses
          dnl       one case. Also: a return statement after a failing
          dnl       assertion
          dnl #177: "function declared but not used": might happen with
          dnl       templates and conditional compilation
 	  dnl #175: "out-of-bounds array indices": the same reason as
	  dnl       for Compaq cxx
	  CXXFLAGSG="$CXXFLAGS -DDEBUG -g --display_error_number --diag_suppress 111 --diag_suppress 177 --diag_suppress 175"
          CXXFLAGSO="$CXXFLAGS -fast -O2 --display_error_number --diag_suppress 111 --diag_suppress 177 --diag_suppress 175"
          CXXFLAGSPIC="-Kpic"
          ;;

      kai_cc)
          CXXFLAGSG="$CXXFLAGS --strict -D__KAI_STRICT --max_pending_instantiations 32 --display_error_number -g +K0 --no_implicit_typename"
          CXXFLAGSO="$CXXFLAGS +K3 -O2 --abstract_float --abstract_pointer -w --display_error_number --max_pending_instantiations 32 --display_error_number"
          CXXFLAGSPIC="-fPIC"
	  ;;

      hp_aCC)
	  dnl ??? disable warning 655 (about all-inlined functions) which
	  dnl triggers for each and every of our DeclExceptionX calls ???
          CXXFLAGSG="$CXXFLAGS -g1 -AA +p"
          CXXFLAGSO="$CXXFLAGS -z +O2 -AA"
          CXXFLAGSPIC="+Z"
	  # for linking shared libs, -b is also necessary...
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
        *3.3*)
  	  AC_MSG_RESULT(F77 compiler is gcc-3.3)
  	  F77_VERSION=gcc3.3
  	  ;;
        *2.4* | *2.5* | *2.6* | *2.7* | *2.8*)
  	  dnl These compilers are too old to support a useful subset
  	  dnl of modern C++, so we don't support them. gcc2.7.2 is 
	  dnl probably the only one that around on reasonably modern
	  dnl systems these times, but maybe someone tries to run
	  dnl deal.II on really old systems?
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
        if test -n "`echo $F77_VERSION_STRING | grep \"WorkShop Compilers\"`" \
                -o \
                -n "`echo $F77_VERSION_STRING | grep \"Sun WorkShop\"`" \
                -o \
                -n "`echo $F77_VERSION_STRING | grep \"Forte Developer\"`" ; then
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
            F77_VERSION="UnknownF77"
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
    egcs-1.1 | gcc2.95 | gcc2.96 | gcc2.97 | gcc3.[[0123]])
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

    UnknownF77)
	dnl Disable unknown FORTRAN compiler.
	dnl Allows configure to finish,
	dnl but disables compiling FORTRAN code
	F77="UnknownF77"
        F77FLAGSG="$FFLAGS -g"
        F77FLAGSO="$FFLAGS -O2"
        F77LIBS="$F77LIBS"  

        F77FLAGSPIC="unknown!"
	AC_MSG_RESULT(Unknown FORTRAN compiler has been disabled!)
        ;;

        dnl Keep this line just in case we change default options
	dnl back to error message
    *)
	AC_MSG_ERROR(No compiler options for F77 compiler
                     "$F77_VERSION" specified: modification of aclocal.m4 necessary)
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
dnl   (don't use this for gcc3.1 or later)
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
  
      CXXFLAGSG="$CXXFLAGSG -D_REENTRANT"
      CXXFLAGSO="$CXXFLAGSO -D_REENTRANT"
      if test "$GXX_VERSION" = "gcc2.95" \
              -o "$GXX_VERSION" = "gcc2.96" \
	      -o "$GXX_VERSION" = "gcc3.0" ; then
        CXXFLAGSG="$CXXFLAGSG -D__USE_MALLOC"
        CXXFLAGSO="$CXXFLAGSO -D__USE_MALLOC"
      fi
    else
      case "$GXX_VERSION" in
	ibm_xlc)
            CXXFLAGSG="$CXXFLAGSG -threaded"  
            CXXFLAGSO="$CXXFLAGSO -threaded"
            ;;

	compaq_cxx)
            CXXFLAGSG="$CXXFLAGSG -pthread"  
            CXXFLAGSO="$CXXFLAGSO -pthread"
	    ;;

        intel_icc)
            CXXFLAGSG="$CXXFLAGSG"  
            CXXFLAGSO="$CXXFLAGSO -parallel"
	    ;;

	*)
            AC_MSG_ERROR(No threading compiler options for this C++ compiler
                         specified at present)
            exit 1
	    ;;
      esac
    fi
  fi
])



dnl -------------------------------------------------------------
dnl Test which library the MT code shall use to support threads.
dnl We support either POSIX or ACE. If requested, then set
dnl $withmultithreading to either POSIX or to the path to ACE.
dnl
dnl Usage:
dnl   DEAL_II_CHECK_USE_MT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_USE_MT, dnl
[
  AC_ARG_WITH(multithreading,
  [  --with-multithreading=name If name==posix, or no name given, then use
                                POSIX threads, otherwise assume use of ACE 
                                and use given argument as path to ACE],
      withmultithreading=$withval,
      withmultithreading=no)

dnl Default (i.e. no arg) means POSIX
  if test "x$withmultithreading" = "xyes" ; then
    withmultithreading=posix
  fi

  if test "x$withmultithreading" != "xno" ; then
    if test "x$withmultithreading" = "xposix" ; then
      DEAL_II_CHECK_POSIX_THREAD_FUNCTIONS
      AC_DEFINE(DEAL_II_USE_MT_POSIX, 1,
                [Defined if multi-threading is to be achieved by using
                 the POSIX functions])
    else
      AC_MSG_CHECKING(for ACE)
      if test -d "$withmultithreading" ; then
        AC_MSG_RESULT(found)
      else
        AC_MSG_RESULT(not found)
        AC_MSG_ERROR(Invalid ACE path)
      fi
      AC_DEFINE(DEAL_II_USE_MT_ACE, 1,
                [Defined if multi-threading is to be achieved by using
                 the ACE library])
      DEAL_II_SET_ACE_FLAGS
    fi
  
    AC_DEFINE(DEAL_II_USE_MT, 1, 
              [Flag indicating whether the library shall be compiled for
               multithreaded applications])
  fi
])



dnl -------------------------------------------------------------
dnl Check whether the POSIX functions needed for multi-threading
dnl are available on this system
dnl 
dnl Usage: DEAL_II_CHECK_POSIX_THREAD_FUNCTIONS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_POSIX_THREAD_FUNCTIONS, dnl
[
  AC_MSG_CHECKING(for posix thread functions)
  AC_LANG(C++)
  AC_TRY_COMPILE(
   [
#	include <pthread.h>
   ],
   [
	pthread_t         p;
	pthread_create (&p, 0, 0, 0);
	pthread_join (p, 0);
   ],
   [
	AC_MSG_RESULT(ok)
   ],
   [
	AC_MSG_ERROR(not found)
   ])

  AC_MSG_CHECKING(for posix thread mutex functions)
  AC_LANG(C++)
  AC_TRY_COMPILE(
   [
#	include <pthread.h>
   ],
   [
	pthread_mutex_t   pm;
	pthread_mutex_init (&pm, 0);
	pthread_mutex_lock (&pm);
	pthread_mutex_unlock (&pm);
	pthread_mutex_destroy (&pm);
   ],
   [
	AC_MSG_RESULT(ok)
   ],
   [
	AC_MSG_ERROR(not found)
   ])

  AC_MSG_CHECKING(for posix thread barrier functions)
  AC_LANG(C++)
  AC_TRY_COMPILE(
   [
#	include <pthread.h>
   ],
   [
	pthread_barrier_t pb;
	pthread_barrier_init (&pb, 0, 1);
	pthread_barrier_wait (&pb);
	pthread_barrier_destroy (&pb);
   ],
   [
	AC_MSG_RESULT(ok)
	x=0
   ],
   [
	AC_MSG_RESULT(not found. barriers will not be supported)
	x=1
   ])
   if test "x$x" = "x1" ; then
     AC_DEFINE(DEAL_II_USE_MT_POSIX_NO_BARRIERS, 1,
	       [Defined if POSIX is supported but not the newer POSIX
                barrier functions. Barriers will then not work in
                the library, but the other threading functionality
                is available.])
  fi
])



dnl -------------------------------------------------------------
dnl If the MT code shall be used through ACE, we might have to 
dnl adjust the compiler flags: Possibly remove -pedantic from
dnl compiler flags again, since ACE yields a lot of error messages
dnl with this flag.
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
dnl ourselves in those files that use this function. The question whether
dnl we have to do so is controlled by a preprocessor variable.
dnl
dnl If the function is not properly declared, then set the preprocessor
dnl variable NO_HAVE_GETRUSAGE
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
	AC_DEFINE(NO_HAVE_GETRUSAGE, 1, 
                  [On SunOS 4.x, the getrusage() function exists, but 
		   is not declared in the respective header file
		   <resource.h>, as one would think when reading the
		   man pages. Then we have to declare this function
		   ourselves in those files that use this function.
		   The question whether we have to do so is controlled
		   by the preprocessor variable.])
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
unsigned int seed = 0;
int i=rand_r(&i);
	],
	[
	  AC_MSG_RESULT(yes)
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
dnl Gcc and some other compilers have __PRETTY_FUNCTION__, showing
dnl an unmangled version of the function we are presently in,
dnl while __FUNCTION__ (or __func__ in ISO C99) simply give the
dnl function name which would not include the arguments of that
dnl function, leading to problems in C++ with overloaded function
dnl names.
dnl
dnl If __PRETTY_FUNCTION__ is not available, try to find out whether
dnl __func__ is available and use the preprocessor to set the first
dnl thing to the second. If this is also not the case, then set it 
dnl to something indicating non-availability.
dnl
dnl Usage: DEAL_II_HAVE_PRETTY_FUNCTION
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_PRETTY_FUNCTION, dnl
[
  AC_MSG_CHECKING(for __PRETTY_FUNCTION__)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#	include <iostream>
    ],
    [
	std::cout << __PRETTY_FUNCTION__ << std::endl;
    ],
    [
      AC_MSG_RESULT(available)
    ],
    [
      AC_MSG_RESULT(not available)
      AC_MSG_CHECKING(for __func__)
      AC_TRY_COMPILE(
        [
#	include <iostream>
        ],
        [
	  std::cout << __func__ << std::endl;
        ],
        [
          AC_MSG_RESULT(available)
	  x=__func__
    	],
        [
          AC_MSG_RESULT(not available)
	  x="\"(not available)\""
        ]) 
      AC_DEFINE_UNQUOTED(__PRETTY_FUNCTION__, $x, 
                [If already available, do not define at all. Otherwise, define
                 to __func__ if that is available. In all other cases,
                 indicate that no information about the present function
                 is available for this compiler.])
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
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(XLC_WORK_AROUND_STD_BUG, 1, 
                [Define if we have to work around a bug in IBM's xlC compiler.
                 See the aclocal.m4 file in the top-level directory for a
                 description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl Sun's Forte compiler (at least up to the Version 7 Early Access
dnl release) has a problem with the following code, when compiling
dnl with debug output:
dnl 
dnl /* ---------------------------------------------------------- */
dnl /* Internal compiler error in abi2_mangler::entity_expression */
dnl /* when compiled with -g.                                     */
dnl template < int dim > struct T {
dnl     typedef T<dim-1> SubT;
dnl     T (SubT);
dnl };
dnl 
dnl template <int dim> T<dim>::T (SubT) {};
dnl 
dnl template class T<3> ;
dnl /* ---------------------------------------------------------- */
dnl
dnl The compiler gets an internal compiler error, so we work around
dnl this problem by a really evil hack in the sources.
dnl
dnl Usage: DEAL_II_CHECK_LOCAL_TYPEDEF_COMP
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_LOCAL_TYPEDEF_COMP, dnl
[
  AC_MSG_CHECKING(for local computed template typedef bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -g"
  AC_TRY_COMPILE(
    [
	template < int dim > struct T {
	    typedef T<dim-1> SubT;
	    T (SubT);
	};

	template <int dim> T<dim>::T (SubT) {};

	template class T<3> ;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_LOCAL_TYPEDEF_COMP_WORKAROUND, 1, 
                [Define if we have to work around a bug in Sun's Forte compiler.
                 See the aclocal.m4 file in the top-level directory for a
                 description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl Sun's Forte compiler (at least up to the Version 7 Early Access
dnl release) have a problem with the following code, when compiling
dnl with debug output:
dnl 
dnl /* ----------------------------------------------- */
dnl /* Problem 14: Access control. Friendship is not   */
dnl /* granted although explicitly declared.           */
dnl template <int N, int M> class T      {    void bar ();  };
dnl 
dnl template <int M>        class T<1,M> { 
dnl   private:
dnl     static int i;
dnl     template <int N1, int N2> friend class T;
dnl };
dnl 
dnl template <int N,int M> void T<N,M>::bar () { 
dnl   T<N-1,M>::i; 
dnl };
dnl 
dnl template class T<2,1> ; 
dnl /* ---------------------------------------------------------- */
dnl
dnl The compiler does not allow access to T<1,1>::i, although the
dnl accessing class is explicitly marked friend.
dnl
dnl Usage: DEAL_II_CHECK_TEMPLATE_SPEC_ACCESS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPLATE_SPEC_ACCESS, dnl
[
  AC_MSG_CHECKING(for partially specialized template access control bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	template <int N, int M> class T      {    void bar ();  };

	template <int M>        class T<1,M> { 
	  private:
	    static int i;
	    template <int N1, int N2> friend class T;
	};

	template <int N,int M> void T<N,M>::bar () { 
	  T<N-1,M>::i; 
	};

	template class T<2,1> ; 
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPLATE_SPEC_ACCESS_WORKAROUND, 1, 
                [Define if we have to work around a bug in Sun's Forte compiler.
                 See the aclocal.m4 file in the top-level directory for a
                 description of this bug.])
    ])
])




dnl -------------------------------------------------------------
dnl Versions of GCC before 3.0 had a problem with the explicit
dnl instantiation of member templates when the member was in fact
dnl an operator. In that case, they needed the "template" keyword,
dnl which is actually not allowed at this place. Test case is
dnl 
dnl /* ----------------------------------------------- */
dnl struct X
dnl {
dnl     template <typename T2>
dnl     X operator = (T2 &){};
dnl };
dnl 
dnl template X X::operator=<float> (float &);
dnl /* ---------------------------------------------------------- */
dnl
dnl The compiler only groks this if the "operator=" is prepended
dnl by "template". We detect this, and either set the 
dnl DEAL_II_MEMBER_OP_TEMPLATE_INST to "template" or nothing, so
dnl that it gets expanded to the right string needed in this place.
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_OP_TEMPLATE_INST
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_OP_TEMPLATE_INST, dnl
[
  AC_MSG_CHECKING(for template member operator instantiation bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	struct X
	{
	    template <typename T2>
	    X operator = (T2 &){};
	};

	template X X::operator=<float> (float &);
    ],
    [],
    [
      AC_MSG_RESULT(no)
      x=""
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      x="template"
    ])
  AC_DEFINE_UNQUOTED(DEAL_II_MEMBER_OP_TEMPLATE_INST, $x, 
                     [Define if we have to work around a bug in gcc with
                      explicitly instantiating template member operators.
                      See the aclocal.m4 file in the top-level directory
                      for a description of this bug.])
])



dnl -------------------------------------------------------------
dnl Versions of GCC before 3.0 had a problem with the following
dnl code:
dnl 
dnl /* ----------------------------------------------- */
dnl namespace NS {
dnl   template <typename T>  class C  {
dnl       template <typename N> friend class C;
dnl   };
dnl };
dnl /* ----------------------------------------------- */
dnl
dnl This is fixed with gcc at least in snapshots before version 3.1,
dnl but the following bug remains:
dnl
dnl /* ----------------------------------------------- */
dnl namespace NS {  template <typename number> class C;  };
dnl
dnl template <typename T> class X {
dnl   template <typename N> friend class NS::C;
dnl };
dnl
dnl template class X<int>;
dnl /* ----------------------------------------------- */
dnl
dnl The compiler gets an internal error for these cases. Since we need this
dnl construct at various places, we check for it and if the compiler
dnl dies, we use a workaround that is non-ISO C++ but works for these
dnl compilers.
dnl
dnl Usage: DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for template friend in namespace bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	namespace NS {
	  template <typename T>  class C  {
	      C(const C<T>&);
	      template <typename N> friend class C;
	  };
	};

	namespace NS2 {  template <typename number> class C;  };

	template <typename T> class X {
	    template <typename N> friend class NS2::C;
	    template <typename N> friend class NS::C;
	};

	template class X<int>;

	namespace NS {
  	  template<typename T>
	  inline C<T>::C(const C<T>&)
          {};
	};

    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE_UNQUOTED(DEAL_II_NAMESP_TEMPL_FRIEND_BUG, 1, 
                         [Define if we have to work around a bug in gcc with
                          marking all instances of a template class as friends
		          to this class if the class is inside a namespace.
                          See the aclocal.m4 file in the top-level directory
                          for a description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl In some cases, we would like to name partial specializations
dnl as friends. However, the standard forbids us to do so. But
dnl then, we can declare the general template as a friend, and
dnl at least gcc extends the friendship to all specializations
dnl of the templates, which is not what the standard says.
dnl
dnl With other compilers, most notably cxx, this does not work.
dnl In this case, we can make individual specializations friends,
dnl which in turn gcc rejects. So check, whether this is possible.
dnl
dnl The respective clause in the standard is 14.5.3.1, which gives
dnl this example:
dnl   template<class T> class task {
dnl     friend class task<int>;
dnl   };
dnl
dnl
dnl Usage: DEAL_II_CHECK_TEMPL_SPEC_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPL_SPEC_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for template specialization friend bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
template <int N, typename T> class X;
template <typename T>        class X<1,T>;

template <typename P> class Y {
    static int i;
    template <int N, typename T> friend class X;
    friend class X<1,P>;
};

template <typename T> class X<1,T> {
    X () { Y<T>::i; };     // access private field
};

template class X<1,int>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE_UNQUOTED(DEAL_II_TEMPL_SPEC_FRIEND_BUG, 1, 
                         [Define if we have to work around a bug with some
			  compilers that will not allow us to specify a
			  fully specialized class of a template as a friend.
                          See the aclocal.m4 file in the top-level directory
                          for a description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl Intel's ICC 5.0.1 compiler has the following problem: it won't
dnl compile this code:
dnl ---------------------------------
dnl template <class Class> void encapsulate (void (Class::*fun_ptr)());
dnl struct X {
dnl  void bar () const;
dnl };
dnl void foo () {
dnl   encapsulate(&X::bar);
dnl };
dnl ---------------------------------
dnl It complains about not finding an instance for the encapsulate function.
dnl The problem is due to the fact that the function we take the address of
dnl is constant, i.e. it should match the template with Class="const X",
dnl but apparently doesn't. Unfortunately, we rely on such code in the
dnl Threads mechanisms. So detect this bug and if present work around it
dnl by providing a second set of encapsulate functions for constant 
dnl functions.
dnl
dnl Usage: DEAL_II_CHECK_TEMPL_CONST_MEM_PTR_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPL_CONST_MEM_PTR_BUG, dnl
[
  AC_MSG_CHECKING(for templates and pointers to const member functions bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	template <class Class> void encapsulate (void (Class::*fun_ptr)());

	struct X {
	  void bar () const;
	};

	void foo () {
	  encapsulate(&X::bar);
	};
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPL_CONST_MEM_PTR_BUG, 1, 
                     [Define if we have to work around a bug in Intel's icc
                      compiler in which the compiler refuses to consider a
                      template when given a pointer to a constant member function.
                      See the aclocal.m4 file in the top-level directory
                      for a description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl More Intel ICC 5.0.1 compiler lossage with pointers to constant
dnl member functions
dnl ---------------------------------
dnl template <class Class> struct Y {
dnl   typedef void (Class::*FunPtr) ();
dnl   FunPtr fun_ptr;
dnl   Class *object;
dnl   void foo () {
dnl     (object->*fun_ptr)();
dnl   };
dnl };
dnl
dnl struct X {};
dnl
dnl template class Y<const X>;
dnl ---------------------------------
dnl Again, it does not store the information that the member function pointer
dnl refers to a constant object, so complains that the object is constant,
dnl but the member function is not.
dnl
dnl Usage: DEAL_II_CHECK_CONST_MEM_FUN_PTR_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CONST_MEM_FUN_PTR_BUG, dnl
[
  AC_MSG_CHECKING(for const member function pointers bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	template <class Class> struct Y {
	  typedef void (Class::*FunPtr) ();
	  FunPtr fun_ptr;
	  Class *object;
	  void foo () {
	    (object->*fun_ptr)();
	  };
	};

	struct X {};

	template class Y<const X>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_CONST_MEM_FUN_PTR_BUG, 1, 
                     [Define if we have to work around another bug in Intel's icc
                      compiler in which the compiler refuses to store the const
                      attribute at a member function pointer with a constant class.
                      See the aclocal.m4 file in the top-level directory
                      for a description of this bug.])
    ])
])



dnl -------------------------------------------------------------
dnl DEC/Compaq's cxx compiler does not want us to implement
dnl virtual functions that were declared abstract before. We do
dnl this with the destructor of the Function class, since we want
dnl to avoid people making objects of that class, but all functions
dnl have default implementations, so the class is not abstract
dnl without that. Since every derived class has a destructor, it
dnl is sufficient to mark the destructor =0.
dnl
dnl Unfortunately, cxx refuses to grok that. It sees the respective
dnl function in the .cc file, but does not instantiate it, leading
dnl to linker errors. Thus, check this misfeature, and if present
dnl simply do not mark the function abstract for this particular
dnl compiler.
dnl
dnl Usage: DEAL_II_CHECK_IMPLEMENTED_PURE_FUNCTION_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_IMPLEMENTED_PURE_FUNCTION_BUG, dnl
[
  AC_MSG_CHECKING(for bug with implementing pure functions)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	template <int dim>
	struct Function
	{
	  public:
	    virtual ~Function () = 0;
	};

	template <int dim>
	Function<dim>::~Function ()
	{};

	template class Function<1>;
	template Function<1>::~Function();
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_IMPLEMENTED_PURE_FUNCTION_BUG, 1, 
                     [Defined if the compiler refuses to compile the definition
		      of a function that was previously declared abstract.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc 2.95 dies on this construct:
dnl -----------------------------
dnl template <int dim> struct TT { typedef int type; };
dnl
dnl template <template <int> class T> struct X {
dnl     typedef typename T<1>::type type;
dnl     void foo (type t);
dnl };
dnl
dnl template <template <int> class T>
dnl void X<T>::foo (type t) {};
dnl
dnl template struct X<TT>;
dnl -----------------------------
dnl
dnl We work around this problem, if we encounter it.
dnl
dnl Usage: DEAL_II_CHECK_TEMPLATE_TEMPLATE_TYPEDEF_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPLATE_TEMPLATE_TYPEDEF_BUG, dnl
[
  AC_MSG_CHECKING(for template template typedef bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	template <int dim> struct TT { typedef int type; };

	template <template <int> class T> struct X {
	    typedef typename T<1>::type type;
	    void foo (type t);
	};

	template <template <int> class T>
	void X<T>::foo (type t) {};

	template struct X<TT>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPLATE_TEMPLATE_TYPEDEF_BUG, 1, 
                     [Defined if the compiler refuses to allow a typedef
                      to a template template class template parameter. For
                      the exact failure mode, look at aclocal.m4 in the
                      top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc 2.95 does not correctly implement the resolution of 
dnl defect report #45 to the C++ standard (see
dnl http://anubis.dkuug.dk/jtc1/sc22/wg21/docs/cwg_active.html#45).
dnl try to detect this, and set a flag correspondingly. in short,
dnl the DR says that this is allowed, since member classes are
dnl implicitly also friends:
dnl -----------------------------
dnl struct X {
dnl   private:
dnl     static int f();
dnl     
dnl     struct Y {
dnl         int g() { return f(); };
dnl     };
dnl };
dnl -----------------------------
dnl
dnl We work around this problem, if we encounter it.
dnl
dnl Usage: DEAL_II_CHECK_NESTED_CLASS_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NESTED_CLASS_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for nested classes are implicit friends bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	struct X {
	  private:
	    static int f();
	    
	    struct Y {
	        int g() { return f(); };
	    };
	};
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_NESTED_CLASS_FRIEND_BUG, 1, 
                     [Defined if the compiler does not properly implement
                      the resolution of defect report #45 to the C++
                      standard, which makes nested types implicit friends
                      of the enclosing class.])
    ])
])



dnl -------------------------------------------------------------
dnl Many compilers get this wrong (see Section 14.7.3.1, number (4)):
dnl ---------------------------------
dnl   template <int dim> struct T {
dnl     static const int i;
dnl   };
dnl
dnl   template <> const int T<1>::i;
dnl   template <> const int T<1>::i = 1;
dnl ---------------------------------
dnl First, by Section 14.7.3.14 of the standard, the first template<>
dnl line must necessarily be the _declaration_ of a specialization,
dnl and the second is then its definition. There is therefore no
dnl reason to report a doubly defined variable (Intel ICC 6.0), or
dnl to choke on these lines at all (Sun Forte)
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_VARIALIZATION_SPEC_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_VAR_SPECIALIZATION_BUG, dnl
[
  AC_MSG_CHECKING(for member variable specialization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	template <int dim> struct T {
	    static const int i;
	};

	template <> const int T<1>::i;
	template <> const int T<1>::i = 1;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG, 1, 
                     [Defined if the compiler refuses to allow the
                      explicit specialization of static member 
                      variables. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc3.1 (and maybe later compilers) has a bug with long double
dnl and optimization (see code below), when compiling on Sparc
dnl machines. Since it affects only one platform and one compiler,
dnl we take the liberty to disable the function in which the problem
dnl occurs (Polynomial::shift in base/source/polynomial.cc), since 
dnl this is a function that is rarely used anyway.
dnl
dnl For more information: the bug just described is reported to
dnl the gcc project under number 7335.
dnl
dnl Usage: DEAL_II_CHECK_LONG_DOUBLE_LOOP_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_LONG_DOUBLE_LOOP_BUG, dnl
[
  AC_MSG_CHECKING(for long double optimization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSO"
  AC_TRY_COMPILE(
    [
	double* copy(long double* first, long double* last, double* result)
	{
	  int n;
	  for (n = last - first; n > 0; --n) {
	    *result = *first;
	    ++first;
	    ++result;
	  }
	  return result;
	}

	void f()
	{
	  long double *p1, *p2;
	  double *p3;
	  copy (p1, p2, p3);
	  p3 = copy (p1, p2, p3);
	};
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. disabling respective functions)
      AC_DEFINE(DEAL_II_LONG_DOUBLE_LOOP_BUG, 1, 
                     [Defined if the compiler gets an internal compiler
                      upon some code involving long doubles, and with
                      optimization. For the details, look at
                      aclocal.m4 in the top-level directory.])
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
dnl When compiling with ACE thread support, there are many constructs
dnl that are not allowed in C++, or that yield warnings when compiling with
dnl -pedantic. Check this, and if that is the case, set the variables
dnl $deal_II_ace_remove_ansi and $deal_II_ace_remove_pedantic to "yes".
dnl
dnl Usage: DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ACE_FORBIDDEN_FLAGS, dnl
[
  AC_MSG_CHECKING(whether compilation with ACE disallows flags)
  AC_LANG(C++)
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

  if test $deal_II_ace_remove_pedantic = "no" ; then 
    AC_MSG_RESULT(no)
  else
    AC_MSG_RESULT(-pedantic)
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



dnl -------------------------------------------------------------
dnl On some Cygwin systems, a system header file includes this
dnl preprocessor define:
dnl   #define quad quad_t
dnl This is of course silly, but beyond that it also hurts as
dnl since we have member functions and variables with that name
dnl and we get compile errors depending or not we have this
dnl particular header file included. 
dnl
dnl Fortunately, the define is only active is _POSIX_SOURCE is
dnl not set, so check for this define, and if necessary set
dnl this flag. We check on all systems, since maybe there are 
dnl other such systems elsewhere...
dnl
dnl Usage: DEAL_II_CHECK_QUAD_DEFINE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_QUAD_DEFINE, dnl
[
  AC_MSG_CHECKING(for quad vs. quad_t define)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <sys/types.h>
#if defined(quad)
    no good system;
#endif
    ],
    [
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. working around)
      CXXFLAGSG="$CXXFLAGSG -D_POSIX_SOURCE"
      CXXFLAGSO="$CXXFLAGSO -D_POSIX_SOURCE"
    ])
])



dnl -------------------------------------------------------------
dnl On DEC OSF1, when we specify "-std strict_ansi", we can include
dnl errno.h, and still not get the definition of the error codes
dnl such as EINTR, EPIPE, etc.
dnl
dnl In this case, use a workaround by explicitly including 
dnl /usr/include/errno.h, instead of just errno.h, which the compiler
dnl maps to one of its own C/C++ compatibility headers, which only
dnl define 3 error codes (for reasons unknown)
dnl
dnl Usage: DEAL_II_CHECK_ERROR_CODES_DEFINITION
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ERROR_CODES_DEFINITION, dnl
[
  AC_MSG_CHECKING(for defintions of error codes in errno.h)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <iostream>
#include <errno.h>
using namespace std;
    ],
    [
    cout << EINTR << endl;
    ],
    [
      AC_MSG_RESULT(yes)
    ],
    [
      AC_MSG_RESULT(no. working around)
      AC_DEFINE(DEAL_II_USE_DIRECT_ERRNO_H, 1, 
                [Define if the compiler provides a <errno.g> header file
		 which does not define all error codes such as EINTR. In
		 that case, use the system include file at /usr/include
		 instead. There is probably a better way to do this, but
		 it is not apparent by looking at the C/C++ compatibility
                 header provided by the compiler.])
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

  if (test -r "$TECPLOT_LIBRARY_PATH" && \
      test -r "$TECPLOT_INCLUDE_PATH/TECIO.h") ; then
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

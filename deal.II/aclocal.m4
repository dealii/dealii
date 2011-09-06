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
dnl Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
dnl
dnl $Id$



dnl -------------------------------------------------------------
dnl Helper macros to add libpaths to LIBS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL, dnl
[
  LIBS="$LIBS $1"
])
AC_DEFUN(DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT, dnl
[
  LIBS="$1 $LIBS"
])
AC_DEFUN(DEAL_II_EXTERNAL_LIBS_SAVE_VAL, dnl   Not reentrant, of course
[
  OLD_LIBS="$LIBS"
])
AC_DEFUN(DEAL_II_EXTERNAL_LIBS_RESTORE_VAL, dnl
[
  LIBS="$OLD_LIBS"
])



dnl -------------------------------------------------------------
dnl Like AC_PATH_PROG, but do not discard arguments given to the
dnl program. In other words, while
dnl    AC_PATH_PROG(CXX, [g++ -pg])
dnl results in CXX=/usr/bin/g++, the result of the current
dnl macro would be CXX="/usr/bin/g++ -pg".
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_PATH_PROG, dnl
[
  dnl First get at the name and arguments of the program in $2. Do
  dnl so by having a loop over all components of $2 and putting the
  dnl components either into $testprog or into $testargs
  testprog=""
  testargs=""
  processingargs="no"
  for i in $2 ; do
    if test "$processingargs" = "no" ; then
      testprog="$i" ;
      processingargs="yes" ;
    else
      testargs="$testargs $i" ;
    fi
  done

  AC_PATH_PROG([$1],[$testprog])
  eval "$1=\"${$1} $testargs\""
])



dnl -------------------------------------------------------------
dnl Determine the C++ compiler in use. Return the name and possibly
dnl version of this compiler in GXX_VERSION.
dnl
dnl Usage: DEAL_II_DETERMINE_CXX_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_DETERMINE_CXX_BRAND, dnl
[
  dnl It used to be the case that AC_PROG_CXX could properly detect
  dnl whether we are working with gcc or not, but then Intel came across
  dnl the brilliant idea to disguise it's compiler as gcc by setting
  dnl GCC preprocessor variables for versions etc. So we have to figure
  dnl out again whether this is really gcc
  if test "$GXX" = "yes" ; then
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    if test "x$GXX_VERSION_STRING" = "x" ; then
      GXX=no
    fi
  fi

  dnl And because it is so convenient, the PathScale compiler also identifies
  dnl itself as GCC...
  if test "$GXX" = "yes" ; then
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "PathScale"`
    if test "x$GXX_VERSION_STRING" != "x" ; then
      GXX=no
    fi
  fi

  if test "$GXX" = yes ; then
    dnl find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`

    full_version=`echo "$GXX_VERSION_STRING" | perl -pi -e 's/.*version (\d\.\d\.\d).*/\1/g;'`
    GXX_BRAND=GNU
    GXX_VERSION=gcc`echo $full_version | perl -pi -e 's/(\d\.\d).*/\1/g;'`
    GXX_VERSION_DETAILED=gcc$full_version

    AC_MSG_RESULT([C++ compiler is $GXX_VERSION (subversion $GXX_VERSION_DETAILED)])


  else
    dnl Check other (non-gcc) compilers

    dnl Check for IBM xlC. For some reasons, depending on some environment
    dnl variables, moon position, and other reasons unknown to me, the
    dnl compiler displays different names in the first line of output, so
    dnl check various possibilities
    is_ibm_xlc="`($CXX -qversion 2>&1) | grep IBM`"
    if test "x$is_ibm_xlc" != "x"  ; then
      dnl Ah, this is IBM's C++ compiler. Unfortunately, we don't presently
      dnl know how to check the version number, so assume that is sufficiently
      dnl high...
      AC_MSG_RESULT(C++ compiler is IBM xlC)
      GXX_BRAND=IBM
      GXX_VERSION=ibm_xlc
      GXX_VERSION_DETAILED="$GXX_VERSION"
    else

      dnl Check whether we are dealing with the MIPSpro C++ compiler
      mips_pro="`($CXX -version 2>&1) | grep MIPSpro`"
      if test "x$mips_pro" != "x" ; then
        GXX_BRAND=MIPSpro
        case "$mips_pro" in
          *7.0* | *7.1* | *7.2* | *7.3*)
            dnl MIPSpro 7.3 does not support standard C++, therefore it is not
            dnl able to compile deal.II. Previous compiler versions neither.
            AC_MSG_RESULT(C++ compiler is $mips_pro)
            AC_MSG_ERROR(This compiler is not supported)
            GXX_VERSION=MIPSpro7.3
            GXX_VERSION_DETAILED="$GXX_VERSION"
            ;;
          *7.4)
            AC_MSG_RESULT(C++ compiler is MIPSpro compiler 7.4)
            AC_MSG_ERROR(This compiler is not supported. Use MIPSPro compiler 7.4x)
            GXX_VERSION=MIPSpro7.4
            GXX_VERSION_DETAILED="$GXX_VERSION"
            ;;
          *7.41* | *7.42* | *7.43* | *7.44*)
            AC_MSG_RESULT(C++ compiler is MIPSpro compiler 7.4x)
            GXX_VERSION=MIPSpro7.4x
            GXX_VERSION_DETAILED="$GXX_VERSION"
            ;;
          *"7.5"*)
            AC_MSG_RESULT(C++ compiler is MIPSpro compiler 7.5)
            GXX_VERSION=MIPSpro7.5
            GXX_VERSION_DETAILED="$GXX_VERSION"
            ;;
          *)
            AC_MSG_RESULT(C++ compiler is unknown version but accepted MIPSpro compiler version)
            GXX_VERSION=MIPSpro-other
            GXX_VERSION_DETAILED="$GXX_VERSION"
	    ;;
        esac
      else

        dnl Intel's ICC C++ compiler? On Linux, it uses -V, on Windows
	dnl it is -help
        is_intel_icc1="`($CXX -V 2>&1) | grep 'Intel'`"
        is_intel_icc2="`($CXX -help 2>&1) | grep 'Intel'`"
	is_intel_icc="$is_intel_icc1$is_intel_icc2"

	dnl When calling the Portland Group compiler, it also
	dnl outputs the string 'Intel' in its help text, so make
	dnl sure we're not confused
	is_pgi="`($CXX -V 2>&1) | grep 'Portland'`"

        if test "x$is_intel_icc" != "x" -a "x$is_pgi" = "x" ; then
	  GXX_BRAND=Intel
	  version_string="`($CXX -V 2>&1) | grep 'Version'` `($CXX -help 2>&1) | grep 'Version'`"
	  version5="`echo $version_string | grep 'Version 5'`"
	  version6="`echo $version_string | grep 'Version 6'`"
	  version7="`echo $version_string | grep 'Version 7'`"
	  version8="`echo $version_string | grep 'Version 8'`"
	  version9="`echo $version_string | grep 'Version 9'`"
	  version10="`echo $version_string | grep 'Version 10'`"
	  version11="`echo $version_string | grep 'Version 11'`"
          if test "x$version5" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-5)
            GXX_VERSION=intel_icc5
          else if test "x$version6" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-6)
            GXX_VERSION=intel_icc6
          else if test "x$version7" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-7)
            GXX_VERSION=intel_icc7
          else if test "x$version8" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-8)
            GXX_VERSION=intel_icc8
          else if test "x$version9" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-9)
            GXX_VERSION=intel_icc9
          else if test "x$version10" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-10)
            GXX_VERSION=intel_icc10
          else if test "x$version11" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-11)
            GXX_VERSION=intel_icc11
          else if test "x$version12" != "x" ; then
            AC_MSG_RESULT(C++ compiler is icc-12)
            GXX_VERSION=intel_icc12
          else
            AC_MSG_RESULT(C++ compiler is icc)
            GXX_VERSION=intel_icc
          fi fi fi fi fi fi fi fi
          GXX_VERSION_DETAILED="$GXX_VERSION"
        else

          dnl Or DEC's cxx compiler?
          is_dec_cxx="`($CXX -V 2>&1) | grep 'Compaq C++'`"
          if test "x$is_dec_cxx" != "x" ; then
            AC_MSG_RESULT(C++ compiler is Compaq-cxx)
            GXX_BRAND=Compaq
            GXX_VERSION=compaq_cxx
            GXX_VERSION_DETAILED="$GXX_VERSION"
          else

      	    dnl Sun Workshop/Studio?
            is_sun_cc_1="`($CXX -V 2>&1) | grep 'Sun WorkShop'`"
            is_sun_cc_2="`($CXX -V 2>&1) | grep 'Sun C++'`"
            if test "x$is_sun_cc_1$is_sun_cc_2" != "x" ; then
              AC_MSG_RESULT(C++ compiler is Sun Workshop compiler)
              GXX_BRAND=SunWorkshop
              GXX_VERSION=sun_workshop
              GXX_VERSION_DETAILED="$GXX_VERSION"
            else

  	      dnl Sun Forte?
              is_sun_forte_cc="`($CXX -V 2>&1) | grep 'Forte'`"
              if test "x$is_sun_forte_cc" != "x" ; then
                AC_MSG_RESULT(C++ compiler is Sun Forte compiler)
                GXX_BRAND=SunForte
                GXX_VERSION=sun_forte
                GXX_VERSION_DETAILED="$GXX_VERSION"
              else

  	        dnl Portland Group C++?
  	        is_pgcc="`($CXX -V 2>&1) | grep 'Portland Group'`"
  	        if test "x$is_pgcc" != "x" ; then
		  GXX_VERSION_STRING=`($CXX -V 2>&1) | grep "^pgCC"`
                  full_version=`echo "$GXX_VERSION_STRING" | perl -pi -e 's/.*pgCC\s+(\S+).*/\1/g;'`
		  GXX_BRAND=PortlandGroup
    		  GXX_VERSION=pgCC`echo $full_version | perl -pi -e 's/(\d\.\d).*/\1/g;'`
    		  GXX_VERSION_DETAILED=pgCC"$full_version"
  	          AC_MSG_RESULT(C++ compiler is Portland Group C++ $full_version)
                else

  	          dnl HP aCC?
  	          is_aCC="`($CXX -V 2>&1) | grep 'aCC'`"
  	          if test "x$is_aCC" != "x" ; then
  	            AC_MSG_RESULT(C++ compiler is HP aCC)
  	            GXX_BRAND=HP
                    GXX_VERSION=hp_aCC
  	            GXX_VERSION_DETAILED="$GXX_VERSION"
                  else

  	            dnl Borland C++
  	            is_bcc="`($CXX -h 2>&1) | grep 'Borland'`"
  	            if test "x$is_bcc" != "x" ; then
  	              AC_MSG_RESULT(C++ compiler is Borland C++)
  	              GXX_BRAND=Borland
		      GXX_VERSION=borland_bcc
  	              GXX_VERSION_DETAILED="$GXX_VERSION"
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
  	                GXX_BRAND=KAI
			GXX_VERSION=kai_cc
   	                GXX_VERSION_DETAILED="$GXX_VERSION"
                      else

                        dnl Maybe PathScale's compiler?
                        is_pathscale="`($CXX -v 2>&1) | grep PathScale`"
                        if test "x$is_pathscale" != "x" ; then
                          AC_MSG_RESULT(C++ compiler is PathScale C++)
  	                  GXX_BRAND=PathScale
			  GXX_VERSION=pathscale_cc
   	                  GXX_VERSION_DETAILED="$GXX_VERSION"
                        else

			  is_clang="`($CXX --version 2>&1) | grep clang`"
			  if test "x$is_clang" != x ; then
                            AC_MSG_RESULT(C++ compiler is clang)
  	                    GXX_BRAND=clang
			    GXX_VERSION=clang
   	                    GXX_VERSION_DETAILED="$GXX_VERSION"
			  else

                            dnl  Aw, nothing suitable found...
                            AC_MSG_RESULT(Unrecognized C++ compiler -- Try to go ahead and get help from dealii@dealii.org)
                            GXX_BRAND=Unknown
			    GXX_VERSION=unknown_cc
   	                    GXX_VERSION_DETAILED="$GXX_VERSION"
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
    fi
  fi
])



dnl -------------------------------------------------------------
dnl See whether the compiler we have has MPI built in (e.g. if it
dnl is actually mpiCC, etc)
dnl
dnl Usage: DEAL_II_DETERMINE_IF_SUPPORTS_MPI
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_DETERMINE_IF_SUPPORTS_MPI, dnl
[
  AC_MSG_CHECKING(if the compiler is built for MPI)
  AC_TRY_LINK(
        [
          #include <mpi.h>
        ],
        [
	  MPI_Init (0,0);
	],
        [
          AC_MSG_RESULT(yes)
	  AC_DEFINE(DEAL_II_COMPILER_SUPPORTS_MPI, 1,
                    [Defined if the compiler supports including <mpi.h>])

	  DEAL_II_COMPILER_SUPPORTS_MPI=1
          DEAL_II_USE_MPI=yes
        ],
        [
          AC_MSG_RESULT(no)
        ])
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
  dnl If CXXFLAGSG or CXXFLAGSO are not set, then initialize them with
  dnl CXXFLAGS
  if test "x$CXXFLAGSG" = "x" ; then
    CXXFLAGSG="$CXXFLAGS" ;
  fi

  if test "x$CXXFLAGSO" = "x" ; then
    CXXFLAGSO="$CXXFLAGS" ;
  fi

  dnl In no case do we want to induce BOOST to use deprecated header files
  CXXFLAGSG="$CXXFLAGSG -DBOOST_NO_HASH -DBOOST_NO_SLIST"
  CXXFLAGSO="$CXXFLAGSO -DBOOST_NO_HASH -DBOOST_NO_SLIST"

  dnl First the flags for gcc compilers
  if test "$GXX" = yes ; then
    CXXFLAGSO="$CXXFLAGSO -O2 -funroll-loops -funroll-all-loops -fstrict-aliasing -Wuninitialized -felide-constructors -ftemplate-depth-128"
    CXXFLAGSG="$CXXFLAGSG -DDEBUG -Wall -W -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Wswitch -ftemplate-depth-128"

    dnl gcc 4.4 has an interesting problem in that it doesn't
    dnl care for one of BOOST signals2's header files and produces
    dnl dozens of pages of error messages of the form
    dnl   warning: invoking macro BOOST_PP_CAT argument 1: empty macro arguments are undefined in ISO C90 and ISO C++98
    dnl This can be avoided by not using -pedantic for this compiler.
    dnl For all other versions, we use this flag, however.
    if test $GXX_VERSION != gcc4.4 ; then
      CXXFLAGS="$CXXFLAGS -pedantic"
    fi

    dnl BOOST uses long long, so don't warn about this
    CXXFLAGSG="$CXXFLAGSG -Wno-long-long"

    dnl See whether the gcc we use already has a flag for C++1x features.
    OLD_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS=-std=c++0x

    AC_MSG_CHECKING(whether compiler has a flag to support C++1x)
    AC_TRY_COMPILE([], [;],
       [
         AC_MSG_RESULT(yes)
         test_cxx1x=yes
       ],
       [
         AC_MSG_RESULT(no)
         test_cxx1x=no
       ])
    CXXFLAGS="${OLD_CXXFLAGS}"

    if test "x$test_cxx1x" = "xyes" ; then
      DEAL_II_CHECK_CXX1X_COMPONENTS("-std=c++0x")
    fi

    dnl On some gcc 4.3 snapshots, a 'const' qualifier on a return type triggers a
    dnl warning. This is unfortunate, since we happen to stumble on this
    dnl in some of our template trickery with iterator classes. If necessary,
    dnl do not use the relevant warning flag
    CXXFLAGS="-Wreturn-type -Werror"
    AC_MSG_CHECKING(whether qualifiers in return types lead to a warning)
    AC_TRY_COMPILE(
          [
            const double foo() { return 1.; }
          ],
          [;],
          [
            AC_MSG_RESULT(no)
          ],
          [
            AC_MSG_RESULT(yes)
            CXXFLAGSG="$CXXFLAGSG -Wno-return-type"
          ])

    dnl Set PIC flags. On some systems, -fpic/PIC is implied, so don't set
    dnl anything to avoid a warning. on AIX make sure we always pass -lpthread
    dnl because this seems to be somehow required to make things work. Likewise
    dnl DEC OSF and linux on x86_64.
    case "$target" in
      *aix* )
	CXXFLAGSPIC=
	LDFLAGSPIC=
        LDFLAGS="$LDFLAGS -lpthread"
	;;

      *dec-osf* )
	CXXFLAGSPIC="-fPIC"
	LDFLAGSPIC="-fPIC"
        LDFLAGS="$LDFLAGS -lpthread"
	;;

      *x86_64*)
	CXXFLAGSPIC="-fPIC"
	LDFLAGSPIC="-fPIC"
	LDFLAGS="$LDFLAGS -lpthread"
	;;

      *cygwin* )
        dnl On Cygwin, when using shared libraries, there might occur
        dnl difficulties when linking libraries for several dimensions,
        dnl as some symbols are defined in all of them. This leads to a
        dnl linker error. We force the linker to ignore multiple symbols,
        dnl but of course this might lead to strange program behaviour if
        dnl you accidentally defined one symbol multiple times...
        dnl (added 2005/07/13, Ralf B. Schulz)
	dnl (modified 2005/12/20, Ralf B. Schulz)
        CXXFLAGSPIC=
	LDFLAGS="$LDFLAGS -Xlinker --allow-multiple-definition"
        SHLIBFLAGS="$SHLIBFLAGS -Xlinker --allow-multiple-definition"
        ;;

      *)
	CXXFLAGSPIC="-fPIC"
	LDFLAGSPIC="-fPIC"
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

    if test "x$DEAL_II_USE_MPI" = "xyes" ; then
      AC_MSG_CHECKING(whether MPI headers have unused parameters)
      CXXFLAGS="-Wunused-parameter -Werror"
      AC_TRY_COMPILE(
        [
          #include <mpi.h>
        ],
        [;],
        [
  	  AC_MSG_RESULT(no)
        ],
        [
          AC_MSG_RESULT(yes)
          CXXFLAGSG="$CXXFLAGSG -Wno-unused-parameter"
          CXXFLAGSO="$CXXFLAGSO -Wno-unused-parameter"
        ])
    fi

    dnl Some system specific things
    case "$target" in
      dnl Use -Wno-long-long on Apple Darwin to avoid some unnecessary
      dnl warnings. However, newer gccs on that platform do not have
      dnl this flag any more, so check whether we can indeed do this
      dnl
      dnl Also, the TBB tries to determine whether the system is
      dnl 64-bit enabled and if so, it builds the object files in 64bit.
      dnl In order to be able to link with these files, we then have to
      dnl link with -m64 as well
      *apple-darwin*)
        OLD_CXXFLAGS="$CXXFLAGS"
        CXXFLAGS=-Wno-long-double

        AC_MSG_CHECKING(whether we can use -Wno-long-double)
        AC_TRY_COMPILE([], [;],
          [
            AC_MSG_RESULT(yes)
	    CXXFLAGSG="$CXXFLAGSG -Wno-long-double"
	    CXXFLAGSO="$CXXFLAGSO -Wno-long-double"
          ],
          [
            AC_MSG_RESULT(no)
          ])

        if test "`/usr/sbin/sysctl -n hw.optional.x86_64`" = "1" ; then
          CXXFLAGS="$CXXFLAGS -m64"
          CXXFLAGSG="$CXXFLAGSG -m64"
          CXXFLAGSO="$CXXFLAGSO -m64"
          LDFLAGS="$LDFLAGS -m64"
        fi

        CXXFLAGS="${OLD_CXXFLAGS}"
        ;;

      dnl On DEC OSF, including both stdio.h and unistd.h causes a warning
      dnl from the preprocessor that cuserid is redefined as a preprocessor
      dnl variable. Suppress this if necessary by switching off warnings
      dnl from the preprocessor
      *dec-osf*)
        AC_MSG_CHECKING(for preprocessor warning with cuserid)
	CXXFLAGS="$CXXFLAGSG -Werror"
	AC_TRY_COMPILE(
          [
#	    include <stdio.h>
#	    include <unistd.h>
          ],
          [;],
          [
            AC_MSG_RESULT(no)
          ],
          [
            AC_MSG_RESULT(yes)
            CXXFLAGSG="$CXXFLAGSG -Wp,-w"
            CXXFLAGSO="$CXXFLAGSO -Wp,-w"
          ])
        ;;
    esac

    SHLIBLD="$CXX"

  else
    dnl Non-gcc compilers. By default, use the C++ compiler also for linking
    dnl shared libraries. If some compiler cannot do that and needs something
    dnl different, then this must be specified in the respective section
    dnl below, overriding this define:
    SHLIBLD="$CXX"

    case "$GXX_VERSION" in
      ibm_xlc)
          CXXFLAGSG="$CXXFLAGSG -DDEBUG -check=bounds -info=all -qrtti=all -qsuppress=1540-2907 -qsuppress=1540-2909"
          CXXFLAGSO="$CXXFLAGSO -O2 -w -qansialias -qrtti=all -qsuppress=1540-2907 -qsuppress=1540-2909"
          CXXFLAGSPIC="-qpic"
          LDFLAGSPIC="-qpic"
          ;;

      MIPSpro*)
          dnl Disable some compiler warnings, as they trigger some warnings in
          dnl system and compiler include files
          dnl cc-1429 CC: WARNING File = /usr/include/internal/stdlib_core.h, Line = 128
          dnl The type "long long" is nonstandard.
          dnl cc-1066 CC: WARNING File = /usr/include/CC/stl_ctype.h, Line = 28
          dnl The indicated enumeration value is out of "int" range.
          dnl cc-1485 CC: WARNING File = /usr/include/CC/iomanip, Line = 122
          dnl This form for taking the address of a member function is nonstandard.
          CXXFLAGSG="$CXXFLAGSG -DDEBUG -D__sgi__ -no_auto_include -ansiW -woff 1429,1066,1485"
          dnl Disable some compiler warnings, that warn about variables
          dnl which are used in Assert templates but not in optimized mode
          dnl cc-1174 CC: full_matrix.templates.h, Line = 1461
          dnl The variable "typical_diagonal_element" was declared but never referenced.
          dnl cc-1552 CC: WARNING File = source/data_out_base.cc, Line = 3493
          dnl The variable "ierr" is set but never used.
          CXXFLAGSO="$CXXFLAGSO -D__sgi__ -O2 -no_auto_include -woff 1174,1552"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-KPIC"
          dnl Avoid output of prelinker (-quiet_prelink)
          dnl Disable compiler warning:
          dnl ld32: WARNING 131: Multiply defined weak symbol:
          dnl   (_M_acquire_lock__Q2_3std15_STL_mutex_lockGv) in auto_derivative_function.g.o
          dnl   and convergence_table.g.o (2nd definition ignored)
          LDFLAGS="$LDFLAGS -quiet_prelink -woff 131"
          dnl
          dnl Always link with math library: The -lm option must be at the end of the
          dnl linker command, therefore it cannot be included into LDFLAGS
          DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lm)
          ;;

      clang*)
          dnl Like many other compilers, clang produces warnings for array
	  dnl accesses out of bounds, even if they are in code that's dead
	  dnl for this dimension. Suppress this.
	  dnl
	  dnl There are a number of other warnings we get that can't easily
	  dnl be worked around and that are definitely not useful. Suppress
	  dnl those too.
          CXXFLAGSG="$CXXFLAGS -DDEBUG -g -Wall -Wno-array-bounds -Wno-parentheses -Wno-delete-non-virtual-dtor -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-variable"
          CXXFLAGSO="$CXXFLAGS -O2 -Wno-array-bounds -Wno-parentheses -Wno-delete-non-virtual-dtor -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-variable"
	  CXXFLAGSPIC="-fPIC"
	  LDFLAGSPIC="-fPIC"
	  ;;

      intel_icc*)
          dnl Earlier icc versions used -Kxxx for flags. Later versions use
          dnl the gcc convention -fxxx. Also, at least since icc11, the
          dnl flag -inline_debug_info has been deprecated, so don't use
	  dnl it any more
	  dnl
	  dnl Exception handling is also standard in later versions, as is rtti
    	  case "$GXX_VERSION" in
	    intel_icc5 | intel_icc6 | intel_icc7 | intel_icc8 | intel_icc9)
	      CXXFLAGSG="$CXXFLAGSG -Kc++eh -Krtti -DDEBUG -inline_debug_info"
              CXXFLAGSO="$CXXFLAGSO -Kc++eh -Krtti -O2 -unroll"
              CXXFLAGSPIC="-KPIC"
              LDFLAGSPIC="-KPIC"
              ;;

	    intel_icc*)
	      CXXFLAGSG="$CXXFLAGSG -DDEBUG"
              CXXFLAGSO="$CXXFLAGSO -O2 -unroll"
              CXXFLAGSPIC="-fPIC"
	      LDFLAGS="$LDFLAGS -lstdc++ -lpthread"
              LDFLAGSPIC="-fPIC"
              ;;
          esac

          dnl Disable some compiler warnings, as they often are wrong on
          dnl our code:
	  dnl  #11: ` unrecognized preprocessing directive" (we use
	  dnl       #warning at one place)
          dnl #175: `subscript out of range' (doesn't take into account that
          dnl       some code is only reachable for some dimensions)
          dnl #327: `NULL reference is not allowed' (this happens when we
          dnl       write "*static_cast<double*>(0)" or some such thing,
          dnl       which we do to create invalid references)
          dnl #424: `extra ";" ignored'
          dnl #525: `type "DataOutBase::DataOutBase" is an inaccessible type
          dnl       (allowed for compatibility)' (I don't understand what the
          dnl       compiler means)
	  dnl #734: `X::X(const X&), required for copy that was eliminated, is
	  dnl       inaccessible'
	  dnl       (valid, but annoying and sometimes hard to work around)
	  dnl #858: `type qualifier on return type is meaningless'
	  dnl       (on conversion operators to types that are already const)
          CXXFLAGSG="$CXXFLAGSG -w1 -wd175 -wd525 -wd327 -wd424 -wd11 -wd734 -wd858"
          CXXFLAGSO="$CXXFLAGSO -w0 -wd424 -wd11"

          dnl To reduce output, use -opt_report_levelmin where possible,
          dnl i.e. post icc5. from icc10 onwards, this flag is called
	  dnl -opt-report, and -vec-report controls output of the
          dnl autovectorizer (to make things simpler, one of the two options
          dnl wants a space between option and level, whereas the other does
          dnl not)
	  dnl
	  dnl Since only the x86 and x86_64 compilers can vectorize, this
	  dnl flag needs to be suppressed on ia64 (itanium)
    	  case "$GXX_VERSION" in
	    intel_icc5)
              ;;
	    intel_icc6 | intel_icc7 | intel_icc8 | intel_icc9)
              CXXFLAGSO="$CXXFLAGSO -opt_report_levelmin"
              ;;
	    *)
	      case "$target" in
		*ia64*)
              	    CXXFLAGSO="$CXXFLAGSO -opt-report 0"
		    ;;
		*)
              	    CXXFLAGSO="$CXXFLAGSO -opt-report 0 -vec-report0"
		    ;;
	      esac
	      ;;
          esac

	  dnl Some versions of icc on some platforms issue a lot of warnings
	  dnl about the unreliability of floating point comparisons. Check
	  dnl whether we can switch that off
	  DEAL_II_ICC_WD_1572

          dnl We would really like to use  -ansi -Xc, since that
	  dnl is _very_ picky about standard C++, and is thus very efficient
          dnl in detecting slight standard violations, but these flags are
          dnl also very efficient in crashing the compiler (it generates a
          dnl segfault), at least with versions prior to 7.0. So only
          dnl use these flags with versions we know are safe
          dnl
          dnl Second thing: icc7 allows using alias information for
          dnl optimization. Use this.
          if test "x$GXX_VERSION" = "xintel_icc7"; then
            CXXFLAGSG="$CXXFLAGSG -Xc -ansi"
            CXXFLAGSO="$CXXFLAGSO -ansi_alias"
          dnl For icc8:
	  dnl    avoid the annoying `LOOP WAS VECTORIZED' remarks
          dnl    use -vec_report0 for reducing output
          else if test "x$GXX_VERSION" = "xintel_icc8" ; then
            CXXFLAGSO="$CXXFLAGSO -ansi_alias -vec_report0"
          fi fi

	  dnl If we are on an x86 platform, add -tpp6 to optimization
	  dnl flags (for version <10), or the equivalent for later
	  dnl processors
	  case "$target" in
            *x86_64*)
	        LDFLAGS="$LDFLAGS -lpthread"
	        ;;

	    *86*)
    	        case "$GXX_VERSION" in
	          intel_icc5)
                    ;;
	          intel_icc6 | intel_icc7 | intel_icc8 | intel_icc9)
                    CXXFLAGSO="$CXXFLAGSO -tpp6"
                    ;;
	          *)
                    CXXFLAGSO="$CXXFLAGSO -mcpu=pentium4"
	            ;;
                esac
		;;
	  esac
          ;;

      compaq_cxx)
          dnl Disable some warning messages:
	  dnl  #11: ` unrecognized preprocessing directive" (we use
	  dnl       #warning at one place)
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
	  dnl #450:`the type "long long" is nonstandard'
	  dnl       BOOST uses long long, unfortunately
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

          CXXFLAGSG="$CXXFLAGSG -model ansi -std strict_ansi -w1 -msg_display_number -timplicit_local -DDEBUG"
          CXXFLAGSO="$CXXFLAGSO -model ansi -std strict_ansi -w2 -msg_display_number -timplicit_local -fast"

          for i in 11 175 236 237 381 487 1136 1156 111 1182 265 450 ; do
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
          LDFLAGSPIC=""
          ;;

      sun_workshop | sun_forte)

          CXXFLAGSG="$CXXFLAGSG -DDEBUG -w"
          CXXFLAGSO="$CXXFLAGSO -w"
          CXXFLAGSPIC="-KPIC"
          LDFLAGSPIC="-G"

	  dnl See if the flag -library=stlport4 is available, and if so use it
	  CXXFLAGS="$CXXFLAGSG -library=stlport4"
	  AC_MSG_CHECKING(whether -library=stlport4 works)
	  AC_TRY_COMPILE(
            [
	    ],
            [
	      ;
            ],
            [
              AC_MSG_RESULT(yes)
              CXXFLAGSG="$CXXFLAGSG -library=stlport4"
              CXXFLAGSO="$CXXFLAGSO -library=stlport4"
            ],
            [
              AC_MSG_RESULT(no)
            ])
 	  ;;

      pgCC*)
	  dnl Suppress warnings:
	  dnl  #68: "integer conversion resulted in a change of sign". This
	  dnl       is what we get every time we use
	  dnl       numbers::invalid_unsigned_int
          dnl #111: "Statement unreachable": we use return statements
	  dnl       occasionally after case-switches where you cannot
 	  dnl       fall though, but other compilers sometimes complain
	  dnl       that the function might not return with a value, if
          dnl       it can't figure out that the function always uses
          dnl       one case. Also: a return statement after a failing
          dnl       assertion
	  dnl #128: a similar case -- code not reachable because there's
	  dnl       a return that's active for a particular space dimension
	  dnl #155: no va_start() seen -- happens alwas in lines 138 and 891
	  dnl       and seems to be spurious
          dnl #177: "function declared but not used": might happen with
          dnl       templates and conditional compilation
 	  dnl #175: "out-of-bounds array indices": the same reason as
	  dnl       for Compaq cxx
	  dnl #185: "dynamic initialization in unreachable code". similar to
	  dnl       the case #128 above
	  dnl #236: "controlling expression is constant". this triggers
	  dnl       somewhere in BOOST with BOOST_ASSERT. I have no idea
	  dnl       what happens here
 	  dnl #284: "NULL references not allowed"
	  CXXFLAGSG="$CXXFLAGSG -DDEBUG -g --display_error_number --diag_suppress 68 --diag_suppress 111 --diag_suppress 128 --diag_suppress 155 --diag_suppress 177 --diag_suppress 175 --diag_suppress 185 --diag_suppress 236 --diag_suppress 284"
          CXXFLAGSO="$CXXFLAGSO -fast -O2 --display_error_number --diag_suppress 68 --diag_suppress 111 --diag_suppress 128 --diag_suppress 155 --diag_suppress 177 --diag_suppress 175 --diag_suppress 185 --diag_suppress 236 --diag_suppress 284"
          CXXFLAGSPIC="-Kpic"
          ;;

      kai_cc)
          CXXFLAGSG="$CXXFLAGSG --strict -D__KAI_STRICT --max_pending_instantiations 32 --display_error_number -g +K0 --no_implicit_typename"
          CXXFLAGSO="$CXXFLAGSO +K3 -O2 --abstract_float --abstract_pointer -w --display_error_number --max_pending_instantiations 32 --display_error_number"
          CXXFLAGSPIC="-fPIC"
	  ;;

      hp_aCC)
	  dnl ??? disable warning 655 (about all-inlined functions) which
	  dnl triggers for each and every of our DeclExceptionX calls ???
          CXXFLAGSG="$CXXFLAGSG -g1 -AA +p"
          CXXFLAGSO="$CXXFLAGSO -z +O2 -AA"
          CXXFLAGSPIC="+Z"
	  # for linking shared libs, -b is also necessary...
          ;;

      borland_bcc)
          CXXFLAGSG="$CXXFLAGSG -q -DDEBUG -w -w-use -w-amp -w-prc"
          CXXFLAGSO="$CXXFLAGSO -q -O2"
          CXXFLAGSPIC=""
          LDFLAGSPIC=""
          AC_MSG_ERROR(Attention! deal.II is not known to work with Borland C++!
		If you intend to port it to Borland C++, please remove this message from aclocal.m4 and call autoconf and configure. If you do not understand this, you will NOT want to do it!)
          ;;

      pathscale_cc)
	  CXXFLAGSG="$CXXFLAGSG -DDEBUG -g"
	  CXXFLAGSO="$CXXFLAGSO -O3"
          ;;

      *)
	  CXXFLAGSG="$CXXFLAGSG -DDEBUG"
	  CXXFLAGSO="$CXXFLAGSO -O2"
          AC_MSG_RESULT(Unknown C++ compiler - using generic options)
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
    case "$target" in
       dnl On Alpha, use the special treatment
       alpha*-osf*)
           AC_MSG_CHECKING(whether -ggdb works for long symbols)
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
               AC_MSG_RESULT(no -- using -gstabs instead)
             ])
           ;;

       dnl For all other systems test whether -ggdb works at all, and if
       dnl not fall back on -g instead. This test is mainly used to
       dnl accomodate for a failure on Mac OS X, where -ggdb leads to
       dnl information the assembler does not understand (see mailing
       dnl list thread on this from mid-October 2002)
       *)
           AC_MSG_CHECKING(whether -ggdb works)
           CXXFLAGS="-ggdb $CXXFLAGSG"
	   AC_TRY_COMPILE(
             [],
             [ ; ],
             [
               CXXFLAGSG="-ggdb $CXXFLAGSG"
               AC_MSG_RESULT(yes)
             ],
             [
               CXXFLAGSG="-g $CXXFLAGSG"
               AC_MSG_RESULT(no -- using -g instead)
             ])
           ;;
    esac

  else
    dnl Non-gcc compilers use -g instead of -ggdb, except for Borland C++
    dnl which wants something entirely different
    case "$GXX_VERSION" in
      borland_bcc)
	CXXFLAGSG="-v -y $CXXFLAGSG"
	;;

      *)
    	CXXFLAGSG="-g $CXXFLAGSG"
	;;
    esac
  fi
])



dnl -------------------------------------------------------------
dnl Given the command line flag specified as argument to this macro,
dnl test whether all components that we need from the C++1X
dnl standard are actually available. If so, add the flag to
dnl CXXFLAGS.g and CXXFLAGS.o, and set a flag in config.h
dnl
dnl Usage: DEAL_II_CHECK_CXX1X_COMPONENTS(cxxflag)
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CXX1X_COMPONENTS, dnl
[
  OLD_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$1"

  all_cxx1x_classes_available=yes

  AC_MSG_CHECKING(for std::array)
  AC_TRY_COMPILE(
       [#include <array>],
       [ std::array<int,3> p; p[0];],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  AC_MSG_CHECKING(for std::condition_variable)
  AC_TRY_COMPILE(
       [#include <condition_variable> ],
       [ std::condition_variable c; c.notify_all()],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  AC_MSG_CHECKING(for std::function and std::bind)
  AC_TRY_COMPILE(
       [#include <functional>
        void f(int, double); ],
       [ std::function<void (int)>
            g = std::bind (f, std::placeholders::_1, 1.1) ;],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  dnl Make sure we don't run into GCC bug 35569
  AC_MSG_CHECKING(for std::bind works with rvalues)
  AC_TRY_COMPILE(
       [#include <functional>
        void f(int); ],
       [ using namespace std;
         using namespace std::placeholders;
         bind(multiplies<int>(),4,_1)(5); ;],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  AC_MSG_CHECKING(for std::shared_ptr)
  AC_TRY_COMPILE(
       [#include <memory>],
       [ std::shared_ptr<int> p(new int(3))],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  AC_MSG_CHECKING(for std::thread)
  AC_TRY_COMPILE(
       [#include <thread>
        void f(int); ],
       [ std::thread t(f,1); t.join();],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  dnl On some systems with gcc 4.5.0, we can compile the code
  dnl above but it will throw an exception when run. So test
  dnl that as well. The test will only be successful if we have
  dnl libpthread available, so link it in for this test. If
  dnl multithreading is requested, it will be added to CXXFLAGS
  dnl later on so there is no need to do this here.
  CXXFLAGS="$1 -lpthread"
  AC_MSG_CHECKING(whether std::thread actually works)
  AC_TRY_RUN(
       [#include <thread>
        void f(int) {}
        int main() { std::thread t(f,1); t.join(); } ],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )
  CXXFLAGS="$1"

  AC_MSG_CHECKING(for std::mutex)
  AC_TRY_COMPILE(
       [#include <mutex> ],
       [ std::mutex m; m.lock();],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  AC_MSG_CHECKING(for std::tuple)
  AC_TRY_COMPILE(
       [#include <tuple>],
       [ std::tuple<int,double,char> p(1,1.1,'a')],
       [ AC_MSG_RESULT(yes) ],
       [ AC_MSG_RESULT(no); all_cxx1x_classes_available=no ]
       )

  CXXFLAGS="${OLD_CXXFLAGS}"

  dnl If the above classes and operations are all defined then we can
  dnl use the corresponding flag in CXXFLAGS* to switch on support and
  dnl correspondingly use the C++ classes instead of the BOOST classes
  AC_MSG_CHECKING(whether C++1x support is complete enough)
  if test "x$all_cxx1x_classes_available" = "xyes" ; then
    AC_MSG_RESULT(yes)

    CXXFLAGSG="$CXXFLAGSG $1"
    CXXFLAGSO="$CXXFLAGSO $1"

    AC_DEFINE(DEAL_II_CAN_USE_CXX1X, 1,
              [Defined if the compiler we use supports the upcoming C++1x standard.])


    dnl Also test for a couple C++1x things that we don't use in the library
    dnl but that users may want to use in their applications and that we
    dnl might want to test in the testsuite
    OLD_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGSG"

    extra_cxx1x_features_available=yes

    AC_MSG_CHECKING(for auto typed variables)
    AC_TRY_COMPILE(
         [#include <vector>],
         [
           std::vector<int> v;
	   auto i = v.begin();
	   *i;
         ],
         [ AC_MSG_RESULT(yes) ],
         [ AC_MSG_RESULT(no); extra_cxx1x_features_available=no ]
         )

    AC_MSG_CHECKING(for range-based for)
    AC_TRY_COMPILE(
         [#include <vector>],
         [
           std::vector<int> v;
	   for (std::vector<int>::iterator i : v)
	     *i;
         ],
         [ AC_MSG_RESULT(yes) ],
         [ AC_MSG_RESULT(no); extra_cxx1x_features_available=no ]
         )

    CXXFLAGS="${OLD_CXXFLAGS}"

    dnl Right now, do nothing with this information since we don't
    dnl have any compiler that gets all this right (gcc 4.5 included)
    dnl and that we could use to test these features.
  else
    AC_MSG_RESULT(no)
  fi
])


dnl -------------------------------------------------------------
dnl Determine the C compiler in use. Return the name and possibly
dnl version of this compiler in CC_VERSION. This function is almost
dnl identifical to DEAL_II_DETERMINE_CXX_BRAND and therefore lacks
dnl a lot of the comments found there to keep it short
dnl
dnl Usage: DEAL_II_DETERMINE_CC_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_DETERMINE_CC_BRAND, dnl
[
  if test "$GCC" = "yes" ; then
    CC_VERSION_STRING=`($CC -v 2>&1) | grep "gcc version"`
    if test "x$CC_VERSION_STRING" = "x" ; then
      GCC=no
    fi
  fi

  if test "$GCC" = yes ; then
    dnl find out the right version
    CC_VERSION_STRING=`($CC -v 2>&1) | grep "gcc version"`

    full_version=`echo "$CC_VERSION_STRING" | perl -pi -e 's/.*version (\d\.\d\.\d).*/\1/g;'`
    CC_VERSION=gcc`echo $full_version | perl -pi -e 's/(\d\.\d).*/\1/g;'`

    AC_MSG_RESULT(C compiler is $CC_VERSION)
  else
    dnl Check other (non-gcc) compilers

    dnl Check for IBM xlC. For some reasons, depending on some environment
    dnl variables, moon position, and other reasons unknown to me, the
    dnl compiler displays different names in the first line of output, so
    dnl check various possibilities
    is_ibm_xlc="`($CC 2>&1) | egrep 'VisualAge C|C Set ++|C for AIX Compiler'`"
    if test "x$is_ibm_xlc" != "x"  ; then
      dnl Ah, this is IBM's C compiler. Unfortunately, we don't presently
      dnl know how to check the version number, so assume that is sufficiently
      dnl high...
      AC_MSG_RESULT(C compiler is IBM xlC)
      CC_VERSION=ibm_xlc
    else

      dnl Check whether we are dealing with the MIPSpro C compiler
      mips_pro="`($CC -version 2>&1) | grep MIPSpro`"
      if test "x$mips_pro" != "x" ; then
        case "$mips_pro" in
          *7.0* | *7.1* | *7.2* | *7.3*)
            dnl MIPSpro 7.3 does not support standard C++, therefore it is not
            dnl able to compile deal.II. Previous compiler versions neither.
            AC_MSG_RESULT(C compiler is $mips_pro)
            AC_MSG_ERROR(This compiler is not supported)
            CC_VERSION=MIPSpro7.3
            ;;
          *7.4)
            AC_MSG_RESULT(C compiler is MIPSpro compiler 7.4)
            AC_MSG_ERROR(This compiler is not supported. Use MIPSPro compiler 7.4x)
            CC_VERSION=MIPSpro7.4
            ;;
          *7.41* | *7.42* | *7.43* | *7.44*)
            AC_MSG_RESULT(C compiler is MIPSpro compiler 7.4x)
            CC_VERSION=MIPSpro7.4x
            ;;
          *"7.5"*)
            AC_MSG_RESULT(C compiler is MIPSpro compiler 7.5)
            CC_VERSION=MIPSpro7.5
            ;;
          *)
            AC_MSG_RESULT(C compiler is unknown version but accepted MIPSpro compiler version)
            CC_VERSION=MIPSpro-other
            ;;
        esac
      else

        dnl Intel's ICC C compiler? On Linux, it uses -V, on Windows
	dnl it is -help
	dnl
	dnl Annoyingly, ecc6.0 prints its version number on a separate
	dnl line (the previous one ends with the string "applications"),
	dnl so join this one to the previous one with a little bit of
	dnl perl.
        is_intel_icc1="`($CC -V 2>&1) | grep 'Intel'`"
        is_intel_icc2="`($CC -help 2>&1) | grep 'Intel'`"
        is_intel_ecc="`($CC -V 2>&1) | perl -pi -e 's/applications\n/\1/g;' | grep 'Intel(R) C++ Itanium(TM) Compiler'`"
	is_intel_icc="$is_intel_icc1$is_intel_icc2$is_intel_ecc"

	dnl When calling the Portland Group compiler, it also
	dnl outputs the string 'Intel' in its help text, so make
	dnl sure we're not confused
	is_pgi="`($CXX -V 2>&1) | grep 'Portland'`"

        if test "x$is_intel_icc" != "x" -a "x$is_pgi" = "x"; then
	  version_string="`($CC -V 2>&1) | grep 'Version'` `($CC -help 2>&1) | grep 'Version'`"
	  version5="`echo $version_string | grep 'Version 5'`"
	  version6="`echo $version_string | grep 'Version 6'`"
	  version7="`echo $version_string | grep 'Version 7'`"
	  version8="`echo $version_string | grep 'Version 8'`"
	  version9="`echo $version_string | grep 'Version 9'`"
	  version10="`echo $version_string | grep 'Version 10'`"
          if test "x$version5" != "x" ; then
            AC_MSG_RESULT(C compiler is icc-5)
            CC_VERSION=intel_icc5
          else if test "x$version6" != "x" ; then
            AC_MSG_RESULT(C compiler is icc-6)
            CC_VERSION=intel_icc6
          else if test "x$version7" != "x" ; then
            AC_MSG_RESULT(C compiler is icc-7)
            CC_VERSION=intel_icc7
          else if test "x$version8" != "x" ; then
            AC_MSG_RESULT(C compiler is icc-8)
            CC_VERSION=intel_icc8
          else if test "x$version9" != "x" ; then
            AC_MSG_RESULT(C compiler is icc-9)
            CC_VERSION=intel_icc9
          else if test "x$version10" != "x" ; then
            AC_MSG_RESULT(C compiler is icc-10)
            CC_VERSION=intel_icc10
          else
            AC_MSG_RESULT(C compiler is icc)
            CC_VERSION=intel_icc
          fi fi fi fi fi fi
        else

          dnl Or DEC's cxx compiler?
          is_dec_cxx="`($CC -V 2>&1) | grep 'Compaq C'`"
          if test "x$is_dec_cxx" != "x" ; then
            AC_MSG_RESULT(C compiler is Compaq cxx)
            CC_VERSION=compaq_cxx
          else

      	    dnl Sun Workshop?
            is_sun_cc="`($CC -V 2>&1) | grep 'Sun WorkShop'`"
            if test "x$is_sun_cc" != "x" ; then
              AC_MSG_RESULT(C compiler is Sun Workshop compiler)
              CC_VERSION=sun_workshop
            else

  	      dnl Sun Forte?
              is_sun_forte_cc="`($CC -V 2>&1) | grep 'Forte'`"
              is_sun_cc="`($CC -V 2>&1) | grep 'Sun C'`"
              if test "x$is_sun_forte_cc$is_sun_cc" != "x" ; then
                AC_MSG_RESULT(C compiler is Sun C compiler)
                CC_VERSION=sun_cc
              else

  	        dnl Portland Group C?
  	        is_pgcc="`($CC -V 2>&1) | grep 'Portland Group'`"
  	        if test "x$is_pgcc" != "x" ; then
		  CC_VERSION_STRING=`($CC -V 2>&1) | grep "^pgcc"`
                  full_version=`echo "$CC_VERSION_STRING" | perl -pi -e 's/.*pgcc\s+(\S+).*/\1/g;'`
    		  CC_VERSION=pgcc`echo $full_version | perl -pi -e 's/(\d\.\d).*/\1/g;'`
  	          AC_MSG_RESULT(C compiler is Portland Group C $full_version)
  	        else

  	          dnl HP aCC?
  	          is_aCC="`($CC -V 2>&1) | grep 'aCC'`"
  	          if test "x$is_aCC" != "x" ; then
  	            AC_MSG_RESULT(C compiler is HP aCC)
  	            CC_VERSION=hp_aCC
  	          else

  	            dnl Borland C
  	            is_bcc="`($CC -h 2>&1) | grep 'Borland'`"
  	            if test "x$is_bcc" != "x" ; then
  	              AC_MSG_RESULT(C compiler is Borland C)
  	              CC_VERSION=borland_bcc
  	            else

  	              dnl KAI C? It seems as if the documented options
		      dnl -V and --version are not always supported, so give
	              dnl the whole thing a second try by looking for /KCC/
	 	      dnl somewhere in the paths that are output by -v. This
	              dnl is risky business, since this combination of
	              dnl characters might appear on other installations
                      dnl of other compilers as well, so put this test to
                      dnl the very end of the detection chain for the
                      dnl various compilers
  	              is_kai_cc="`($CC --version 2>&1) | grep 'KAI C'`"
  	              is_kai_cc="$is_kai_cc`($CC -v 2>&1) | grep /KCC/`"
  	              if test "x$is_kai_cc" != "x" ; then
  	                AC_MSG_RESULT(C compiler is KAI C)
  	                CC_VERSION=kai_cc
  	              else

                        dnl  Aw, nothing suitable found...
                        AC_MSG_RESULT([Unrecognized compiler -- still trying])
			CC_VERSION=unknown_cc
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
dnl Set C compiler flags to their default values. They will be
dnl modified according to other options in later steps of
dnl configuration
dnl
dnl CFLAGS.o  : flags for optimized mode
dnl CFLAGS.g  : flags for debug mode
dnl
dnl Usage: DEAL_II_SET_CC_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_SET_CC_FLAGS, dnl
[
  dnl First the flags for gcc compilers
  if test "$GCC" = yes ; then
    CFLAGSO="$CFLAGS -O3 -funroll-loops -funroll-all-loops -fstrict-aliasing"
    CFLAGSG="$CFLAGS -g"
    dnl Set PIC flags. On some systems, -fpic/PIC is implied, so don't set
    dnl anything to avoid a warning. on AIX make sure we always pass -lpthread
    dnl because this seems to be somehow required to make things work. Likewise
    dnl DEC OSF.
    case "$target" in
      *aix* )
	CFLAGSPIC=
	;;

      *dec-osf* )
	CFLAGSPIC="-fPIC"
	;;

      *cygwin*)
        CFLAGSPIC=
        ;;

      *apple-darwin*)
        dnl The TBB tries to determine whether the system is
        dnl 64-bit enabled and if so, it builds the object files in 64bit.
        dnl In order to be able to link with these files, we then have to
        dnl link with -m64 as well
        if test "`/usr/sbin/sysctl -n hw.optional.x86_64`" = "1" ; then
          CFLAGS="$CFLAGS -m64"
        fi
        ;;

      *)
	CFLAGSPIC="-fPIC"
	;;
    esac

  else
    dnl Non-gcc compilers. By default, use the C compiler also for linking
    dnl shared libraries. If some compiler cannot do that and needs something
    dnl different, then this must be specified in the respective section
    dnl below, overriding this define:
    SHLIBLD="$CC"

    case "$CC_VERSION" in
      ibm_xlc)
          CFLAGSO="$CFLAGS -O2"
          CFLAGSPIC="-fPIC"
	  SHLIBLD="$CXX"
          ;;

      MIPSpro*)
          CFLAGSO="$CFLAGS -O2"
          CFLAGSPIC="-KPIC"
          ;;

      intel_icc*)
          CFLAGSO="$CFLAGS -O2 -unroll"
    	  case "$CC_VERSION" in
	    intel_icc5 | intel_icc6 | intel_icc7 | intel_icc8 | intel_icc9)
                CFLAGSPIC="-KPIC"
		;;

	    intel_icc*)
		CFLAGSPIC="-fPIC"
		;;
	  esac

          dnl To reduce output, use -opt_report_levelmin where possible,
          dnl i.e. post icc5. from icc10 onwards, this flag is called
	  dnl -opt-report, and -vec-report controls output of the
          dnl autovectorizer (to make things simpler, one of the two options
          dnl wants a space between option and level, whereas the other does
          dnl not)
    	  case "$CC_VERSION" in
	    intel_icc5)
              ;;
	    intel_icc6 | intel_icc7 | intel_icc8 | intel_icc9)
              CFLAGSO="$CFLAGSO -opt_report_levelmin"
              ;;
	    *)
              CFLAGSO="$CFLAGSO -opt-report 0 -vec-report0"
	      ;;
          esac

          CFLAGSO="$CFLAGSO -ansi_alias -vec_report0"

	  dnl If we are on an x86 platform, add -tpp6 to optimization
	  dnl flags
	  case "$target" in
	    *86*)
    	        case "$CC_VERSION" in
	          intel_icc5)
                    ;;
	          intel_icc6 | intel_icc7 | intel_icc8 | intel_icc9)
                    CFLAGSO="$CFLAGSO -tpp6"
                    ;;
	          *)
                    CFLAGSO="$CFLAGSO -mcpu=pentium4"
	            ;;
                esac
		;;
	  esac

	  dnl Check whether we can switch off the annoying 1572 warning
	  dnl message about unreliable floating point comparisons
	  DEAL_II_ICC_C_WD_1572

          ;;

      sun_cc*)
          CFLAGS="$CFLAGS -g"
          CFLAGSO="$CFLAGS -fast -O2"
	  CFLAGSPIC="-fPIC"
          ;;

      *)
          AC_MSG_RESULT(Unknown C compiler - using generic options)
	  CFLAGSO="$CFLAGSO -O2"
          ;;
    esac
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

    dnl Get version string of the compiler. Some compilers, most
    dnl notably the IBM compilers have the bad habit of dumping
    dnl all of their helptexts here, so only consider the first
    dnl 10 lines. Otherwise we'll have a problem later on when
    dnl we do things like "echo $F77_VERSION_STRING | grep ..." and
    dnl the shell says that we exceeded the limit for the length of
    dnl command lines :-(
    F77_VERSION_STRING="`($F77 -v 2>&1) | head -n 10`"
    if test -n "`echo $F77_VERSION_STRING | grep \"GNU F77\"`" -o \
	    -n "`echo $F77_VERSION_STRING | grep \"gcc version\"`" ; then
      dnl Yes, this is a GNU g77 version. find out the right version
      G77_VERSION_STRING="`($F77 -v 2>&1) | grep \"gcc version\"`"

      full_version=`echo "$G77_VERSION_STRING" | perl -pi -e 's/.*version (\d\.\d\.\d).*/\1/g;'`
      F77_VERSION=gcc`echo $full_version | perl -pi -e 's/(\d\.\d).*/\1/g;'`

      AC_MSG_RESULT(F77 compiler is $CC_VERSION)

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

            F77_VERSION_STRING=`($F77 -V 2>&1)`
            is_intel_ifc="`echo $F77_VERSION_STRING | grep 'Intel(R) Fortran'`"
	    if test "x$is_intel_ifc" != "x" ; then
  	      AC_MSG_RESULT(F77 compiler is Intel Fortran)
              F77_VERSION=INTEL_F77

            else


              dnl Now, this is a hard case, we have no more clues...
              F77_VERSION="UnknownF77"
  	      AC_MSG_RESULT(F77 compiler is unknown. no flags set!)
            fi
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
    gcc*)
        F77FLAGSG="$FFLAGS -ggdb -DDEBUG -pedantic -W -Wall"
        F77FLAGSO="$FFLAGS -O2 -funroll-loops -funroll-all-loops -fstrict-aliasing"

        case "$target" in
           *cygwin* )
                F77FLAGSPIC=
                ;;

           *apple-darwin*)
		dnl Add -m64 to flags for the same TBB-related reason as
		dnl above when setting CXXFLAGS
                if test "`/usr/sbin/sysctl -n hw.optional.x86_64`" = "1" ; then
         	  F77FLAGSG="$F77FLAGSG -m64"
                  F77FLAGSO="$F77FLAGSO -m64"
                fi
		;;

           * )
	        F77FLAGSPIC="-fPIC"
                ;;
        esac

	dnl Make sure we link both possible libraries here. Shame on gfortran!
	AC_CHECK_LIB(g2c, e_wsfe, [ F77LIBS="$F77LIBS -lg2c" ])

	dnl Check whether libgfortran contains certain symbols. We used
	dnl to use _gfortran_allocate but that appears to have disappeared
	dnl at one point in time so if we can't find it check other symbols
	dnl till we find one we recognize
	AC_CHECK_LIB(gfortran, _gfortran_allocate, [found=1])
        AC_CHECK_LIB(gfortran, _gfortran_st_write_done, [found=1])

	dnl libgfortran appears to exist. Link with it.
        if test "x$found" = "x1" ; then
          F77LIBS="$F77LIBS -lgfortran"
        fi
	;;

    AIXF77)
	dnl Set flags for AIX's xlf compiler. -qextname instructs the compiler
	dnl to append an underscore to external function names, which is what
	dnl we expect when linking with FORTRAN functions
        F77FLAGSG="$FFLAGS -g -qextname"
        F77FLAGSO="$FFLAGS -O3 -w -qextname"
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

    INTEL_F77*)
            F77FLAGSG="$FFLAGS"
            F77FLAGSO="$FFLAGS -O3"
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
	AC_MSG_ERROR(No compiler options for F77 compiler "$F77_VERSION" specified: modification of aclocal.m4 necessary)
        ;;
  esac
])


dnl -------------------------------------------------------------
dnl
dnl Check whether we can use -rpath for linking. If so, we can use
dnl -rpath instead of -L in LDFLAGS below, so the dynamic linker
dnl finds all libraries.
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_RPATH,
[ OLD_LDFLAGS="$LDFLAGS"
  LDFLAGS="-Wl,-rpath,$DEAL_II_PATH/lib $LDFLAGS"
  AC_MSG_CHECKING([whether compiler understands option -Wl,-rpath])
  AC_LINK_IFELSE(
   [ AC_LANG_PROGRAM([[]],[[]])],
   [ AC_MSG_RESULT(yes)
     LD_PATH_OPTION="-Wl,-rpath,"
     DEAL_II_LD_UNDERSTANDS_RPATH=yes
   ],
   [ AC_MSG_RESULT(no)
     LDFLAGS="$OLD_LDFLAGS"
     LD_PATH_OPTION="no"
   ])
  AC_SUBST(DEAL_II_LD_UNDERSTANDS_RPATH)
])



dnl -------------------------------------------------------------
dnl
dnl Check whether we can use -soname for linking in the shared
dnl library version. On Mac OS X, -soname is called
dnl -dylib_install_name
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_LINK_SONAME,
[
  dnl First try -soname
  OLD_LDFLAGS="$LDFLAGS"
  LDFLAGS="-Wl,-soname,libbase.so.6.2.1 $LDFLAGS $LDFLAGSPIC $SHLIBFLAGS"
  AC_MSG_CHECKING([whether compiler understands option -Wl,-soname])
  AC_LINK_IFELSE(
   [ AC_LANG_PROGRAM([[]],[[]])],
   [
     AC_MSG_RESULT(yes)
     DEAL_II_LD_UNDERSTANDS_SONAME="yes"
   ],
   [
     AC_MSG_RESULT(no)
     DEAL_II_LD_UNDERSTANDS_SONAME="no"
   ]
  )
  LDFLAGS="$OLD_LDFLAGS"
  AC_SUBST(DEAL_II_LD_UNDERSTANDS_SONAME)

  dnl Now try the -dylib_install_name thing
  OLD_LDFLAGS="$LDFLAGS"
  LDFLAGS="-Wl,-dynamic,-install_name -Wl,libbase.so.6.2.1 $LDFLAGS -shared"
  AC_MSG_CHECKING([whether compiler understands option -Wl,-dynamic,-install_name])
  AC_LINK_IFELSE(
   [ AC_LANG_PROGRAM([[]],[[]])],
   [
     AC_MSG_RESULT(yes)
     DEAL_II_LD_UNDERSTANDS_DYLIB_INSTALL_NAME="yes"
   ],
   [
     AC_MSG_RESULT(no)
     DEAL_II_LD_UNDERSTANDS_DYLIB_INSTALL_NAME="no"
   ]
  )
  LDFLAGS="$OLD_LDFLAGS"
  AC_SUBST(DEAL_II_LD_UNDERSTANDS_DYLIB_INSTALL_NAME)
])



dnl -------------------------------------------------------------
dnl
dnl Check for a problem with dynamic cast on dynamic libs on
dnl Mac OS X Snow Leopard. The test is only called on Mac OS X.
dnl The test is due to Scott Miller <scott.miller@psu.edu>
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_DYNAMIC_CAST_BUG,
[
  AC_MSG_CHECKING(for dynamic_cast problem with shared libs)

  if (cd contrib/config/tests/darwin-dynamic-cast ; \
      $CXX -dynamiclib BaseClass.cpp -o libDynamicCastTestLib.dylib ; \
      $CXX -L. -lDynamicCastTestLib main.cc -o main ; \
      ./main) ; then
    AC_MSG_RESULT(no) ;
  else
    AC_MSG_RESULT(yes) ;
    AC_DEFINE(DEAL_II_HAVE_DARWIN_DYNACAST_BUG, 1,
              [Defined if the compiler has a bug with dynamic casting
               and dynamic libraries])
    CXXFLAGSG="$CXXFLAGSG -mmacosx-version-min=10.4"
    CXXFLAGSO="$CXXFLAGSO -mmacosx-version-min=10.4"
    LDFLAGS="$LDFLAGS -mmacosx-version-min=10.4"
  fi
  rm -f contrib/config/tests/darwin-dynamic-cast/libDynamicCastTestLib.dylib
  rm -f contrib/config/tests/darwin-dynamic-cast/main.o
  rm -f contrib/config/tests/darwin-dynamic-cast/main
])



dnl -------------------------------------------------------------
dnl Check whether user wants optimization for a certain type of
dnl CPU. If so, then set some flags, dependent on what he
dnl wants and the compiler. Not very many CPUS are listed here,
dnl but this is simple to expand if desired.
dnl
dnl To use this feature, use --with-cpu=something
dnl
dnl Usage: DEAL_II_CHECK_CPU_OPTIMIZATIONS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CPU_OPTIMIZATIONS, dnl
[
  AC_ARG_WITH(cpu,
              [AS_HELP_STRING([--with-cpu=cpu],
              [Optimize specifically for the given CPU type, rather than just generating code for this processor family.])],
      withcpu="$withval",
      withcpu="")
  AC_MSG_CHECKING(for CPU to optimize for)
  case "$withcpu" in
    PowerPC64)
        AC_MSG_RESULT(PowerPC64)
	case "$GXX_VERSION" in
	  gcc*)
	      dnl Tune for this processor
	      CXXFLAGSG="$CXXFLAGSG -maix64"
	      CXXFLAGSO="$CXXFLAGSO -maix64 -mpowerpc64 -mcpu=powerpc64 -mtune=powerpc64"

	      dnl On this stupid system, we get TOC overflows if we use the
	      dnl standard flags, so restrict TOC entries to the absolute minimal
	      CXXFLAGSG="$CXXFLAGSG -mminimal-toc"
	      CXXFLAGSO="$CXXFLAGSO -mminimal-toc"

	      dnl Also set 64-bit mode for f77 compiler, assuming we use IBM's xlf
	      F77FLAGSG="$F77FLAGSG -q64"
	      F77FLAGSO="$F77FLAGSO -q64"

	      dnl When generating 64-bit code, we need to pass respective flags when
	      dnl linking (also for static libs)
	      AR="$AR -X 64"
	      LDFLAGS="$LDFLAGS -maix64"

	      dnl And we must always link with pthreads
	      DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lpthread)
              ;;
        esac
	;;
    athlon* | pentium* | i386 | i486 | i586 | i686 | k6* | winchip* | opteron)
        AC_MSG_RESULT(x86 derivate ($withcpu))
	case "$GXX_VERSION" in
	  gcc*)
	      dnl Tune for this processor, but only in optimized mode
              dnl (to prevent the effects of possible compiler bugs to affect
              dnl both debug as well as optimized versions)
	      CXXFLAGSO="$CXXFLAGSO -march=$withcpu"

	      dnl Also set the mode for f77 compiler
	      F77FLAGSO="$F77FLAGSO -march=$withcpu"
          ;;
        esac
        ;;

    native)
        AC_MSG_RESULT(native processor variant)
	case "$GXX_VERSION" in
	  gcc*)
	      dnl Tune for this processor, but only in optimized mode
              dnl (to prevent the effects of possible compiler bugs to affect
              dnl both debug as well as optimized versions)
	      CXXFLAGSO="$CXXFLAGSO -march=native"

	      dnl Also set the mode for f77 compiler
	      F77FLAGSO="$F77FLAGSO -march=native"
              ;;

          intel_icc*)
              dnl Same, but for the icc compiler
              CXXFLAGSO="$CXXFLAGSO -xhost"
              ;;
        esac
        ;;

    *)
        AC_MSG_RESULT(none given or not recognized)
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




dnl -------------------------------------------------------------
dnl Versions of gcc on different platforms use a multitude of flags to
dnl denote thread safe libraries and the like. They are, for example
dnl -threads/-pthread/-mthread, sometimes *thread, sometimes *threads.
dnl On some other platforms, such as Mac OS X, no flags are necessary
dnl and none are understood by the compiler.
dnl
dnl Find out which flag the right one on the present platform
dnl
dnl Usage: DEAL_II_FIND_THREAD_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_GET_THREAD_FLAGS, dnl
[
  case "$target" in
    *apple-darwin*)
	dnl Mac OS X is special in that the compiler generates thread-safe
	dnl code by default, apparently.
	;;

    *)
    	dnl Everything else needs the following setup:
	AC_MSG_CHECKING(for platform specific thread flags)
  	AC_LANG(C++)
  	for i in threads mt pthread pthreads mthreads Kthread kthread invalid_last_entry; do
    	  CXXFLAGS="$CXXFLAGSG -$i"
    	  DEAL_II_TRY_COMPILER_FLAG(
	    [
	     thread_flag="$i"
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
	;;
  esac
])



dnl -------------------------------------------------------------
dnl Test whether multithreading support is requested.
dnl
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
dnl   DEAL_II_CHECK_MULTITHREADING
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MULTITHREADING, dnl
[
  AC_ARG_ENABLE(threads,
              [AS_HELP_STRING([--enable-threads],
              [Use multiple threads inside deal.II])],
    [
      enablethreads="$enableval"
    ],
    [
      dnl Set default to yes unless on cygwin where the TBB is
      dnl currently not yet ported to
      case "$target" in
        *cygwin* )
          enablethreads=no
	  ;;

	* )
          enablethreads=yes
          ;;
      esac
    ])

  if test "$enablethreads" = yes ; then
    dnl On cygwin, the TBB does not compile. Error out.
    case "$target" in
      *cygwin* )
        AC_MSG_ERROR(Multithreading is not supported on CygWin)
        ;;
    esac

    if test "$GXX" = yes ; then
      DEAL_II_GET_THREAD_FLAGS
      DEAL_II_THREAD_CPPFLAGS

      CXXFLAGSG="$CXXFLAGSG -D_REENTRANT"
      CXXFLAGSO="$CXXFLAGSO -D_REENTRANT"
    else
      case "$GXX_VERSION" in
	ibm_xlc)
            CXXFLAGSG="$CXXFLAGSG -qthreaded"
            CXXFLAGSO="$CXXFLAGSO -qthreaded"
            ;;

	compaq_cxx)
            CXXFLAGSG="$CXXFLAGSG -pthread"
            CXXFLAGSO="$CXXFLAGSO -pthread"
	    ;;

        intel_icc*)
	    LDFLAGS="$LDFLAGS -lpthread"
	    ;;

        clang*)
	    LDFLAGS="$LDFLAGS -lpthread"
	    ;;

	pgCC*)
	    LDFLAGS="$LDFLAGS -lpthread"
	    ;;

	*)
            AC_MSG_ERROR(No threading compiler options for this C++ compiler specified at present)
            exit 1
	    ;;
      esac
    fi

    DEAL_II_CHECK_POSIX_THREAD_FUNCTIONS
    AC_DEFINE(DEAL_II_USE_MT_POSIX, 1,
              [Defined if multi-threading is to be achieved by using
               the POSIX functions])

    dnl In any case, we need to link with libdl since the Threading
    dnl Building Blocks require this
    LDFLAGS="$LDFLAGS -ldl"
  fi

  if test "x$enablethreads" != "xno" ; then
    DEAL_II_USE_MT_VAL=1
  else
    DEAL_II_USE_MT_VAL=0
  fi

  AC_DEFINE_UNQUOTED(DEAL_II_USE_MT, $DEAL_II_USE_MT_VAL,
                     [Flag indicating whether the library shall be
                      compiled for multithreaded applications. If so,
		      then it is set to one, otherwise to zero.])
])



dnl -------------------------------------------------------------
dnl On some systems, notably AIX and SUN Solaris, using threads
dnl leads to warnings since the POSIX_MUTEX_INITIALIZER preprocessor
dnl variable used to initialize POSIX mutex objects does not contain
dnl initializers for all elements of the mutex. This is not wrong,
dnl but leads to the message "warning: aggregate has a partly
dnl bracketed initializer", which is annoying since it shows up
dnl _very_ often in our files, although this is something that
dnl happens inside gcc systems headers. So avoid the warning if
dnl necessary
dnl
dnl Usage:
dnl   DEAL_II_CHECK_CHECK_PARTLY_BRACKETED_INITIALIZER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_PARTLY_BRACKETED_INITIALIZER, dnl
[
  if test "$enablethreads" = yes ; then
    case "$GXX_VERSION" in
      gcc*)
  	AC_MSG_CHECKING(for only partly bracketed mutex initializer)
	AC_LANG(C++)
	CXXFLAGS="$CXXFLAGSG -Werror"
	AC_TRY_COMPILE(
   	[
#	include <vector>
	],
        [;],
        [
	  AC_MSG_RESULT(no)
        ],
        [
	  AC_MSG_RESULT(yes)
	  CXXFLAGSG="$CXXFLAGSG -Wno-missing-braces"
	  CXXFLAGSO="$CXXFLAGSO -Wno-missing-braces"
        ])
	;;
      *)
        ;;
    esac
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
  CXXFLAGS="$CXXFLAGSG"
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

  AC_MSG_CHECKING(for posix thread condition variable functions)
  AC_LANG(C++)
  AC_TRY_COMPILE(
   [
#	include <pthread.h>
   ],
   [
	pthread_cond_t   pc;
	pthread_cond_init (&pc, 0);
	pthread_cond_signal (&pc);
	pthread_cond_broadcast (&pc);

        pthread_mutex_t pm;
        pthread_cond_wait (&pc, &pm);
	pthread_cond_destroy (&pc);
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
dnl Check whether some backward compatibility features are disabled
dnl Usage:
dnl   DEAL_II_CHECK_COMPAT_BLOCKER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_COMPAT_BLOCKER, dnl
[
  AC_ARG_ENABLE(compat-blocker,
              [AS_HELP_STRING([--enable-compat-blocker=mapping],
              [Block functions that implicitely assume a Q1 mapping])],
      enable_compat_blocker="$enableval",
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
  CXXFLAGS="$CXXFLAGSG"
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
dnl Check if the declared prototype of abort() has a throw()
dnl specification. We overload abort() in our testsuite, so have
dnl to make sure that we match the exception specification
dnl correctly.
dnl
dnl Usage: DEAL_II_CHECK_ABORT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ABORT, dnl
[
  AC_MSG_CHECKING([for exception specifications on abort()])
  AC_LANG(C++)
  CXXFLAGS="-Werror $CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <cstdlib>
extern "C" void abort () { for(;;) ; }
    ],
    [
    ],
    [
	AC_MSG_RESULT(none)
    ],
    [
	AC_MSG_RESULT(yes)
	AC_DEFINE(DEAL_II_ABORT_NOTHROW_EXCEPTION, 1,
                  [Defined if the prototype of abort() has a no-throw
                   exception specification.])
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
  CXXFLAGS="$2"
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
dnl See if the isfinite function is available in namespace std
dnl (this function is C99 and therefore not part of C++98, but
dnl some implementations provide it nevertheless).
dnl
dnl Usage: DEAL_II_CHECK_ISFINITE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ISFINITE, dnl
[
  AC_MSG_CHECKING(for std::isfinite)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <cmath>
    ],
    [
	double d=0;
	std::isfinite (d);
    ],
    [
	AC_MSG_RESULT(yes)
	AC_DEFINE(DEAL_II_HAVE_ISFINITE, 1,
                  [Defined if std::isfinite is available])
    ],
    [
        AC_MSG_RESULT(no)
    ])
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
  AC_MSG_CHECKING(for rand_r)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
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
  CXXFLAGS="$2"
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
}

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
}

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
}
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
      AC_MSG_RESULT(yes)
    ],
    [
      AC_MSG_RESULT(no)
      AC_MSG_CHECKING(for __func__)
      AC_TRY_COMPILE(
        [
#	include <iostream>
        ],
        [
	  std::cout << __func__ << std::endl;
        ],
        [
          AC_MSG_RESULT(yes)
	  x=__func__
    	],
        [
          AC_MSG_RESULT(no)
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

	template <int dim> T<dim>::T (SubT) {}

	template class T<3>;
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
dnl template <int N, int M> class T      {    int bar ();  };
dnl
dnl template <int M>        class T<1,M> {
dnl   private:
dnl     static int i;
dnl     template <int N1, int N2> friend class T;
dnl };
dnl
dnl template <int N,int M> int T<N,M>::bar () {
dnl   return T<N-1,M>::i;
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
	template <int N, int M> struct T      {    int bar ();  };

	template <int M>        struct T<1,M> {
            T ();
	  private:
	    static int i;
	    template <int N1, int N2> friend class T;
	};

	template <int N,int M> int T<N,M>::bar () {
	  return T<N-1,M>::i;
	}

	template class T<2,1>;
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
dnl     X operator = (T2 &) { return X(); };
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
	    X operator = (T2 &) { return X(); }
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
dnl Some compiler versions, notably ICC, have trouble with the
dnl following code in which we explicitly call a destructor.
dnl This has to be worked around with a typedef. The problem is
dnl that the workaround fails with some other compilers, so that
dnl we can not unconditionally use the workaround...
dnl
dnl Usage: DEAL_II_CHECK_EXPLICIT_DESTRUCTOR_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_EXPLICIT_DESTRUCTOR_BUG, dnl
[
  AC_MSG_CHECKING(for explicit destructor call bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      namespace dealii
      {
        namespace FEValuesViews
        {
          template <int dim, int spacedim> struct Scalar {};
        }

        template <int dim, int spacedim>
        struct X
        {
            FEValuesViews::Scalar<dim,spacedim> scalars[dim*spacedim];

            void f()
              {
                scalars[0].dealii::FEValuesViews::Scalar<dim,spacedim>::~Scalar ();
              }
        };

        template struct X<2,2>;
      }
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE_UNQUOTED(DEAL_II_EXPLICIT_DESTRUCTOR_BUG, 1,
                         [Define if we have to work around a bug where the
			  compiler doesn't accept an explicit destructor call.
                          See the aclocal.m4 file in the top-level directory
                          for a description of this bug.])
    ])
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
  AC_MSG_CHECKING(for 1st template friend in namespace bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
	namespace NS {
	  template <typename T>  class C  {
	      C(const C<T>&);
	      template <typename N> friend class C;
	  };
	}

	namespace NS2 {  template <typename number> class C;  }

	template <typename T> class X {
	    template <typename N> friend class NS2::C;
	    template <typename N> friend class NS::C;
	};

	template class X<int>;

	namespace NS {
  	  template<typename T>
	  inline C<T>::C(const C<T>&)
          {}
	}
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
dnl Another bug in gcc with template and namespaces (fixed since 3.2,
dnl but present in 3.0):
dnl
dnl /* ----------------------------------------------- */
dnl namespace NS {
dnl   template <typename> struct Foo;
dnl }
dnl
dnl class Bar {
dnl     template <typename Y> friend struct NS::Foo;
dnl };
dnl
dnl namespace NS {
dnl   template <typename> struct Foo { Foo (); };
dnl }
dnl
dnl template struct NS::Foo<int>;
dnl /* ----------------------------------------------- */
dnl
dnl gcc2.95 provides a really unhelpful error message, 3.0 gets an
dnl internal compiler error.
dnl
dnl Usage: DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG2
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NAMESP_TEMPL_FRIEND_BUG2, dnl
[
  AC_MSG_CHECKING(for 2nd template friend in namespace bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -Werror"
  AC_TRY_COMPILE(
    [
	namespace NS {
	  template <typename> struct Foo;
	}

	class Bar {
	    template <typename Y> friend struct NS::Foo;
	};

	namespace NS {
	  template <typename> struct Foo { Foo (); };
	}

	template struct NS::Foo<int>;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE_UNQUOTED(DEAL_II_NAMESP_TEMPL_FRIEND_BUG2, 1,
                         [Define if we have to work around another bug in gcc with
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
    int f () { return Y<T>::i; };     // access private field
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
	{}

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
	void X<T>::foo (type) {}

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
dnl gcc 2.95 as well as some other compilers do not correctly implement
dnl the resolution of defect report #45 to the C++ standard (see
dnl http://anubis.dkuug.dk/jtc1/sc22/wg21/docs/cwg_active.html#45).
dnl try to detect this, and set a flag correspondingly. in short,
dnl the DR says that this is allowed, since member classes are
dnl implicitly also friends:
dnl -----------------------------
dnl struct X {
dnl     X ();
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
            X ();
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
dnl gcc up to 3.3 won't accept the following code:
dnl -----------------------------
dnl template <typename> class X {
dnl     template <typename> class Y {};
dnl
dnl     template <typename T>
dnl     template <typename>
dnl     friend class X<T>::Y;
dnl };
dnl
dnl X<int> x;
dnl -----------------------------
dnl
dnl They don't accept the X<T>::Y here, probably because the class is
dnl not complete at this point. gcc3.4 gets it right, though. One
dnl can work around by simply saying
dnl   template <typename> friend class Y;
dnl but then icc doesn't understand this :-( So everyone's got a bug
dnl here. Also, note that the standard says that Y is an implicit friend
dnl of X, but again, many compiler don't implement this correctly, which
dnl is why we have to do something like the above in the first place...
dnl
dnl Usage: DEAL_II_CHECK_NESTED_CLASS_TEMPL_FRIEND_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_NESTED_CLASS_TEMPL_FRIEND_BUG, dnl
[
  AC_MSG_CHECKING(for nested template class friends bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <typename> class X {
          template <typename> class Y {};

          template <typename T>
          template <typename>
          friend class X<T>::Y;
      };

      X<int> x;
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_NESTED_CLASS_TEMPL_FRIEND_BUG, 1,
                     [Defined if the compiler does not understand friend
	              declarations for nested member classes when giving
                      a full class specification.])
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
dnl Usage: DEAL_II_CHECK_MEMBER_VAR_SPECIALIZATION_SPEC_BUG
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
dnl This is a variant of the previous test. Some icc 11.0
dnl builds (sub-releases) on Windows apparently don't allow
dnl the declaration of an explicit specialization of member
dnl arrays of templates:
dnl ---------------------------------
dnl template <int dim>
dnl struct X
dnl {
dnl    static const int N = 2*dim;
dnl    static const int x[N];
dnl };
dnl template <> const int X<2>::x[N];
dnl ---------------------------------
dnl That version of icc requests that there be an initialization,
dnl i.e. it thinks that this is the *definition*, not merely a
dnl *declaration* of an explicit specialization. This is wrong,
dnl however.
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_ARRAY_SPECIALIZATION_SPEC_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_ARRAY_SPECIALIZATION_BUG, dnl
[
  AC_MSG_CHECKING(for member array specialization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <int dim>
      struct X
      {
        static const int N = 2*dim;
        static const int x[N];
      };

      template <> const int X<2>::x[N];
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_MEMBER_ARRAY_SPECIALIZATION_BUG, 1,
                     [Defined if the compiler refuses to allow the
                      explicit specialization of static member
                      arrays. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl gcc 2.95 doesn't like it if we have a member template function
dnl and define it as a template while specializing the outer class
dnl template. This is a nasty bug that is hard to work around...
dnl
dnl Usage: DEAL_II_CHECK_MEMBER_TEMPLATE_SPECIALIZATION_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MEMBER_TEMPLATE_SPECIALIZATION_BUG, dnl
[
  AC_MSG_CHECKING(for template member function specialization bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <int dim> struct X
      {
        template <typename T> void f(T);
      };

      template <>
      template <typename T>
      void X<1>::f (T)
      {}

      template void X<1>::f(int);
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_MEMBER_TEMPLATE_SPECIALIZATION_BUG, 1,
                     [Defined if the compiler refuses to specialize
                      an outer class template while keeping a member
                      as a template. For the exact failure mode, look at
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
	  long double *p1=0, *p2=0;
	  double *p3=0;
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
dnl gcc2.95 (but not later compilers) has a bug with taking the
dnl address of a function with template template parameters (or
dnl with calling this function by specifying explicitly the template
dnl arguments). This requires some working around that in turn does
dnl not work with later compilers.
dnl
dnl Usage: DEAL_II_CHECK_FUNPTR_TEMPLATE_TEMPLATE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_FUNPTR_TEMPLATE_TEMPLATE_BUG, dnl
[
  AC_MSG_CHECKING(for address of template template function bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <int> struct X {};

      template <int dim, template <int> class T>
      void f(T<dim>);

      template <int dim, template <int> class T>
      void* g()
      {
        void (*p) (T<dim>) = &f<dim,T>;
        return (void*)p;
      }

      template void* g<2,X> ();
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_FUNPTR_TEMPLATE_TEMPLATE_BUG, 1,
                     [Defined if the compiler needs a workaround for
		      certain problems with taking the address of
                      template template functions. For the details, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl We compile all the files in the deal.II subdirectory multiple
dnl times, for the various space dimensions. When a program
dnl needs the libs for more than one dimension, the symbols in
dnl these libs must be either static or weak. Unfortunately, some
dnl compilers don't mark functions in anonymous namespaces as
dnl static, and also mangle their names the same way each time,
dnl so we get into trouble with duplicate symbols. Check
dnl whether this is so.
dnl
dnl The check is a little more complicated than the usual one,
dnl since we can't use AC_TRY_COMPILE (we have to compile and link
dnl more than one file).
dnl
dnl Usage: DEAL_II_CHECK_ANON_NAMESPACE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ANON_NAMESPACE_BUG, dnl
[
  AC_MSG_CHECKING(for anonymous namespace and name mangling bug)

  dnl Create the testfile
  echo "namespace { int SYMBOL() {return 1;}; }" >  conftest.cc
  echo "static int f() { return SYMBOL(); }"     >> conftest.cc

  dnl Then compile it twice...
  $CXX -c conftest.cc -o conftest.1.$ac_objext
  $CXX -c conftest.cc -o conftest.2.$ac_objext

  dnl Create a file with main() and also compile it
  echo "int main () {}" > conftest.cc
  $CXX -c conftest.cc -o conftest.3.$ac_objext

  dnl Then try to link everything
  if (($CXX conftest.1.$ac_objext conftest.2.$ac_objext \
            conftest.3.$ac_objext -o conftest  2>&1 ) > /dev/null );\
  then
      AC_MSG_RESULT(no)
  else
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_ANON_NAMESPACE_BUG, 1,
                     [Defined if the compiler needs to see the static
	              keyword even for functions in anonymous namespaces,
                      to avoid duplicate symbol errors when linking.
                      For the details, look at aclocal.m4 in the
                      top-level directory.])
  fi
  rm -f conftest.1.$ac_objext conftest.2.$ac_objext conftest
])



dnl -------------------------------------------------------------
dnl A second test in this direction: if the name of a function is
dnl not mangled differently for each compiler invokation, then
dnl it should at least result in a weak symbol. Test this.
dnl
dnl Note that this is not a problem in itself if the name is
dnl mangled differently each time a file is compiled.
dnl
dnl Usage: DEAL_II_CHECK_ANON_NAMESPACE_BUG2
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ANON_NAMESPACE_BUG2, dnl
[
  AC_MSG_CHECKING(for anonymous namespace and weak linkage bug)

  dnl Create the testfile
  echo "namespace { int SYMBOL() {return 1;}; }" >  conftest.cc
  echo "static int f() { return SYMBOL(); }"     >> conftest.cc

  dnl Compile it
  $CXX -c conftest.cc -o conftest.$ac_objext

  dnl Then look for lines in the output of "nm" that have the name of
  dnl SYMBOL in them. Then make sure that we don't find a line with
  dnl " T " in it, i.e. a text symbol with strong linkage
  check="`nm conftest.$ac_objext | grep SYMBOL | grep ' T '`"

  dnl Then try to link everything
  if test "x$check" = "x" ;
  then
      AC_MSG_RESULT(no)
  else
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG, 1,
                     [Another test if the compiler needs to see the static
	              keyword even for functions in anonymous namespaces,
                      to avoid duplicate symbol errors when linking.
                      For the details, look at aclocal.m4 in the
                      top-level directory.])
  fi
  rm -f conftest.$ac_objext
])



dnl -------------------------------------------------------------
dnl A third test: MIPSpro gives a bogus warning about unused code
dnl on code that is clearly used, when inline member functions
dnl are in classes in anonymous namespaces. Check for this to allow
dnl us to work around this problem
dnl
dnl Usage: DEAL_II_CHECK_ANON_NAMESPACE_BUG3
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ANON_NAMESPACE_BUG3, dnl
[
  AC_MSG_CHECKING(for bogus warning with anonymous namespaces)

  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  case "$GXX_VERSION" in
    gcc*)
	CXXFLAGS="$CXXFLAGSG -Werror"
	;;

    MIPSpro*)
	CXXFLAGS="$CXXFLAGSG -diag_error 1174"
	;;
  esac

  AC_TRY_COMPILE(
   [
     template<int dim> class Point
     {
       public:
         Point() {};
     };

     class GridTools {
       public:
         template <int dim> static void shift (const Point<dim> &s);

         template <typename Predicate>
         static void transform (const Predicate    &predicate);
     };

     namespace {
       template <int dim> class ShiftPoint {
         public:
           ShiftPoint (const Point<dim>   &) {};
           void g() const;
       };

       template<int dim> void ShiftPoint<dim>::g() const {}
     }

     template <typename Predicate>
     void GridTools::transform (const Predicate &predicate)
     { predicate.g(); }

     template <int dim>
     void GridTools::shift (const Point<dim>   &s)
     { transform (ShiftPoint<dim>(s)); }

     template void GridTools::shift (const Point<1>   &);
   ],
   [
	;
   ],
   [
	AC_MSG_RESULT(no)
   ],
   [
	AC_MSG_RESULT(yes)
        AC_DEFINE_UNQUOTED(DEAL_II_ANON_NAMESPACE_BOGUS_WARNING, 1,
                     [Flag indicating whether there is a bug in the
	              compiler that leads to bogus warnings for
	              inline class members in anonymous namespaces])
   ])
])



dnl -------------------------------------------------------------
dnl A test to identify Apples buggy gcc3.3 version. If the
dnl explicit instantiations are placed at the end of a source
dnl file, sometimes only weak symbols are generated, which
dnl lead to linker problems.
dnl For more information, see
dnl http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24331
dnl As this bug can be easily avoided by applying the
dnl "november2004gccupdate" Patch which can be found
dnl on http://www.apple.com/developer, this check will
dnl only produce a warning at the end of the configuration
dnl process.
dnl
dnl Usage: DEAL_II_CHECK_WEAK_LINKAGE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_WEAK_LINKAGE_BUG, dnl
[
  DARWIN_GCC_WEAK_LINKAGE_BUG="no"
  case "$target" in
    *-apple-darwin* )

    AC_MSG_CHECKING(for weak linkage bug (Apple gcc3.3) )

    dnl Create 1st testfile
    echo "template <int> void SYMBOL() {}" >  conftest.cc
    echo "int g() { SYMBOL<1> (); }" >>  conftest.cc
    echo "template void SYMBOL<1> ();"     >> conftest.cc

    dnl Compile it
    $CXX -c conftest.cc -o conftest.$ac_objext

    dnl and write symbols to file 1
    nm conftest.$ac_objext | grep SYMBOL > symb1

    dnl Create 2nd testfile
    echo "template <int> void SYMBOL() {}" >  conftest.cc
    echo "template void SYMBOL<1> ();"     >> conftest.cc
    echo "int g() { SYMBOL<1> (); }" >>  conftest.cc

    dnl Compile it
    $CXX -c conftest.cc -o conftest.$ac_objext

    dnl and write symbols to file 2
    nm conftest.$ac_objext | grep SYMBOL > symb2

    dnl Compare the relevant symbols. According to the C++
    dnl standard, both source codes should produce the
    dnl same symbols.
    check="`diff symb1 symb2`"

    dnl Then try to link everything
    if test "x$check" = "x" ;
    then
        AC_MSG_RESULT(no)
    else
        AC_MSG_RESULT(yes)
        DARWIN_GCC_WEAK_LINKAGE_BUG="yes"
        AC_DEFINE(DEAL_II_WEAK_LINKAGE_BUG, 1,
                       [This error appears in the Apple edition of the
  		        gcc 3.3, which ships with Darwin7.9.0 and
 		        probably previous version. It leads to
		        problems during linking.
                        For the details, look at aclocal.m4 in the
                        top-level directory.])
    fi
    rm -f conftest.$ac_objext
    rm -f symb1 symb2
    ;;
  esac
])



dnl -------------------------------------------------------------
dnl We have so many templates in deal.II that sometimes we need
dnl to make it clear with which types a template parameter can
dnl be instantiated. There is a neat trick to do this: SFINAE
dnl (substitution failure is not an error). The idea is this: the
dnl C++ standard prescribes that a template function is only
dnl considered in a call, if all parts of its signature can be
dnl instantiated with the template parameter replaced by the
dnl respective types/values in this particular call. Example:
dnl   template <typename T>
dnl   typename T::type  foo(T) {...};
dnl   ...
dnl   foo(1);
dnl
dnl The compiler should detect that in this call, the template
dnl parameter T must be identified with the type "int". However,
dnl the return type T::type does not exist. The trick now is
dnl that this is not considered an error: this template is simply
dnl not considered, the compiler keeps on looking for another
dnl possible function foo.
dnl
dnl That allows for a neat trick to rule out a template for certain
dnl template arguments without changing the function signature at
dnl all: Make the return type un-instantiable if the template
dnl type presently considered for instantiation does not qualify.
dnl An example of this is shown below.
dnl
dnl Unfortunately, older compilers do not support this trick: they
dnl issue an error if the return type cannot be installed. In this
dnl case, we back out and just drop the constraint and allow for
dnl all template types. In that case, you'll simply get compiler
dnl error when it tries to compile this template (in case it's in
dnl the .h file), or linker errors for missing functions in case
dnl it's in the .cc file and explicitly instantiated.
dnl
dnl Usage: DEAL_II_CHECK_SFINAE_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_SFINAE_BUG, dnl
[
  AC_MSG_CHECKING(for SFINAE bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <typename> struct int_or_double { static const bool value = false;};
      template <> struct int_or_double<int>    { static const bool value = true; };
      template <> struct int_or_double<double> { static const bool value = true; };

      template <bool, typename> struct constraint_and_return_value {};

      template <typename T> struct constraint_and_return_value<true,T>
      {
          typedef T type;
      };

      // deduction for T=char should file, since return type cannot be
      // instantiated...
      template <typename T>
      typename constraint_and_return_value<int_or_double<T>::value,void>::type
      f (T);

      // ...however, this is not an error. rather, the compiler should
      // instead choose the following function:
      void f(int);
    ],
    [
      f('c');
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. disabling template constraints)
      AC_DEFINE(DEAL_II_SFINAE_BUG, 1,
                     [Defined if the compiler does not support the
		      substitution-failure-is-not-an-error paradigm.
                      For the details, look at aclocal.m4 in the
                      top-level directory.])
    ])
])




dnl -------------------------------------------------------------
dnl Old versions of gcc had a problem with arrays inside ?:
dnl expressions: they decayed too quickly to pointers. This then
dnl leads to erroneous warnings :-(
dnl
dnl Usage: DEAL_II_CHECK_ARRAY_CONDITIONAL_DECAY_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ARRAY_CONDITIONAL_DECAY_BUG, dnl
[
  AC_MSG_CHECKING(for array assignment in conditional bug)
  AC_LANG(C++)
  CXXFLAGS="-W -Wall -Werror"
  AC_TRY_COMPILE(
    [
    ],
    [
  const int x[2][2] = {{1,1},{1,1}};
  const int (&y)[2] = (1 ? x[0] : x[1]);
  return y[0];
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_ARRAY_CONDITIONAL_DECAY_BUG, 1,
                     [Defined if the compiler has a problem with
	              assigning arrays in conditionals])
    ])
])



dnl -------------------------------------------------------------
dnl Old versions of gcc (in particular gcc 3.3.3) has a problem
dnl if an argument to a member function is an array with bounds
dnl that depend on a static const variable inside that class.
dnl
dnl Usage: DEAL_II_CHECK_ARRAY_CONDITIONAL_DECAY_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_ARRAY_ARG_BUG, dnl
[
  AC_MSG_CHECKING(for array argument bug)
  AC_LANG(C++)
  CXXFLAGS="-W -Wall -Werror"
  AC_TRY_COMPILE(
    [
      template <int dim> struct X {
        static const unsigned int N = 1<<dim;
        void f(int (&)[N]);
      };

      template <int dim> void X<dim>::f(int (&)[N]) {}
    ],
    [
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_ARRAY_ARG_BUG, 1,
                     [Defined if the compiler has a problem with
	              using arrays as arguments in functions])
    ])
])



dnl -------------------------------------------------------------
dnl Some versions of gcc get this example wrong:
dnl ---------------------------------
dnl struct X
dnl {
dnl     template <typename T> void operator << (T);
dnl };
dnl
dnl void f()
dnl {
dnl   X x;
dnl   x.operator << <double> (1);
dnl }
dnl ---------------------------------
dnl They want to see a "template" for disambiguation in
dnl    x.template operator << <double> (1);
dnl which shouldn't be necessary since the left hand side of the
dnl dot operator is not template dependent. Surprisingly, this is
dnl only the case for operators, not if operator<< were a regular
dnl function. Annoyingly, other compilers barf on seeing the
dnl disambiguating "template" keyword.
dnl
dnl Usage: DEAL_II_CHECK_TEMPL_OP_DISAMBIGUATION_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TEMPL_OP_DISAMBIGUATION_BUG, dnl
[
  AC_MSG_CHECKING(for template operator disambiguation bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      struct X
      {
          template <typename T> void operator << (T);
      };

      void f()
      {
        X x;
        x.operator << <double> (1);
      }
    ],
    [],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes. using workaround)
      AC_DEFINE(DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG, 1,
                     [Defined if the compiler requires the use of the
                      template keyword for disambiguation keyword in
                      certain contexts in which it is not supposed to
                      do so. For the exact failure mode, look at
                      aclocal.m4 in the top-level directory.])
    ])
])



dnl -------------------------------------------------------------
dnl Some older versions of gcc compile this, despite the 'explicit'
dnl keyword:
dnl ---------------------------------
dnl struct X {
dnl     template <typename T>
dnl     explicit X(T);
dnl };
dnl
dnl void f(X);
dnl
dnl int main () { f(1); }
dnl ---------------------------------
dnl
dnl Check for this misfeature.
dnl
dnl Usage: DEAL_II_CHECK_EXPLICIT_CONSTRUCTOR_BUG
dnl
dnl --------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_EXPLICIT_CONSTRUCTOR_BUG, dnl
[
  AC_MSG_CHECKING(for explicit template constructor bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      struct X {
        template <typename T>
        explicit X(T);
      };

      void f(X);
    ],
    [
      f(1);
    ],
    [
      AC_MSG_RESULT(yes. disabling some functions)
      AC_DEFINE(DEAL_II_EXPLICIT_CONSTRUCTOR_BUG, 1,
                     [Defined if the compiler does not honor the explicit
                      keyword on template constructors.])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Some older versions of gcc deduce pointers to const functions in
dnl template contexts to pointer-to-function of const objects.
dnl This is not correct
dnl
dnl Check for this misfeature.
dnl
dnl Usage: DEAL_II_CHECK_CONST_MEMBER_DEDUCTION_BUG
dnl
dnl --------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CONST_MEMBER_DEDUCTION_BUG, dnl
[
  AC_MSG_CHECKING(for const member deduction bug)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
      template <typename T> struct identity { typedef T type; };

      template <typename C> void new_thread (void (C::*fun_ptr)(),
				       typename identity<C>::type &c);
      template <typename C> void new_thread (void (C::*fun_ptr)() const,
				       const typename identity<C>::type &c);
      struct X { void f() const; };
    ],
    [
      X x;
      new_thread (&X::f, x);
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_CONST_MEMBER_DEDUCTION_BUG, 1,
                     [Defined if the compiler has a bug in deducing
		      the type of pointers to const member functions.])
    ])
])



dnl -------------------------------------------------------------
dnl Check for GCC bug 36052, see
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=36052
dnl
dnl Usage: DEAL_II_CHECK_TYPE_QUALIFIER_BUG
dnl
dnl --------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_TYPE_QUALIFIER_BUG, dnl
[
  case "$GXX_VERSION" in
    gcc*)
      AC_MSG_CHECKING(for warning bug with type qualifiers)
      AC_LANG(C++)
      CXXFLAGS="$CXXFLAGSG -Werror"
      AC_TRY_COMPILE(
        [
          struct S {
              typedef double value_type;
          };

          template <typename T> struct Traits {
              typedef const typename T::value_type dereference_type;
          };

          template <class BlockVectorType> struct ConstIterator {
              typedef typename Traits<BlockVectorType>::dereference_type dereference_type;

              dereference_type operator * () const  { return 0; }
          };

          template class ConstIterator<S>;
        ],
        [
        ],
        [
          AC_MSG_RESULT(no)
        ],
        [
          AC_MSG_RESULT(yes)
          CXXFLAGSG="$CXXFLAGSG -Wno-ignored-qualifiers"
        ])
    ;;
  esac
])




dnl -------------------------------------------------------------
dnl In the gcc libstdc++ headers for std::complex, there is
dnl no defined default copy constructor, but a templated copy
dnl constructor. So when using using a normal assignment
dnl between identical types, the compiler synthesizes the
dnl default operator, rather than using the template.
dnl The code will still be ok, though.
dnl
dnl With -Wsynth in gcc we then get a warning. So if we find that
dnl this is still the case, disable -Wsynth, i.e. remove it from
dnl the list of warning flags.
dnl
dnl This is gcc bug 18644:
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=18644
dnl
dnl Usage: DEAL_II_CHECK_WSYNTH_AND_STD_COMPLEX
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_WSYNTH_AND_STD_COMPLEX, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for problem with -Wsynth and std::complex)
    AC_LANG(C++)
    CXXFLAGS="-Wsynth -Werror"
    AC_TRY_COMPILE(
      [
#       include <complex>
      ],
      [
	std::complex<double> x;
  	x = std::complex<double>(1,0);
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
	CXXFLAGSG="`echo $CXXFLAGSG | perl -pi -e 's/-Wsynth\s*//g;'`"
      ])
  fi
])



dnl -------------------------------------------------------------
dnl Older gcc versions appear to frown upon the way we write the
dnl IsBlockMatrix<MatrixType> template. If that's the case,
dnl remove the -Wctor-dtor-privacy flag.
dnl
dnl This is gcc bug 18644:
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=18644
dnl
dnl Usage: DEAL_II_CHECK_CTOR_DTOR_PRIVACY
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CTOR_DTOR_PRIVACY, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for problem with -Wctor-dtor-privacy)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSG -Werror"
    AC_TRY_COMPILE(
      [
        template <typename T>
        struct IsInt
        {
          private:
            struct yes_type { char c[1]; };
            struct no_type  { char c[2]; };

            static yes_type check_for_int (const int *);

            static no_type check_for_int (...);

          public:
            static const bool value = (sizeof(check_for_int((T*)0))
        			       ==
        			       sizeof(yes_type));
        };

        const bool x = IsInt<double>::value;
      ],
      [
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
	CXXFLAGSG="$CXXFLAGSG -Wno-ctor-dtor-privacy"
      ])
  fi
])



dnl -------------------------------------------------------------
dnl On Mac OS X, gcc appears to have a bug that prevents us from
dnl compiling a bit of code that involves boost::bind. Check for
dnl that.
dnl
dnl Usage: DEAL_II_CHECK_BOOST_BIND_COMPILER_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BOOST_BIND_COMPILER_BUG, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for boost::bind compiler internal error)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSO $BOOST_INCLUDE_DIR"
    AC_TRY_COMPILE(
      [
#include <complex>
#include <iostream>
#include <boost/bind.hpp>

template<typename number>
void bug_function (number test)
{
  std::cout << test << std::endl;
}
      ],
      [
  std::complex<float> float_val (1., 2.);
  boost::bind(&bug_function<std::complex<float> >,
	      float_val)();
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
	AC_DEFINE(DEAL_II_BOOST_BIND_COMPILER_BUG, 1,
	          [Defined if the compiler gets an internal error compiling
		   some code that involves boost::bind])
      ])
  fi
])



dnl -------------------------------------------------------------
dnl Some older versions of GCC (3.x series) have trouble with a
dnl certain part of the BOOST GRAPH library, yielding errors of
dnl the kind "sorry, unimplemented: use of `enumeral_type' in
dnl template type unification".
dnl
dnl Usage: DEAL_II_CHECK_BOOST_GRAPH_COMPILER_BUG
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BOOST_GRAPH_COMPILER_BUG, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for boost::graph compiler internal error)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSO $BOOST_INCLUDE_DIR"
    AC_TRY_COMPILE(
      [
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

namespace types
{
  using namespace boost;

  typedef adjacency_list<vecS, vecS, undirectedS,
			 property<vertex_color_t, default_color_type,
				  property<vertex_degree_t,int> > > Graph;
}

void create_graph (types::Graph)
{}
      ],
      [
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
	AC_DEFINE(DEAL_II_BOOST_GRAPH_COMPILER_BUG, 1,
	          [Defined if the compiler gets an internal error compiling
		   some code that involves boost::graph])
      ])
  fi
])



dnl -------------------------------------------------------------
dnl Check for boost option and find pre-installed boost
dnl -------------------------------------------------------------

AC_DEFUN(DEAL_II_CHECK_BOOST, dnl
[
  AC_ARG_WITH(boost,
              [AS_HELP_STRING([--with-boost=/path/to/boost],
              [Use an installed boost library instead of the contributed one. The optional argument points to the directory containing the boost subdirectory for header files.])],
  DEAL_II_WITH_BOOST($withval))

  CPPFLAGS="$BOOST_INCLUDE_DIR $CPPFLAGS"
  LDFLAGS="$BOOST_LIB_DIR $LDFLAGS"

  AC_CHECK_HEADER(boost/shared_ptr.hpp,,
    [AC_MSG_ERROR([Your boost installation is incomplete!])])
  AC_CHECK_HEADER(boost/type_traits.hpp,,
    [AC_MSG_ERROR([Your boost installation is incomplete!])])
  AC_CHECK_HEADER(boost/tuple/tuple.hpp,,
    [AC_MSG_ERROR([Your boost installation is incomplete!])])
  DEAL_II_CHECK_BOOST_SHARED_PTR_ASSIGNMENT
])

dnl -------------------------------------------------------------
dnl The boost::shared_ptr class has a templated assignment operator
dnl but no assignment operator matching the default operator
dnl signature (this if for boost 1.29 at least). So when using
dnl using a normal assignment between identical types, the
dnl compiler synthesizes teh default operator, rather than using
dnl the template. It doesn't matter here, but is probably an
dnl oversight on behalf of the operators. It should be fixed in newer
dnl versions of boost.
dnl
dnl With -Wsynth in gcc we then get a warning. So if we find that
dnl this is still the case, disable -Wsynth, i.e. remove it from
dnl the list of warning flags.
dnl
dnl Usage: DEAL_II_CHECK_BOOST_SHARED_PTR_ASSIGNMENT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BOOST_SHARED_PTR_ASSIGNMENT, dnl
[
  if test "x$GXX" = "xyes" ; then
    AC_MSG_CHECKING(for boost::shared_ptr assignment operator= template buglet)
    AC_LANG(C++)
    CXXFLAGS="-Wsynth -Werror"
    AC_TRY_COMPILE(
      [
#       include <boost/shared_ptr.hpp>
      ],
      [
        boost::shared_ptr<int> a,b;
        a = b;
      ],
      [
        AC_MSG_RESULT(no)
      ],
      [
        AC_MSG_RESULT(yes)
	CXXFLAGSG="`echo $CXXFLAGSG | perl -pi -e 's/-Wsynth\s*//g;'`"
      ])
  fi
])



dnl -------------------------------------------------------------
dnl Some versions of icc on some platforms issue a lot of warnings
dnl about the unreliability of floating point comparisons. Check
dnl whether we can switch that off by checking whether the compiler
dnl allows -wdXXXX for this warning:
dnl #1572: `floating-point equality and inequality comparisons are
dnl        unreliable'
dnl        (while true, this warning also triggers on comparisons
dnl        with zero, or comparing two variables for which one is
dnl        greater; there is about no way to write numeric code
dnl        without triggering this warning many times over)
dnl
dnl Usage: DEAL_II_ICC_WD_1572
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_ICC_WD_1572, dnl
[
  AC_MSG_CHECKING(whether -wd1572 is allowed)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -wd1572"
  AC_TRY_COMPILE( [], [],
      [
        AC_MSG_RESULT(yes)
	CXXFLAGSG="$CXXFLAGSG -wd1572"
      ],
      [
        AC_MSG_RESULT(no)
      ])
])


dnl -------------------------------------------------------------
dnl Some test, but for the icc C compiler.
dnl
dnl Usage: DEAL_II_ICC_C_WD_1572
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_ICC_C_WD_1572, dnl
[
  AC_MSG_CHECKING(whether -wd1572 is allowed for the C compiler)
  AC_LANG(C)
  OLDCFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -wd1572"
  AC_TRY_COMPILE( [], [],
      [
        AC_MSG_RESULT(yes)
	dnl Keep -wd1572 in CFLAGS
      ],
      [
        AC_MSG_RESULT(no)

	dnl Remove -wd1572 again from flags
	CFLAGS="$OLDCFLAGS"
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
dnl Annoyingly, the Portland Group compiler compiles code with
dnl __builtin_expect just fine, but then doesn't want to link it,
dnl saying it doesn't know this function. So simply not test for
dnl __builtin_expect with that compiler, and to be extrasure make
dnl sure we're doing the right thing also link .
dnl
dnl Usage: DEAL_II_HAVE_BUILTIN_EXPECT
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_BUILTIN_EXPECT, dnl
[
  if test ! "x$GXX_BRAND" = "PortlandGroup" ; then
    AC_MSG_CHECKING(for __builtin_expect)
    AC_LANG(C++)
    CXXFLAGS="$CXXFLAGSG"
    AC_TRY_LINK(
      [
	  bool f() {}
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
  fi
])



dnl -------------------------------------------------------------
dnl Newer versions of gcc have a very nice feature: you can set
dnl a verbose terminate handler, that not only aborts a program
dnl when an exception is thrown and not caught somewhere, but
dnl before aborting it prints that an exception has been thrown,
dnl and possibly what the std::exception::what() function has to
dnl say. Since many people run into the trap of not having a
dnl catch clause in main(), they wonder where that abort may be
dnl coming from.  The terminate handler then at least says what is
dnl missing in their program.
dnl
dnl This test checks whether this feature is available.
dnl
dnl Usage: DEAL_II_HAVE_VERBOSE_TERMINATE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_VERBOSE_TERMINATE, dnl
[
  AC_MSG_CHECKING(for __verbose_terminate_handler)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_LINK(
    [
#include <exception>

namespace __gnu_cxx
{
  extern void __verbose_terminate_handler ();
}

struct preload_terminate_dummy
{
    preload_terminate_dummy()
      { std::set_terminate (__gnu_cxx::__verbose_terminate_handler); }
};

static preload_terminate_dummy dummy;
    ],
    [
	throw 1;
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_VERBOSE_TERMINATE, 1,
                [Define if the compiler provides __verbose_terminate_handler])
    ],
    [
      AC_MSG_RESULT(no)
    ])
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
  AC_MSG_CHECKING(for std::advance warning)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -Werror"
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
dnl
dnl Usage: DEAL_II_CHECK_MIN_VECTOR_CAPACITY
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_MIN_VECTOR_CAPACITY, dnl
[
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG -Werror"

  AC_MSG_CHECKING(for minimal std::vector<T> capacity)
  AC_TRY_RUN(
	[
#include <vector>
int main () {
    std::vector<int> v(1);
    v.reserve (1);
    v.resize (1);
    return v.capacity();
}
	],
	[
	  dnl That's impossible: the return value can't be zero!
	  AC_MSG_ERROR(impossible result -- aborting)
	],
	[
	  result="$?"
	  AC_MSG_RESULT($result)
          AC_DEFINE_UNQUOTED(DEAL_II_MIN_VECTOR_CAPACITY, $result,
                   [Set to the minimal number of elements a std::vector<T> can
	            always hold, i.e. its minimal capacity.])
	]
  )

  dnl Do same thing with std::vector<bool>
  AC_MSG_CHECKING(for minimal std::vector<bool> capacity)
  AC_TRY_RUN(
	[
#include <vector>
int main () {
    std::vector<bool> v(1);
    v.reserve (1);
    v.resize (1);
    return v.capacity();
}
	],
	[
	  dnl That's impossible: the return value can't be zero!
	  AC_MSG_ERROR(impossible result -- aborting)
	],
	[
	  result="$?"
	  AC_MSG_RESULT($result)
          AC_DEFINE_UNQUOTED(DEAL_II_MIN_BOOL_VECTOR_CAPACITY, $result,
                   [Set to the minimal number of elements a std::vector<bool> can
	            always hold, i.e. its minimal capacity.])
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
	return std::numeric_limits<unsigned int>::min();
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
dnl Check whether the <iosfwd> header is available
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
dnl Check whether the std::vector::iterator is just a plain pointer
dnl
dnl Usage: DEAL_II_CHECK_VECTOR_ITERATOR_IS_POINTER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_VECTOR_ITERATOR_IS_POINTER, dnl
[
  AC_MSG_CHECKING(whether vector iterators are plain pointers)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <vector>
template <typename T> void f(T) {}

template void f(int *);
template void f(std::vector<int>::iterator);
    ],
    [
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(DEAL_II_VECTOR_ITERATOR_IS_POINTER, 1,
                [Define if vector iterators are just plain pointers])
    ])
])



dnl -------------------------------------------------------------
dnl Check whether glibc-like stacktrace information is available
dnl for the Exception class. If it is, then try to also determine
dnl whether the compiler accepts the -rdynamic flag, since that is
dnl recommended for linking if one wants to have meaningful
dnl backtraces.
dnl
dnl Usage: DEAL_II_HAVE_GLIBC_STACKTRACE
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_GLIBC_STACKTRACE, dnl
[
  AC_MSG_CHECKING(for glibc-like stacktrace information)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <execinfo.h>
#include <stdlib.h>
    ],
    [
        void * array[25];
        int nSize = backtrace(array, 25);
        char ** symbols = backtrace_symbols(array, nSize);
        free(symbols);
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_GLIBC_STACKTRACE, 1,
                [Define if deal.II is linked against a libc that
                 provides stacktrace debug information that can be
                 printed out in the exception class])

      dnl On Mac OS X, -rdynamic is accepted by the compiler (i.e.
      dnl it doesn't produce an error) but we always get a warning
      dnl that it isn't supported. That's pretty stupid because
      dnl we can't test for it. Consequently, only run the test
      dnl if not on OS X.
      case "$target" in
        *apple-darwin*)
	  ;;

        *)
          AC_MSG_CHECKING(whether compiler accepts -rdynamic)

          CXXFLAGS="$CXXFLAGSG -rdynamic"
          AC_TRY_LINK(
            [],
            [;],
            [
              AC_MSG_RESULT(yes)
              LDFLAGS="$LDFLAGS -rdynamic"
            ],
            [
              AC_MSG_RESULT(no)
            ])
          ;;
      esac
    ],
    [
      AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Check whether the compiler offers a way to demangle symbols
dnl from within the program. Used inside the exception stacktrace
dnl mechanism.
dnl
dnl The example code is taken from
dnl   http://gcc.gnu.org/onlinedocs/libstdc++/18_support/howto.html#6
dnl
dnl Usage: DEAL_II_HAVE_DEMANGLER
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_HAVE_DEMANGLER, dnl
[
  AC_MSG_CHECKING(for libstdc++ demangler)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_LINK(
    [
#include <exception>
#include <iostream>
#include <cxxabi.h>
#include <cstdlib>

struct empty { };

template <typename T, int N>
  struct bar { };
    ],
    [
  int     status;
  char   *realname;

  // exception classes not in <stdexcept>, thrown by the implementation
  // instead of the user
  std::bad_exception  e;
  realname = abi::__cxa_demangle(e.what(), 0, 0, &status);
  std::cout << e.what() << "\t=> " << realname << "\t: " << status << '\n';
  free(realname);


  // typeid
  bar<empty,17>          u;
  const std::type_info  &ti = typeid(u);

  realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
  std::cout << ti.name() << "\t=> " << realname << "\t: " << status << '\n';
  free(realname);

  return 0;
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_LIBSTDCXX_DEMANGLER, 1,
                [Define if the std c++ library provides a demangler conforming
                 to the GCC libstdc++ interface.])
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
  AC_MSG_CHECKING(for definitions of error codes in errno.h)
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



dnl -------------------------------------------------------------
dnl Check whether socket functions (i.e. functions using the
dnl network) affect the numerical accuracy. To spice up things,
dnl on some versions of cygwin (<=1.5.13), calling the windows
dnl socket function (e.g. the functions defined in unistd.h)
dnl will set the internal accuracy of the FPU to 64 bits (double)
dnl instead of 80 bits. Consequently, all long double operations
dnl are performed with double accuracy only. The compiler will
dnl not notice this fact and still report long double limits in
dnl the linits<>-class. We have to check for this, otherwise
dnl some functions like the QGauss constructor will go into an
dnl infinite loop.
dnl
dnl Usage: DEAL_II_CHECK_BROKEN_SOCKETS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BROKEN_SOCKETS, dnl
[
  AC_MSG_CHECKING(for bad socket functions/FPU interaction)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_RUN(
    [
#include <unistd.h>
#include <limits>
#include <iostream>

int main()
{
   char buf[100];
   gethostname(buf,99);

   volatile long double x=1.0;
   x += std::numeric_limits<long double>::epsilon();

   if (x == 1.0)
     return 1;  // this shouldn't happen...
   else
     return 0;
}
    ],
    [
      AC_MSG_RESULT(no)
    ],
    [
      AC_MSG_RESULT([yes. disabling socket functions])
      AC_DEFINE(DEAL_II_BROKEN_SOCKETS, 1,
                [Define if the use of socket functionality leads to strange
                 results with floating point computations on cygwin
                 systems.])
    ])
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
  AC_MSG_CHECKING(for consistency of CXXFLAGSG flags)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
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

  AC_MSG_CHECKING(for consistency of CXXFLAGSO flags)
  CXXFLAGS="$CXXFLAGSO"
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
dnl     Postface: with autoconf 2.63, AC_LANG_PROGRAM(Fortran 77) is called
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
dnl      AC_MSG_ERROR([invalid combination of flags!])
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
dnl      AC_MSG_ERROR([invalid combination of flags!])
dnl      exit 1;
dnl    ])
])




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



dnl -------------------------------------------------------------
dnl Check whether NetCDF is installed. If so we will be able to read
dnl from and write to NetCDF binary or ascii files.
dnl
dnl The NetCDF installation directory (including the lib and include
dnl directory) is given to the --with-netcdf configure option or to
dnl the NETCDF_DIR environment variable
dnl
dnl Usage: DEAL_II_CONFIGURE_NETCDF
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CONFIGURE_NETCDF, dnl
[
  AC_ARG_WITH(netcdf,
              [AS_HELP_STRING([--with-netcdf=/path/to/netcdf],
              [Specify the path to the NetCDF installation, of which the include and library directories are subdirs; use this if you want to override the NETCDF_DIR environment variable.])],
     [
        DEAL_II_NETCDF_DIR="$withval"
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$NETCDF_DIR" != "x" ; then
  	  DEAL_II_NETCDF_DIR="$NETCDF_DIR"
        else
	  DEAL_II_NETCDF_DIR=""
        fi
     ])

  AC_ARG_WITH(netcdf-include,
              [AS_HELP_STRING([--with-netcdf-include=/path/to/netcdf],
              [Specify the path to the NetCDF headers file; use this if you want to override the NETCDF_INCLUDE_DIR environment variable.])],
      [
	 NETCDF_INCLUDE_DIR="$withval"
      ])

  AC_ARG_WITH(netcdf-libs,
              [AS_HELP_STRING([--with-netcdf-libs=/path/to/netcdf],
              [Specify the path to the NetCDF libraries; use this if you want to override the NETCDF_LIBDIR environment variable.])],
      [
	 NETCDF_LIBDIR="$withval"
      ])


  if test "x$DEAL_II_NETCDF_DIR" != "x" ; then
    if test "x$NETCDF_INCLUDE_DIR" != "x" ; then
      CPPFLAGS="-I$NETCDF_INCLUDE_DIR $CPPFLAGS"
      CXXFLAGS="-I$NETCDF_INCLUDE_DIR $CXXFLAGS"
    else
      CPPFLAGS="-I$DEAL_II_NETCDF_DIR/include $CPPFLAGS"
      CXXFLAGS="-I$DEAL_II_NETCDF_DIR/include $CXXFLAGS"
    fi
    if test "x$NETCDF_LIBDIR" != "x" ; then
      LDFLAGS="$LDFLAGS -L$NETCDF_LIBDIR"
    else
      LDFLAGS="-L$DEAL_II_NETCDF_DIR/lib $LDFLAGS"
      if test "$LD_PATH_OPTION" != "no" ; then
	LDFLAGS="$LD_PATH_OPTION$DEAL_II_NETCDF_DIR/lib $LDFLAGS"
      fi
    fi
  else
    if test "x$NETCDF_INCLUDE_DIR" != "x" ; then
      CPPFLAGS="-I$NETCDF_INCLUDE_DIR $CPPFLAGS"
      CXXFLAGS="-I$NETCDF_INCLUDE_DIR $CXXFLAGS"
    fi
    if test "x$NETCDF_LIBDIR" != "x" ; then
      LDFLAGS="$LDFLAGS -L$NETCDF_LIBDIR"
    fi
  fi

  dnl Check for header, if found check for C library,
  dnl if found check for C++ library,
  dnl if successful, HAVE_LIBNETCDF will be set
  AC_CHECK_HEADER(netcdfcpp.h)

  AC_CHECK_LIB(netcdf, nc_open,
  [ DEAL_II_EXTERNAL_LIBS_SAVE_VAL()
    DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT(-lnetcdf_c++ -lnetcdf)
    AC_MSG_CHECKING([for NcFile::NcFile in -lnetcdf_c++])
    AC_LINK_IFELSE(
    [  AC_LANG_PROGRAM([[#include <netcdfcpp.h>
	             ]],
                     [[NcFile test("test")]])
    ],
    [ AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_LIBNETCDF,1,[Define to 1 if you have the `NetCDF' library (-lnetcdf).])
    ],
    [ AC_MSG_RESULT(no)
      DEAL_II_EXTERNAL_LIBS_RESTORE_VAL()
    ])
  ])

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
      2.3*) dnl
        if test ! -d $DEAL_II_PETSC_DIR/lib/$DEAL_II_PETSC_ARCH \
           ; then
          AC_MSG_ERROR([PETSc has not been compiled for the architecture specified with --with-petsc-arch])
        fi
        ;;
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
    AC_DEFINE([DEAL_II_USE_PETSC_VERSION_DEV], [1],
              [Defined if a PETSc installation was found and is not a release])
  else
    DEAL_II_USE_PETSC_VERSION_DEV=no
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
      2.3*) dnl
        AC_TRY_COMPILE(
        [#include "$DEAL_II_PETSC_DIR/bmake/$DEAL_II_PETSC_ARCH/petscconf.h"
        ],
        [#ifdef PETSC_HAVE_MPIUNI
           compile error;
         #endif
        ],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR([PETSc was not built for MPI, but deal.II is!]
        )])
      ;;
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
      2.3*) dnl
        AC_TRY_COMPILE(
        [#include "$DEAL_II_PETSC_DIR/bmake/$DEAL_II_PETSC_ARCH/petscconf.h"
        ],
        [#ifndef PETSC_HAVE_MPIUNI
           compile error;
         #endif
        ],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR([PETSc was built for MPI, but deal.II is not!]
        )])
      ;;
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
    2.3*) dnl
      if test -f $DEAL_II_PETSC_DIR/lib/$DEAL_II_PETSC_ARCH/libmpiuni.a ; then
        DEAL_II_PETSC_MPIUNI_LIB="$DEAL_II_PETSC_DIR/lib/$DEAL_II_PETSC_ARCH/libmpiuni.a"
      fi
      ;;
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
    3.1*)
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
  dnl time. We'll simply grep through it
  DEAL_II_TRILINOS_LIBS="`grep Trilinos_LIBRARIES $DEAL_II_TRILINOS_INCDIR/TrilinosConfig.cmake \
    | perl -pi -e 's/.*\"(.*)\".*/\1/g; s/;/ /g;'`"
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
              [Specify the path to the Metis installation, of which the include and library directories are   subdirs; use this if you want to override the METIS_DIR environment variable.])],
     [
	USE_CONTRIB_METIS=yes
        DEAL_II_METIS_DIR="$withval"
	AC_MSG_RESULT($DEAL_II_METIS_DIR)

        dnl Make sure that what was specified is actually correct
        if test ! -d $DEAL_II_METIS_DIR/Lib ; then
          AC_MSG_ERROR([Path to Metis specified with --with-metis does not point to a complete Metis installation])
	fi

	DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR"
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$METIS_DIR" != "x" ; then
  	  USE_CONTRIB_METIS=yes
          DEAL_II_METIS_DIR="$METIS_DIR"
	  AC_MSG_RESULT($DEAL_II_METIS_DIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_METIS_DIR/Lib ; then
            AC_MSG_ERROR([The path to Metis specified in the METIS_DIR environment variable does not point to a complete Metis installation])
	  fi
	  DEAL_II_METIS_LIBDIR="$DEAL_II_METIS_DIR"
        else
	  USE_CONTRIB_METIS=no
          DEAL_II_METIS_DIR=""
        fi
     ])

  AC_ARG_WITH(metis-libs,
              [AS_HELP_STRING([--with-metis-libs=/path/to/metis],
              [Specify the path to the METIS libraries; use this if you want to override the METIS_LIBDIR environment variable.])],
     [
	USE_CONTRIB_METIS=yes
        DEAL_II_METIS_LIBDIR="$withval"
	AC_MSG_RESULT($DEAL_II_METIS_LIBDIR)

        dnl Make sure that what was specified is actually correct
        if test ! -d $DEAL_II_METIS_LIBDIR ; then
          AC_MSG_ERROR([Path to Metis specified with --with-metis does not point to a complete Metis installation])
	fi
     ],
     [
        dnl Take something from the environment variables, if it is there
        if test "x$METIS_LIBDIR" != "x" ; then
  	  USE_CONTRIB_METIS=yes
          DEAL_II_METIS_LIBDIR="$METIS_LIBDIR"
	  AC_MSG_RESULT($DEAL_II_METIS_LIBDIR)

          dnl Make sure that what this is actually correct
          if test ! -d $DEAL_II_METIS_LIBDIR ; then
            AC_MSG_ERROR([The path to Metis specified in the METIS_DIR environment variable does not point to a complete Metis installation])
	  fi
        else
          dnl Unless --with-metis has been set before, declare that METIS
	  dnl is not desired.
          if test "x$USE_CONTRIB_METIS" != "xyes" ; then
  	    USE_CONTRIB_METIS=no
            DEAL_II_METIS_LIBDIR=""
          fi
        fi
     ])

  if test "x$USE_CONTRIB_METIS" = "xyes" ; then
    AC_DEFINE(DEAL_II_USE_METIS, 1,
              [Defined if a Metis installation was found and is going
               to be used])
    LDFLAGS="$LDFLAGS -L$DEAL_II_METIS_LIBDIR -lmetis"

    dnl AC_MSG_CHECKING(for Metis version)
    dnl DEAL_II_METIS_VERSION=`cat $DEAL_II_METIS_DIR/VERSION`
    dnl AC_MSG_RESULT($DEAL_II_METIS_VERSION)
  fi
])


dnl --------------------------------------------------
dnl What to do if UMFPack is selected
dnl --------------------------------------------------
AC_DEFUN(DEAL_II_WITH_UMFPACK, dnl
[
  AC_ARG_WITH(umfpack,
              [AS_HELP_STRING([--with-umfpack=umfpack-root-directory],
              [Use installed UMFPack version. 'umfpack-root-directory' should be the directory containing directories AMD and UMFPACK. The contributed UMFPack library is used if no argument is given. Default is not to use UMFPack.])],
     [
        AC_MSG_CHECKING(UMFPACK library)
        if test "x$withval" = "xyes" ; then
          USE_CONTRIB_UMFPACK='yes'
          UMFPACK_DIR="`pwd`/contrib/umfpack"
          UMFPACK_INCLUDE_DIR="-I`pwd`/contrib/umfpack/UMFPACK/Include"
          AC_MSG_RESULT(using included version)
          DEAL_II_USE_INTERNAL_UMFPACK=yes
          AC_DEFINE(HAVE_LIBUMFPACK,1,[UMFPACK is $1])
        else
          if test "x$withval" != "xno" ; then
            USE_CONTRIB_UMFPACK='yes'
            UMFPACK_DIR="$withval"
            UMFPACK_INCLUDE_DIR="-I$withval/Include"
            AC_MSG_RESULT(trying version at $withval)
            AC_DEFINE(HAVE_LIBUMFPACK,1,[UMFPACK is $1])
	  else
	    AC_MSG_RESULT(explicitly disabled)
          fi
        fi

     ])

  acx_umfpack=no

  AC_ARG_WITH(umfpack-include,
              [AS_HELP_STRING([--with-umfpack-include=/path/to/UMFPACK],
              [Specify the path to the UMFPACK headers file; use this if you want to override the UMFPACK_INCDIR environment variable.])],
     [
	UMFPACK_INCDIR="$withval"
        acx_umfpack=yes
     ])

  AC_ARG_WITH(umfpack-libs,
              [AS_HELP_STRING([--with-umfpack-libs=/path/to/UMFPACK],
              [Specify the path to the UMFPACK libraries; use this if you want to override the UMFPACK_LIBDIR environment variable.])],
     [
	UMFPACK_LIBDIR="$withval"
        acx_umfpack=yes
     ])

  if test "x$acx_umfpack" = "xyes" ; then
    AC_DEFINE(HAVE_LIBUMFPACK,1,[UMFPACK is $1])
  fi

  if test "x$UMFPACK_DIR" != "x" -a "x$acx_umfpack" = "xno" ; then
    dnl A pathname has been given to --with-umfpack but nothing
    dnl has been specified through the other two flags

    dnl Check whether the libraries are there (unless we use the
    dnl internal version, which we will first have to compile
    dnl before the libs are there)
    if test "x$DEAL_II_USE_INTERNAL_UMFPACK" != "xyes" ; then
      dnl Try old naming scheme for umfpack/amd libraries (before
      dnl Tim Davis incorporated everything into SuiteSparse)
      OLD_LDFLAGS="$LDFLAGS"
      LDFLAGS="-L${UMFPACK_DIR}/UMFPACK $LDFLAGS"
      AC_CHECK_LIB(
        [umfpack],
        [umfpack_di_defaults],
        [
          DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lumfpack)
          if test "$LD_PATH_OPTION" != "no"; then
            LDFLAGS="$LD_PATH_OPTION${UMFPACK_DIR}/UMFPACK $LDFLAGS"
          fi
        ],
        [
          dnl Old naming scheme failed, try the new one
          LDFLAGS="-L${UMFPACK_DIR}/lib $OLD_LDFLAGS"
          AC_CHECK_LIB(
            [umfpack],
            [umfpack_di_defaults],
            [
              DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lumfpack)
              if test "$LD_PATH_OPTION" != "no"; then
                LDFLAGS="$LD_PATH_OPTION${UMFPACK_DIR}/lib $LDFLAGS"
              fi
            ],
            [
              AC_MSG_ERROR(installation of UMFPACK could not be determined)
            ]
          )
        ]
      )

      dnl Now do the same for amd. this one comes second since on the linker
      dnl line, -lumfpack has to preceded -lamd
      OLD_LDFLAGS="$LDFLAGS"
      LDFLAGS="-L${UMFPACK_DIR}/AMD/Lib $LDFLAGS"
      AC_CHECK_LIB(
        [amd],
        [amd_info],
        [
          DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lamd)
          if test "$LD_PATH_OPTION" != "no"; then
            LDFLAGS="$LD_PATH_OPTION${UMFPACK_DIR}/AMD/Lib $LDFLAGS"
          fi
        ],
        [
          dnl Old naming scheme failed, try the new one
          LDFLAGS="-L${UMFPACK_DIR}/lib $OLD_LDFLAGS"
          AC_CHECK_LIB(
            [amd],
            [amd_info],
            [
              DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lamd)
              if test "$LD_PATH_OPTION" != "no"; then
                LDFLAGS="$LD_PATH_OPTION${UMFPACK_DIR}/lib $LDFLAGS"
              fi
            ],
            [
              AC_MSG_ERROR(installation of AMD could not be determined)
            ]
          )
        ]
      )
    fi

    dnl Try old and new naming scheme for header files
    AC_CHECK_FILE(
      [${UMFPACK_DIR}/UMFPACK/Include/umfpack.h],
      [
        UMFPACK_INCLUDE_DIR=-I${UMFPACK_DIR}/UMFPACK/Include
      ],
      [
        AC_CHECK_FILE(
          [${UMFPACK_DIR}/include/suitesparse/umfpack.h],
          [
            UMFPACK_INCLUDE_DIR=-I${UMFPACK_DIR}/include/suitesparse
          ],
          [
            AC_MSG_ERROR(installation of UMFPACK could not be determined)
          ]
        )
      ]
    )

    AC_CHECK_FILE(
      [${UMFPACK_DIR}/AMD/Include/amd.h],
      [
        UMFPACK_INCLUDE_DIR="${UMFPACK_INCLUDE_DIR} -I${UMFPACK_DIR}/AMD/Include"
      ],
      [
        AC_CHECK_FILE(
          [${UMFPACK_DIR}/include/suitesparse/umfpack.h],
          [
            UMFPACK_INCLUDE_DIR="${UMFPACK_INCLUDE_DIR} -I${UMFPACK_DIR}/include/suitesparse"
          ],
          [
            AC_MSG_ERROR(installation of UMFPACK could not be determined)
          ]
        )
      ]
    )

    dnl We also have to see whether the UFconfig.h file can be found
    dnl somewhere. This is not of importance if we use the
    dnl built-in version
    if test "x$DEAL_II_USE_INTERNAL_UMFPACK" != "xyes" ; then
      AC_CHECK_FILE(
        [${UMFPACK_DIR}/UFconfig/UFconfig.h],
        [
          UMFPACK_INCLUDE_DIR="${UMFPACK_INCLUDE_DIR} -I${UMFPACK_DIR}/UFconfig"
        ],
        [
          AC_MSG_ERROR(not found)
        ]
      )
    fi

  else

    if test "x$UMFPACK_INCDIR" != "x" ; then
       dnl Something has been passed to --with-umfpack-include

       UMFPACK_INCLUDE_DIR="-I${UMFPACK_INCDIR}"
       OLD_CXXFLAGS="$CXXFLAGS"
       OLD_CPPFLAGS="$CPPFLAGS"
       CXXFLAGS="-I$UMFPACK_INCDIR $CXXFLAGS"
       CPPFLAGS="-I$UMFPACK_INCDIR $CPPFLAGS"
       AC_CHECK_HEADER(
               [umfpack.h],
               [],
               [AC_MSG_ERROR(installation of UMFPACK could not be determined)]
       )
       AC_CHECK_HEADER([amd.h], [],
                       [AC_MSG_ERROR(installation of UMFPACK could not be determined)])
       AC_CHECK_HEADER([UFconfig.h], [],
                       [AC_MSG_ERROR(installation of UMFPACK could not be determined)])

       if test "x$UMFPACK_LIBDIR" != "x" ; then
          AC_CHECK_LIB(
               [amd],
               [amd_info],
               [
                    DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lamd)
                    if test "$LD_PATH_OPTION" != "no"; then
                        LDFLAGS="$LD_PATH_OPTION${UMFPACK_LIBDIR} $LDFLAGS"
                    fi
               ],
               [AC_MSG_ERROR(installation of AMD could not be determined)]
          )
          AC_CHECK_LIB(
               [umfpack],
               [umfpack_di_defaults],
               [
                    DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lumfpack)
                    if test "$LD_PATH_OPTION" != "no"; then
                        LDFLAGS="$LD_PATH_OPTION${UMFPACK_LIBDIR} $LDFLAGS"
                    fi
               ],
               [AC_MSG_ERROR(installation of UMFPACK could not be determined)]
          )
       else
          AC_CHECK_LIB(
               [amd],
               [amd_info],
               [
                    DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lamd)
               ],
               [AC_MSG_ERROR(installation of AMD could not be determined)]
          )
          AC_CHECK_LIB(
               [umfpack],
               [umfpack_di_defaults],
               [
                    DEAL_II_ADD_EXTERNAL_LIBS_AT_TAIL(-lumfpack)
               ],
               [AC_MSG_ERROR(installation of UMFPACK could not be determined)]
          )
       fi
    fi
  fi
])


dnl --------------------------------------------------
dnl Include the LAPACK library
dnl --------------------------------------------------
AC_DEFUN(DEAL_II_WITH_LAPACK, dnl
[
  if test "x$1" != "xno" ; then
    if test "x$1" != "xyes" ; then lapack="$1"; else lapack="lapack"; fi
    AC_CHECK_LIB($lapack, dgbsv_,
      [ DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT(-l$lapack)
        AC_DEFINE([HAVE_LIBLAPACK], [1],
                  [Defined if deal.II was configured with LAPACK support])
        AC_SUBST(DEAL_II_USE_LAPACK, "yes")
      ],
      [AC_MSG_ERROR([LAPACK library $lapack not found])]
    )
  fi
])

dnl --------------------------------------------------
dnl Print error in BLAS detection
dnl --------------------------------------------------
AC_DEFUN(ABORT_BLAS_ON_ERROR, dnl
[
  AC_MSG_ERROR([[Configuring BLAS library failed although it was requested.
  Most common is one of the following reasons:
    1. A library with name lib$1.a or lib$1.so is not in your library path.
    2. Neither -lgfortran nor -lg2c is in \$F77LIBS.
    3. BLAS requires -lg2c, but only -lgfortran was found.
    4. BLAS requires -lgfortran, but only -lg2c was found.

  In cases 2-4, you may not have the right version of GNU FORTRAN installed
  (<g2c> is the GNU F77 support library, <gfortran> the library for the
  GNU Fortran95 compiler shipped since gcc 4.0) or BLAS has to be recompiled.
]])
])

dnl --------------------------------------------------
dnl Make sure we can link with blas. If something was
dnl given as an argument to --with-blas=xxx, then use
dnl that library, otherwise use 'blas' as library name.
dnl
dnl On Mac OS X, first try something different. OS X has
dnl the concept of a 'framework', which is something like
dnl a collection of libraries. One links with a framework
dnl using '-framework name' on the linker line. On OS X
dnl 10.3, there is a framework called 'vecLib', later
dnl versions have 'Accelerate', both of which include
dnl both the standard blas and cblas libraries. So if no
dnl name is given to --with-blas, and we are on Mac OS X,
dnl try first to link with the Accelerate framework, and
dnl if that fails with the vecLib framework. Only if
dnl that too fails, fall back to the usual rules.
dnl
dnl For more information on frameworks, see here:
dnl http://www.macresearch.org/performance_tutorial_part_i_introducing_accelerate
dnl http://www.macresearch.org/using_the_accelerate_framework
dnl --------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_BLAS_FRAMEWORK, dnl
[
  dnl No argument to --with-blas has been given.
  dnl Try the special arguments for Mac OS X if
  dnl we are on that sort of system
  if (echo $target | grep apple-darwin > /dev/null) ; then
    OLD_LDFLAGS="$LDFLAGS"

    dnl Try to use -framework Accelerate
    AC_MSG_CHECKING(-framework Accelerate)
    LDFLAGS="$OLD_LDFLAGS -framework Accelerate"
    AC_LINK_IFELSE(
              [  AC_LANG_PROGRAM([[extern "C" void daxpy(int,double,double*,int,double*,int);]],
                                 [[daxpy(0,0.,0,0,0,0);]])
              ],
              [
                AC_MSG_RESULT(yes)
                framework_works=yes
              ],
              [
                AC_MSG_RESULT(no)
              ])

    dnl If that didn't work, try to use -framework vecLib
    if test "x$framework_works" != "xyes"; then
      AC_MSG_CHECKING(-framework vecLib)
      LDFLAGS="$OLD_LDFLAGS -framework vecLib"
      AC_LINK_IFELSE(
                [  AC_LANG_PROGRAM([[extern "C" void daxpy(int,double,double*,int,double*,int);]],
                                   [[daxpy(0,0.,0,0,0,0);]])
                ],
                [
                  AC_MSG_RESULT(yes)
                  framework_works=yes
                ],
                [
                  AC_MSG_RESULT(no)
                ])
    fi

    dnl If none of the above worked, revert LDFLAGS to their previous
    dnl values
    if test "x$framework_works" != "xyes"; then
      LDFLAGS="$OLD_LDFLAGS"
    fi
  fi
])

AC_DEFUN(DEAL_II_WITH_BLAS, dnl
[
  if test "x$1" != "xno" ; then
    if test "x$1" != "xyes" ; then
      blas="$1"
      AC_CHECK_LIB($blas, daxpy_,
                   [
                     DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT(-l$blas)
                     AC_DEFINE([HAVE_LIBBLAS], [1],
                               [Defined if deal.II was configured with BLAS support])
                   ],,$F77LIBS)
      AC_SUBST(DEAL_II_USE_BLAS, "yes")
      AC_SUBST(NEEDS_F77LIBS, "yes")
    else
      DEAL_II_CHECK_BLAS_FRAMEWORK
      if test "x$framework_works" != "xyes"; then
        blas="blas";
        AC_CHECK_LIB($blas, daxpy_,
                     [
                       DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT(-l$blas)
                       AC_DEFINE([HAVE_LIBBLAS], [1],
                                 [Defined if deal.II was configured with BLAS support])
                     ],,$F77LIBS)

        AC_SUBST(DEAL_II_USE_BLAS, "yes")
        AC_SUBST(NEEDS_F77LIBS, "yes")
      fi
    fi
  fi
])



dnl --------------------------------------------------
dnl What to do if external boost is selected
dnl --------------------------------------------------
AC_DEFUN(DEAL_II_WITH_BOOST, dnl
[
  if test "x$1" != "xyes" ; then
    BOOST_INCLUDE_DIR="-I$1"
    if test "$LD_PATH_OPTION" != "no" ; then
      BOOST_LIB_DIR="$LD_PATH_OPTION$1/lib"
    fi
  else
    BOOST_INCLUDE_DIR=''
    BOOST_LIB_DIR=''
  fi
])


dnl --------------------------------------------------
dnl Include the GSL library
dnl --------------------------------------------------
AC_DEFUN(DEAL_II_WITH_ZLIB, dnl
[
  if test "x$1" != "xyes" ; then
    zlib=$1
  else
    zlib=z
  fi

  dnl See if we can find the function crc32 in libz.so
  AC_CHECK_LIB($zlib, crc32,
    [
      dnl Yes, we can. Now also check whether we can do
      dnl #include <zlib.h>
      AC_CHECK_HEADER(zlib.h,
        [
          DEAL_II_ADD_EXTERNAL_LIBS_AT_FRONT(-l$zlib)
          AC_DEFINE(HAVE_LIBZ,[],"")
        ])
    ])
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

    if test ! -d "${use_p4est}/DEBUG" -o ! -d "${use_p4est}/FAST" ; then
    echo "${use_p4est}/DEBUG"
      AC_MSG_ERROR([p4est directories $use_p4est/DEBUG or $use_p4est/FAST not found])
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
    export USE_CONTRIB_P4EST

    DEAL_II_P4EST_DIR=${use_p4est}
    export DEAL_II_P4EST_DIR

    CXXFLAGSG="$CXXFLAGSG -I$use_p4est/DEBUG/include"
    CXXFLAGSO="$CXXFLAGSO -I$use_p4est/FAST/include"
  else
    AC_MSG_RESULT(no)
  fi
])


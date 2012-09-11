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

  dnl Then icc came along and started to identify itself as
  dnl    icpc version 12.1.0 (gcc version 4.2.1 compatibility)
  dnl which also doesn't help...
  if test "$GXX" = "yes" ; then
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "icpc"`
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

                          dnl Maybe CLang?
                          is_clang="`($CXX --version 2>&1) | grep clang`"
                          if test "x$is_clang" != x ; then
                            AC_MSG_RESULT(C++ compiler is clang)
                            GXX_BRAND=clang
                            GXX_VERSION=clang
                            GXX_VERSION_DETAILED="$GXX_VERSION"
                          else

                            dnl Maybe Cray C++?
                            is_cray="`($CXX -V 2>&1) | grep Cray`"
                            if test "x$is_cray" != x ; then
                              AC_MSG_RESULT(C++ compiler is Cray C++)
                              GXX_BRAND=cray
                              GXX_VERSION=cray
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
    CXXFLAGSO="$CXXFLAGSO -O2 -funroll-loops -funroll-all-loops -fstrict-aliasing -Wuninitialized -felide-constructors"
    CXXFLAGSG="$CXXFLAGSG -DDEBUG -Wall -W -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Wswitch"

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

    dnl See whether the gcc we use already has a flag for C++2011 features.
    DEAL_II_CHECK_CXX1X


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

    dnl Newer versions of gcc can pass a flag to the assembler to
    dnl compress debug sections. At the time of writing this test,
    dnl this can save around 230 MB of disk space on the object
    dnl files we produce (810MB down to 570MB for the debug versions
    dnl of object files). Unfortunately, the sections have to be
    dnl unpacked again when they are put into the shared libs, so
    dnl no savings there.
    dnl
    dnl The flag also doesn't appear to be working on Cygwin, as
    dnl per email by John Fowkes on the mailing list in Feb 2012,
    dnl so don't run the test on cygwin.
    case "$target" in
      *cygwin* )
         ;;

      * )
         CXXFLAGS="-Wa,--compress-debug-sections"
         AC_MSG_CHECKING([whether the assembler understands -Wa,--compress-debug-sections])
         AC_TRY_LINK(
          [
          ],
          [;],
          [
            AC_MSG_RESULT(yes)
            CXXFLAGSG="$CXXFLAGSG -Wa,--compress-debug-sections"
          ],
          [
            AC_MSG_RESULT(no)
          ])
         ;;
    esac

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

          CFLAGS="$CFLAGS -m64"
          CFLAGSG="$CFLAGSG -m64"
          CFLAGSO="$CFLAGSO -m64"

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
#           include <stdio.h>
#           include <unistd.h>
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
          dnl Set flags for IBM's xlC compiler. As of version 11.1, it still
          dnl doesn't grok all of deal.II, but it's getting closer. The
          dnl -qxflag=EnableIssue214PartialOrdering flag is necessary to
          dnl resolve code like this that we have in the MemoryConsumption
          dnl namespace:
          dnl ----------------------------------------------
          dnl template <typename T> void f (T* const);
          dnl template <typename T> void f (const T &);
          dnl
          dnl void g() {
          dnl  int *p;
          dnl  f(p);
          dnl }
          dnl ----------------------------------------------
          dnl
          dnl Similarly, -qxflag=IgnoreCVOnTopOfFunctionTypes is necessary for
          dnl ----------------------------------------------
          dnl struct S {};
          dnl S & foo (const S&);
          dnl
          dnl struct X {
          dnl   template <typename T>
          dnl   void f (const T &t) const;
          dnl };
          dnl
          dnl void
          dnl print_summary ()
          dnl {
          dnl   X x;
          dnl   x.f (foo);
          dnl }
          dnl ----------------------------------------------
          dnl We have this kind of code with S=std::ostream, foo=std::endl
          dnl and X=LogStream when we do things like
          dnl     deallog << std::endl;
          CXXFLAGSG="$CXXFLAGSG -DDEBUG -check=bounds -info=all -qrtti=all -qsuppress=1540-2907 -qsuppress=1540-2909 -qxflag=EnableIssue214PartialOrdering -qxflag=IgnoreCVOnTopOfFunctionTypes"
          CXXFLAGSO="$CXXFLAGSO -O2 -w -qansialias -qrtti=all -qsuppress=1540-2907 -qsuppress=1540-2909 -qxflag=EnableIssue214PartialOrdering -qxflag=IgnoreCVOnTopOfFunctionTypes"
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
          dnl #1565: attributes are ignored on a class declaration that is not
	  dnl        also a definition (this happens in BOOST in a number of
	  dnl        places)
          CXXFLAGSG="$CXXFLAGSG -w1 -wd175 -wd525 -wd327 -wd424 -wd11 -wd734 -wd858 -wd1565"
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

	  dnl Intel's MPI implementation based on ICC requires that
	  dnl mpi.h be included *before* things like <stdio.h> or
	  dnl <iostream>. How they envision this to work is beyond
	  dnl me because they may be indirectly included from other
	  dnl header files. Besides this, autoconf generates tests
	  dnl that don't follow this rule and so fail at ./configure
	  dnl time. There is nothing we can do about it. However,
	  dnl there is a workaround described here:
	  dnl   http://software.intel.com/en-us/articles/intel-mpi-library-for-linux-running-list-of-known-issues/#A3
	  dnl We switch this on if we use Intel's ICC + MPI
	  if test "x$DEAL_II_USE_MPI" = "xyes" ; then
	    CXXFLAGSG="$CXXFLAGSG -DMPICH_IGNORE_CXX_SEEK"
	    CXXFLAGSO="$CXXFLAGSO -DMPICH_IGNORE_CXX_SEEK"
	  fi

          dnl Finally, see if the compiler supports C++0x
          DEAL_II_CHECK_CXX1X
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

          dnl pgCC can't (as of writing this, with version 12.5 in mid-2012) compile a part of BOOST.
          dnl Fortunately, BOOST provides a workaround by setting a specific preprocessor
          dnl symbol that can be set. Do so.
          CXXFLAGSG="$CXXFLAGSG -DBOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE"
          CXXFLAGSO="$CXXFLAGSO -DBOOST_MPL_CFG_NO_HAS_XXX_TEMPLATE"
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
dnl See if there is a compiler flag to enable C++1x features
dnl
dnl Usage: DEAL_II_CHECK_CXX1X
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_CXX1X, dnl
[
  AC_MSG_CHECKING(whether compiler has a flag to support C++2011)

  dnl See if -std=c++0x exists
  OLD_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS=-std=c++0x
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

  dnl If the flag above works, then see if the compiler is complete
  dnl enough in this area
  if test "x$test_cxx1x" = "xyes" ; then
    DEAL_II_CHECK_CXX1X_COMPONENTS("-std=c++0x")
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

  AC_MSG_CHECKING(for std::type_traits)
  AC_TRY_COMPILE(
       [#include <type_traits>],
       [ const bool m0 = std::is_trivial<double>::value;
         const bool m1 = std::is_standard_layout<double>::value;
         const bool m2 = std::is_pod<double>::value; ],
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
              [Defined if the compiler we use supports the C++2011 standard
               well enough to allow using the standard library classes instead
               of the corresponding BOOST classes.])


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

    dnl Intel icc 12.1 has this crazy behavior where it needs -std=c++0x
    dnl to compile BOOST, but it fails every single one of the header
    dnl file tests above. So we end up here. Work around this by using
    dnl the flag even though we can't use a single piece of functionality.
    if test "x$GXX_VERSION" = "xintel_icc12" ; then
      CXXFLAGSG="$CXXFLAGSG $1"
      CXXFLAGSO="$CXXFLAGSO $1"
    fi
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
    dnl Verify that we indeed have a compiler that identifies
    dnl itself as GCC
    CC_VERSION_STRING=`($CC -v 2>&1) | grep "gcc version"`
    if test "x$CC_VERSION_STRING" = "x" ; then
      GCC=no
    fi

    dnl Then icc came along and started to identify itself as
    dnl    icc version 12.1.0 (gcc version 4.2.1 compatibility)
    dnl which also doesn't help...
    CC_VERSION_STRING=`($CC -v 2>&1) | grep "icc"`
    if test "x$CC_VERSION_STRING" != "x" ; then
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

                      is_clang="`($CC --version 2>&1) | grep clang`"
                      if test "x$is_clang" != x ; then
                        AC_MSG_RESULT(C compiler is clang)
                        CC_VERSION=clang
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

      clang)
          CFLAGS="$CFLAGS -g"
          CFLAGSO="$CFLAGS -fast -O2"
          CFLAGSPIC="-fPIC"
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
              dnl Tune for this processor
              CXXFLAGSG="$CXXFLAGSG -march=$withcpu"
              CXXFLAGSO="$CXXFLAGSO -march=$withcpu"

              dnl Also set the mode for f77 compiler
              F77FLAGSG="$F77FLAGSG -march=$withcpu"
              F77FLAGSO="$F77FLAGSO -march=$withcpu"
          ;;
        esac
        ;;

    native)
        AC_MSG_RESULT(native processor variant)
        case "$GXX_VERSION" in
          gcc*)
              dnl Tune for this processor
              CXXFLAGSG="$CXXFLAGSG -march=native"
              CXXFLAGSO="$CXXFLAGSO -march=native"

              dnl Also set the mode for f77 compiler
              F77FLAGSG="$F77FLAGSG -march=native"
              F77FLAGSO="$F77FLAGSO -march=native"
              ;;

          intel_icc*)
              dnl Same, but for the icc compiler
              CXXFLAGSO="$CXXFLAGSO -xhost"
              CXXFLAGSG="$CXXFLAGSG -xhost"
              ;;
        esac
        ;;

    *)
        AC_MSG_RESULT(none given or not recognized)
        ;;
  esac
])



dnl -------------------------------------------------------------
dnl Check whether the compiler allows for vectorization and that
dnl vectorization actually works. For this test, we use compiler
dnl intrinsics similar to what is used in the deal.II library and
dnl check whether the arithmetic operations are correctly performed
dnl on examples where all numbers are exactly represented as
dnl floating point numbers.
dnl
dnl Usage: DEAL_II_COMPILER_VECTORIZATION_LEVEL
dnl 0 means no vectorization, 1 support for SSE2, 2 support for AVX
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_DETECT_VECTORIZATION_LEVEL, dnl
[
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  dnl SSE2 check in debug mode
  AC_MSG_CHECKING(whether CPU supports SSE2)
  AC_TRY_RUN(
    [
#include <emmintrin.h>
#include <mm_malloc.h>
        int main()
        {
        __m128d a, b;
        const unsigned int vector_bytes = sizeof(__m128d);
        const int n_vectors = vector_bytes/sizeof(double);
        __m128d * data =
          reinterpret_cast<__m128d*>(_mm_malloc (2*vector_bytes, vector_bytes));
        double * ptr = reinterpret_cast<double*>(&a);
        ptr[0] = (volatile double)(1.0);
        for (int i=1; i<n_vectors; ++i)
          ptr[i] = 0.0;
        b = _mm_set1_pd ((volatile double)(2.25));
        data[0] = _mm_add_pd (a, b);
        data[1] = _mm_mul_pd (b, data[0]);
        ptr = reinterpret_cast<double*>(&data[1]);
        unsigned int return_value = 0;
        if (ptr[0] != 7.3125)
          return_value = 1;
        for (int i=1; i<n_vectors; ++i)
          if (ptr[i] != 5.0625)
            return_value = 1;
        _mm_free (data);
        return return_value;
        }
    ],
    [
      AC_MSG_RESULT(yes)
      dnl AVX check in debug mode
      AC_MSG_CHECKING(whether CPU supports AVX)
      AC_TRY_RUN(
      [
#include <immintrin.h>
#include <mm_malloc.h>
        int main()
        {
        __m256d a, b;
        const unsigned int vector_bytes = sizeof(__m256d);
        const int n_vectors = vector_bytes/sizeof(double);
        __m256d * data =
          reinterpret_cast<__m256d*>(_mm_malloc (2*vector_bytes, vector_bytes));
        double * ptr = reinterpret_cast<double*>(&a);
        ptr[0] = (volatile double)(1.0);
        for (int i=1; i<n_vectors; ++i)
          ptr[i] = 0.0;
        b = _mm256_set1_pd ((volatile double)(2.25));
        data[0] = _mm256_add_pd (a, b);
        data[1] = _mm256_mul_pd (b, data[0]);
        ptr = reinterpret_cast<double*>(&data[1]);
        unsigned int return_value = 0;
        if (ptr[0] != 7.3125)
          return_value = 1;
        for (int i=1; i<n_vectors; ++i)
          if (ptr[i] != 5.0625)
            return_value = 1;
        _mm_free (data);
        }
      ],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(DEAL_II_COMPILER_VECTORIZATION_LEVEL, 2,
                  [Equal to 0 in the generic case,
                   equal to 1 if CPU compiled for supports SSE2,
                   equal to 2 if CPU compiled for supports AVX])
      ],
      [
        AC_MSG_RESULT(no)
        AC_DEFINE(DEAL_II_COMPILER_VECTORIZATION_LEVEL, 1,
                  [Equal to 0 in the generic case,
                   equal to 1 if CPU compiled for supports SSE2,
                   equal to 2 if CPU compiled for supports AVX])
      ])
    ],
    [
        AC_DEFINE(DEAL_II_COMPILER_VECTORIZATION_LEVEL, 0,
                  [Equal to 0 in the generic case,
                   equal to 1 if CPU compiled for supports SSE2,
                   equal to 2 if CPU compiled for supports AVX])
        AC_MSG_RESULT(no)
    ])
])



dnl -------------------------------------------------------------
dnl Check whether the compiler allows to use arithmetic operations
dnl +-*/ on vectorized data types or whether we need to use
dnl _mm_add_pd for addition and so on. +-*/ is preferred because
dnl it allows the compiler to choose other optimizations like
dnl fused multiply add, whereas _mm_add_pd explicitly enforces the
dnl assembler command.
dnl
dnl Usage: DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DEAL_II_CHECK_VECTOR_ARITHMETICS, dnl

[
  AC_MSG_CHECKING(whether compiler supports vector arithmetics)
  AC_LANG(C++)
  CXXFLAGS="$CXXFLAGSG"
  AC_TRY_COMPILE(
    [
#include <emmintrin.h>
    ],
    [
        __m128d a, b;
        a = _mm_set_sd (1.0);
        b = _mm_set1_pd (2.1);
        __m128d c = a + b;
        __m128d d = b - c;
        __m128d e = c * a + d;
        __m128d f = e/a;
        (void)f;
    ],
    [
        AC_MSG_RESULT(yes)
        AC_DEFINE(DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS, 1,
                  [Defined if the compiler can use arithmetic operations on
                  vectorized data types])
    ],
    [
        AC_MSG_RESULT(no)
    ])
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
#       include <iostream>
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
#       include <iostream>
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

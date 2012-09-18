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
          dnl #497: "declaration of "dim" hides template parameter" (while
	  dnl       theoretically useful, pgCC unfortunately gets this one
	  dnl       wrong on legitimate code where no such parameter is
	  dnl       hidden, see the email by ayaydemir on 9/3/2012)
          CXXFLAGSG="$CXXFLAGSG -DDEBUG -g --display_error_number --diag_suppress 68 --diag_suppress 111 --diag_suppress 128 --diag_suppress 155 --diag_suppress 177 --diag_suppress 175 --diag_suppress 185 --diag_suppress 236 --diag_suppress 284 --diag_suppress 497"
          CXXFLAGSO="$CXXFLAGSO -fast -O2 --display_error_number --diag_suppress 68 --diag_suppress 111 --diag_suppress 128 --diag_suppress 155 --diag_suppress 177 --diag_suppress 175 --diag_suppress 185 --diag_suppress 236 --diag_suppress 284 --diag_suppress 497"
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
          CFLAGSO="$CFLAGS -O2"
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

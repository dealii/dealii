Shall we still support DEC OSF? *mhm* :-]




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




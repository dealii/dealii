# $Id$
# Copyright by the deal.II authors, 1998, 1999, 2000, 2001, 2002, 2010, 2011, 2012


ifneq ($(CLEAN),yes)

# deal.II top level directory
D=..


# get options which hold for all files of the project
include $D/common/Make.global_options

# top level targets
default all: debug optimized


# list the directories and the various kinds of files, the ones with the
# largest compile jobs first so that 'make -j N' saturates many processors
# also towards the end of compiling (e.g. numerics/vectors.cc takes several
# minutes to compile...)
all-dirs      := numerics fe dofs lac base grid hp multigrid distributed

cc-files      := $(shell for i in $(all-dirs) ; do echo $D/source/$$i/*.cc ; done)
h-files       := $(sort $(shell echo $D/include/deal.II/*/*.h))

# build unique object file names
tmp1     := $(shell echo $(cc-files) | $(PERL) -pi -e 's,$D/source/,,g; s,/,_,g;')
o-files  := $(addprefix $(LIBDIR)/optimized/, $(tmp1:.cc=.$(OBJEXT)) )
go-files := $(addprefix $(LIBDIR)/debug/, $(tmp1:.cc=.$(OBJEXT)))


debug: $(LIBDIR)/libdeal_II.g$(lib-suffix)
optimized: $(LIBDIR)/libdeal_II$(lib-suffix)



# rules how to generate object files from source files.
$(LIBDIR)/debug/%.$(OBJEXT) :
	@echo "======================debug======$(MT)== $(<:$D/source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
$(LIBDIR)/optimized/%.$(OBJEXT) :
	@echo "======================optimized==$(MT)== $(<:$D/source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@



%.inst : %.inst.in $D/common/template-arguments
	@echo "======================instantiate==== $(@F)"
	@$D/common/scripts/expand_instantiations $D/common/template-arguments < $< > $@ || exit 1

%.inst : %.inst.pl $D/common/dealiitemplates.pm
	@echo "======================instantiate==== $(@F)"
	@perl -I $D/common/ $< > $@ || exit 1

#Todo: remove this target after renaming files
%.inst : %.pl $D/common/dealiitemplates.pm
	@echo "======================instantiate==== $(@F)"
	@perl -I $D/common/ $< > $@ || exit 1



# there are extra .o file for the function parser and umfpack that we
# need to link into out library
ifeq ($(enable-parser),yes)
  extra-o-files   := $(LIBDIR)/contrib/functionparser/fparser.$(OBJEXT)
  extra-g.o-files := $(LIBDIR)/contrib/functionparser/fparser.$(OBJEXT)
endif

ifeq ($(USE_CONTRIB_UMFPACK),yes)
  umfpack-files   := $(sort $(shell echo $(LIBDIR)/contrib/umfpack/*.$(OBJEXT)))
  extra-o-files   += $(umfpack-files)
  extra-g.o-files += $(umfpack-files)
endif



# in general, when linking shared libraries, we will want to link with
# the underlying libs as well. AIX even requires this. on the other
# hand, on DEC OSF with CXX this will fail with duplicate symbols, so
# make sure that it doesn't happen on that system and that compiler
ifneq ($(GXX-VERSION),compaq_cxx)
  # if we use HSL, then we need to link against libhsl.so.
  ifeq ($(USE_CONTRIB_HSL),yes)
    deplibs.g += $(lib-contrib-hsl)
    deplibs.o += $(lib-contrib-hsl)
  endif

  # same with PETSC
  ifeq ($(USE_CONTRIB_PETSC),yes)
    deplibs.g += $(lib-contrib-petsc.g) \
	         $(call DEAL_II_RPATH,$(DEAL_II_PETSC_LIBDIR))
    deplibs.o += $(lib-contrib-petsc.o) \
	         $(call DEAL_II_RPATH,$(DEAL_II_PETSC_LIBDIR))
  endif

  # same with SLEPC
  ifeq ($(USE_CONTRIB_SLEPC),yes)
    deplibs.g += $(lib-contrib-slepc) \
	         $(call DEAL_II_RPATH,$(DEAL_II_SLEPC_LIBDIR))
    deplibs.o += $(lib-contrib-slepc) \
	         $(call DEAL_II_RPATH,$(DEAL_II_SLEPC_LIBDIR))
  endif

  # same with TRILINOS
  ifeq ($(USE_CONTRIB_TRILINOS),yes)
    deplibs.g += $(lib-contrib-trilinos) \
	         $(call DEAL_II_RPATH,$(DEAL_II_TRILINOS_LIBDIR))
    deplibs.o += $(lib-contrib-trilinos) \
	         $(call DEAL_II_RPATH,$(DEAL_II_TRILINOS_LIBDIR))
  endif

  # same with MUMPS
  ifeq ($(USE_CONTRIB_MUMPS),yes)
    deplibs.g += $(lib-contrib-mumps)
    deplibs.o += $(lib-contrib-mumps)
  endif

  # and METIS
  ifeq ($(USE_CONTRIB_METIS),yes)
    deplibs.g += $(lib-contrib-metis)
    deplibs.o += $(lib-contrib-metis)
  endif

  # also link with TBB
  ifeq ($(enable-threads),yes)
    deplibs.g += $D/lib/libtbb_debug$(shared-lib-suffix)
    deplibs.o += $D/lib/libtbb$(shared-lib-suffix)
  endif

  # and p4est
  ifeq ($(USE_CONTRIB_P4EST),yes)
    deplibs.g += $(DEAL_II_P4EST_DIR)/DEBUG/$(DEAL_II_P4EST_LIBDIR_NAME)/libp4est$(shared-lib-suffix) \
                 $(DEAL_II_P4EST_DIR)/DEBUG/$(DEAL_II_P4EST_LIBDIR_NAME)/libsc$(shared-lib-suffix) \
		 $(call DEAL_II_RPATH,$(DEAL_II_P4EST_DIR)/DEBUG/$(DEAL_II_P4EST_LIBDIR_NAME))
    deplibs.o += $(DEAL_II_P4EST_DIR)/FAST/$(DEAL_II_P4EST_LIBDIR_NAME)/libp4est$(shared-lib-suffix) \
                 $(DEAL_II_P4EST_DIR)/FAST/$(DEAL_II_P4EST_LIBDIR_NAME)/libsc$(shared-lib-suffix) \
		 $(call DEAL_II_RPATH,$(DEAL_II_P4EST_DIR)/FAST/$(DEAL_II_P4EST_LIBDIR_NAME))
  endif

  # for cygwin and darwin, also lapack and such are necessary
  # we better link all libraries mentioned in $(LIBS)
  ifneq ($(findstring cygwin,$(TARGET))$(findstring darwin,$(TARGET)),)
    deplibs.g += $(LIBS)
    deplibs.o += $(LIBS)
  endif

else
  deplibs.g =
  deplibs.o =
endif


# rules for static libraries
$(LIBDIR)/libdeal_II.g$(static-lib-suffix): $(go-files) $(extra-g.o-files)
	@echo "=====deal.II==========debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files) $(extra-g.o-files)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II$(static-lib-suffix): $(o-files) $(extra-o-files)
	@echo "=====deal.II==========optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files) $(extra-o-files)
	@$(RANLIB) $@


# Rules for shared libraries. these are much more difficult to link
# because the linker has to update relocations etc. This requires
# seeking back and forth in memory mapped files and is excruciatingly
# slow on remotely mounted file systems. Consequently, we link in a
# local file in the $TMPDIR directory, and only after the fact move
# them to their final location
$(LIBDIR)/libdeal_II.g$(shared-lib-suffix): $(go-files) $(extra-g.o-files)
	@echo "=====deal.II==========debug======$(MT)== Linking library:   $(@F)"
	@rm -f $@
	@WORKFILE=`mktemp $${TMPDIR:-/tmp}/tmp.XXXXXXXXXX` && \
		$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $$WORKFILE $(call DEAL_II_ADD_SONAME,deal_II.g) $(go-files) $(extra-g.o-files) $(deplibs.g) && \
	 	mv $$WORKFILE $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II.g) $@

$(LIBDIR)/libdeal_II$(shared-lib-suffix): $(o-files) $(extra-o-files)
	@echo "=====deal.II==========optimized==$(MT)== Linking library:   $(@F)"
	@rm -f $@
	@WORKFILE=`mktemp $${TMPDIR:-/tmp}/tmp.XXXXXXXXXX` && \
		$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $$WORKFILE $(call DEAL_II_ADD_SONAME,deal_II) $(o-files) $(extra-o-files) $(deplibs.o) && \
		mv $$WORKFILE $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II) $@



# generate a file that contains prototypes for BLAS and LAPACK functions
../include/deal.II/lac/lapack_templates.h: ../include/deal.II/lac/lapack_templates.h.in \
						$D/common/scripts/lapack_templates.pl
	@echo ===================================== Remaking $@
	@$(PERL) $D/common/scripts/lapack_templates.pl ../include/deal.II/lac/lapack_templates.h.in > $@


# A rule to make sure the various */Makefile.dep files have the right
# dependencies on their source files. Use the opportunity to also ensure that
# .inst and the lapack file are built before everything else
Makefile.dep: ../include/deal.II/lac/lapack_templates.h \
		$(cc-files)
	@echo "===================================== Remaking $@"
	@for i in $(all-dirs) ; do \
		echo $$i/Makefile.dep: $$i/*.cc  \
					$(h-files) \
					Makefile \
					$D/common/Make.global_options ; \
	 done >> $@
	@echo "lac/Makefile.dep: $D/include/deal.II/lac/lapack_templates.h.in \
				 $D/include/deal.II/lac/lapack_templates.h.in" \
	 >> $@
include Makefile.dep


# Rule to generate the dependency files, one for each source
# directory. These file are automagically remade whenever needed,
# i.e. whenever one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion below.
#
# If the command fails, then remove Makefile.dep again and fail
%/Makefile.dep:
	@echo "===================================== Remaking $@"
	@$D/common/scripts/make_dependencies -n $(INCLUDE) "-B\$$(LIBDIR)" \
			$(filter $D/source/$(dir $@)%, $(cc-files)) \
	 | $(PERL) -pe 's!LIBDIR\)/(debug|optimized)/(.*):!LIBDIR)/$$1/$(@:%/Makefile.dep=%)_$$2:!g;' \
	 > $@
	@if test ! -s $@ ; then \
	  echo "*** make_dependencies failed to produce a valid file $@" ; \
	  rm $@ ; \
	  false ; \
	fi


# include all the dependencies
include $(addsuffix /Makefile.dep, $(all-dirs))

endif



################### some clean-up rules

clean:
	-rm -f Makefile.dep *~ */*~ */*/*~ */Makefile.dep source/*/*.inst \
		$D/lib/{debug,optimized}/*.o


################### coverage
# A target that creates a visual representation of which files and
# lines are covered by the testsuite, the example programs, and
# possible other executables. To use this rule, configure deal.II
# as usual, then add "-fprofile-arcs -ftest-coverage" to the following
# variables in common/Make.global_options:
# - SHLIBFLAGS
# - LDFLAGS
# - CXXFLAGS.g
# Then compile the library and run all the executables for which you
# want it reported that. Note, however, that you can't run them in
# parallel or they will step on each other's toes when writing output.
# In other words, the testsuite has to be run with -j1, however long
# that may take.
#
# When all of these programs have been run, then do
#   cd source
#   make coverage
# This may take a long time: both running lcov and, in particular,
# genhtml can take hours. In the latter case, there isn't even an
# indication that anything is happening, but it is -- just be
# patient, let it run overnight.
coverage:
	for dir in $(all-dirs) ; do \
	  cd $$dir ; \
	  for file in *.cc ; do \
	    basename=`echo $$file | perl -pi -e 's/\\.cc$$//;'` ; \
	    for ext in o gcda gcno ; do \
	      rm -f $${basename}.$$ext ; \
	      if test -f $D/lib/debug/$${dir}_$${basename}.$$ext ; then \
	        ln $D/lib/debug/$${dir}_$${basename}.$$ext $${basename}.$$ext ; \
	      fi ; \
	    done ; \
	  done ; \
	  cd .. ; \
	done
	lcov --base-directory . --directory . -c -o d2.info
	lcov --remove d2.info "/usr*" -o d2.info
	rm -rf ../coverage
	genhtml --no-branch-coverage -o ../coverage -t "deal.II test coverage" --num-spaces 8 d2.info


###################

.PHONY: default all debug optimized TAGS
.PHONY: clean

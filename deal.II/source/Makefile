# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002, 2010


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
all-dirs      := numerics fe dofs grid hp multigrid lac base distributed

cc-files      := $(shell for i in $(all-dirs) ; do echo $D/source/$$i/*.cc ; done)
h-files       := $(sort $(shell echo $D/include/deal.II/*/*.h))

inst-in-files := $(shell echo */*.inst.in)
inst-files    := $(inst-in-files:%.in=%)


# build unique object file names
tmp1               := $(shell echo $(cc-files) | $(PERL) -pi -e 's,$D/source/,,g; s,/,_,g;')
tmp-o-files-indep  := $(addprefix $(LIBDIR)/optimized/, $(tmp1:.cc=.$(OBJEXT)) )
tmp-go-files-indep := $(addprefix $(LIBDIR)/debug/, $(tmp1:.cc=.$(OBJEXT)))

o-files  := $(tmp-o-files-indep)
go-files := $(tmp-go-files-indep)


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
	@echo "===================================== $(@F)"
	@$D/common/scripts/expand_instantiations $D/common/template-arguments < $< > $@



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
    deplibs.g += $(lib-contrib-petsc.g)
    deplibs.o += $(lib-contrib-petsc.o)
  endif

  # same with SLEPC
  ifeq ($(USE_CONTRIB_SLEPC),yes)
    deplibs.g += $(lib-contrib-slepc)
    deplibs.o += $(lib-contrib-slepc)
  endif

  # same with TRILINOS
  ifeq ($(USE_CONTRIB_TRILINOS),yes)
    deplibs.g += $(lib-contrib-trilinos)
    deplibs.o += $(lib-contrib-trilinos)
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
    deplibs.g += $(DEAL_II_P4EST_DIR)/DEBUG/lib/libp4est$(shared-lib-suffix) \
                 $(DEAL_II_P4EST_DIR)/DEBUG/lib/libsc$(shared-lib-suffix) \
		 -Wl,-rpath,$(DEAL_II_P4EST_DIR)/DEBUG/lib
    deplibs.o += $(DEAL_II_P4EST_DIR)/FAST/lib/libp4est$(shared-lib-suffix) \
                 $(DEAL_II_P4EST_DIR)/FAST/lib/libsc$(shared-lib-suffix) \
		 -Wl,-rpath,$(DEAL_II_P4EST_DIR)/FAST/lib
  endif

  # for cygwin and darwin, also llapack and such are necessary
  # we better link all libraries mentioned in $(LIBS)
  ifneq ($(findstring cygwin,$(TARGET))$(findstring darwin,$(TARGET)),)
    deplibs.g += $(LIBS)
    deplibs.o += $(LIBS)
  endif

else
  deplibs.g =
  deplibs.o =
endif


# rules how to make the libraries themselves
$(LIBDIR)/libdeal_II.g$(static-lib-suffix): $(go-files)
	@echo "=====deal.II==========debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files) $(extra-g.o-files)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II$(static-lib-suffix): $(o-files)
	@echo "=====deal.II==========optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files) $(extra-o-files)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II.g$(shared-lib-suffix): $(go-files)
	@echo "=====deal.II==========debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II.g) $(call DEAL_II_ADD_SONAME,deal_II.g) $(go-files) $(extra-g.o-files) $(deplibs.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II.g) $@

$(LIBDIR)/libdeal_II$(shared-lib-suffix): $(o-files)
	@echo "=====deal.II==========optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II) $(call DEAL_II_ADD_SONAME,deal_II) $(o-files) $(extra-o-files) $(deplibs.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II) $@



# generate a file that contains prototypes for BLAS and LAPACK functions
../include/deal.II/lac/lapack_templates.h: ../include/deal.II/lac/lapack_templates.h.in \
						$D/common/scripts/lapack_templates.pl
	@echo ===================================== Remaking $@
	@$(PERL) $D/common/scripts/lapack_templates.pl ../include/deal.II/lac/lapack_templates.h.in > $@


# A rule to make sure the various */Makefile.dep files have the right
# dependencies on their source files. Use the opportunity to also ensure that
# .inst and the lapack file are built before everything else
Makefile.dep: $(inst-files) ../include/deal.II/lac/lapack_templates.h \
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
	@(($D/common/scripts/make_dependencies -n $(INCLUDE) "-B\$$(LIBDIR)" \
			$(filter $D/source/$(dir $@)%, $(cc-files)) \
		| $(PERL) -pe 's!LIBDIR\)/(debug|optimized)/(.*):!LIBDIR)/$$1/$(@:%/Makefile.dep=%)_$$2:!g;' \
	  ) > $@) \
	 || (rm -f $@ ; false)


# include all the dependencies
include $(addsuffix /Makefile.dep, $(all-dirs))

endif



################### some clean-up rules

clean:
	-rm -f Makefile.dep *~ */*~ */*/*~ */Makefile.dep source/*/*.inst \
		$D/lib/{debug,optimized}/*.o


###################

.PHONY: default all debug optimized TAGS
.PHONY: clean

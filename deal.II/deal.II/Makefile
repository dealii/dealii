# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# deal.II top level directory
D=..


# get options which hold for all files of the project
include $D/common/Make.global_options



# top level targets
default: 2dg 2d
all: 2dg 2d 3dg 3d 1dg 1d
debug: 2dg 3dg 1dg


1dg: $(LIBDIR)/libdeal_II_1d.g$(lib-suffix)
1d : $(LIBDIR)/libdeal_II_1d$(lib-suffix)
2dg: $(LIBDIR)/libdeal_II_2d.g$(lib-suffix)
2d : $(LIBDIR)/libdeal_II_2d$(lib-suffix)
3dg: $(LIBDIR)/libdeal_II_3d.g$(lib-suffix)
3d : $(LIBDIR)/libdeal_II_3d$(lib-suffix)

# create lists of file names. compile vectors.cc first because it is
# one of the files that takes the longest
cc-files     = source/numerics/vectors.cc $(shell echo source/*/*.cc)
h-files      = $(sort $(shell echo include/*/*.h))

inst-in-files= $(shell echo source/*/*.inst.in)
inst-files   = $(inst-in-files:%.in=%)


# replace subdirectories in cc file names by a prefix, preprend a "lib/[g]o"
# and change the suffix
tmp1        = $(shell echo $(cc-files) | $(PERL) -pi -e 's,source/,,g; s,/,_,g;')
o-files-1d  = $(addprefix $(LIBDIR)/1d/, $(tmp1:.cc=_1d.$(OBJEXT)) )
go-files-1d = $(addprefix $(LIBDIR)/1d/, $(tmp1:.cc=_1d.g.$(OBJEXT)))
o-files-2d  = $(addprefix $(LIBDIR)/2d/, $(tmp1:.cc=_2d.$(OBJEXT)) )
go-files-2d = $(addprefix $(LIBDIR)/2d/, $(tmp1:.cc=_2d.g.$(OBJEXT)))
o-files-3d  = $(addprefix $(LIBDIR)/3d/, $(tmp1:.cc=_3d.$(OBJEXT)) )
go-files-3d = $(addprefix $(LIBDIR)/3d/, $(tmp1:.cc=_3d.g.$(OBJEXT)))

o-files  = $(o-files-1d)  $(o-files-2d)  $(o-files-3d)
go-files = $(go-files-1d) $(go-files-2d) $(go-files-3d)


# rules how to generate object files from source files. note that
# there are some files which needs exceptional rules; these are
# listed immediately below
$(LIBDIR)/1d/%.g.$(OBJEXT) :
	@echo "=====deal.II====1d====debug======$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=1 -c $< -o $@
$(LIBDIR)/1d/%.$(OBJEXT) :
	@echo "=====deal.II====1d====optimized==$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=1 -c $< -o $@

$(LIBDIR)/2d/%.g.$(OBJEXT) :
	@echo "=====deal.II====2d====debug======$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=2 -c $< -o $@
$(LIBDIR)/2d/%.$(OBJEXT) :
	@echo "=====deal.II====2d====optimized==$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=2 -c $< -o $@

$(LIBDIR)/3d/%.g.$(OBJEXT) :
	@echo "=====deal.II====3d====debug======$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=3 -c $< -o $@
$(LIBDIR)/3d/%.$(OBJEXT) :
	@echo "=====deal.II====3d====optimized==$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=3 -c $< -o $@

source/%.inst : source/%.inst.in $D/common/template-arguments
	@echo "=====deal.II========================= $(@F)"
	@$D/common/scripts/expand_instantiations $D/common/template-arguments < $< > $@



# in general, when linking shared libraries, we will want to link with
# the underlying libs as well. AIX even requires this. on the other
# hand, on DEC OSF with CXX this will fail with duplicate symbols, so
# make sure that it doesn't happen on that system and that compiler
ifneq ($(GXX-VERSION),compaq_cxx)
  deplibs.g = $(lib-lac.g) $(lib-base.g)
  deplibs.o = $(lib-lac.o) $(lib-base.o)

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

  # and METIS
  ifeq ($(USE_CONTRIB_METIS),yes)
    deplibs.g += $(lib-contrib-metis)
    deplibs.o += $(lib-contrib-metis)
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
$(LIBDIR)/libdeal_II_1d.g$(static-lib-suffix): $(go-files-1d)
	@echo "=====deal.II====1d====debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files-1d)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II_1d$(static-lib-suffix): $(o-files-1d)
	@echo "=====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files-1d)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II_2d.g$(static-lib-suffix): $(go-files-2d)
	@echo "=====deal.II====2d====debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files-2d)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II_2d$(static-lib-suffix): $(o-files-2d)
	@echo "=====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files-2d)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II_3d.g$(static-lib-suffix): $(go-files-3d)
	@echo "=====deal.II====3d====debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files-3d)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II_3d$(static-lib-suffix): $(o-files-3d)
	@echo "=====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files-3d)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II_1d.g$(shared-lib-suffix): $(go-files-1d)
	@echo "=====deal.II====1d====debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II_1d.g) $(call DEAL_II_ADD_SONAME,deal_II_1d.g) $(go-files-1d) $(deplibs.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II_1d.g) $@

$(LIBDIR)/libdeal_II_1d$(shared-lib-suffix): $(o-files-1d)
	@echo "=====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II_1d) $(call DEAL_II_ADD_SONAME,deal_II_1d) $(o-files-1d) $(deplibs.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II_1d) $@


$(LIBDIR)/libdeal_II_2d.g$(shared-lib-suffix): $(go-files-2d)
	@echo "=====deal.II====2d====debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II_2d.g) $(call DEAL_II_ADD_SONAME,deal_II_2d.g) $(go-files-2d) $(deplibs.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II_2d.g) $@

$(LIBDIR)/libdeal_II_2d$(shared-lib-suffix): $(o-files-2d)
	@echo "=====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II_2d) $(call DEAL_II_ADD_SONAME,deal_II_2d) $(o-files-2d) $(deplibs.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II_2d) $@


$(LIBDIR)/libdeal_II_3d.g$(shared-lib-suffix): $(go-files-3d)
	@echo "=====deal.II====3d====debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II_3d.g) $(call DEAL_II_ADD_SONAME,deal_II_3d.g) $(go-files-3d) $(deplibs.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II_3d.g) $@

$(LIBDIR)/libdeal_II_3d$(shared-lib-suffix): $(o-files-3d)
	@echo "=====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,deal_II_3d) $(call DEAL_II_ADD_SONAME,deal_II_3d) $(o-files-3d) $(deplibs.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,deal_II_3d) $@



# Rule to generate the dependency file. This file is
# automagically remade whenever needed, i.e. whenever
# one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom
# of this file.
#
# If the command fails, then remove Makefile.dep again and fail
Makefile.dep: $(cc-files) $(h-files) \
		$(inst-files) \
		Makefile $D/common/Make.global_options
	@echo "===================================== Remaking deal.II/Makefile"
	@(for dir in dofs fe grid hp multigrid numerics ; do $D/common/scripts/make_dependencies $(INCLUDE) "-B\$$(LIBDIR)" `echo source/$$dir/*cc` \
		| $(PERL) -p -e 's!LIBDIR\)/(.*):!LIBDIR)/DIM_PLACEHOLDER/$$1:!g;' \
		| $(PERL) -pe 's!((\.g)?.$(OBJEXT)):!_DIM_PLACEHOLDER$$1:!g;' \
	        | $(PERL) -pe 's!^(.*)/DIM_PLACEHOLDER/(.*)_DIM_PLACEHOLDER(\..*):!$$1/1d/DIRPLACEHOLDER_$$2_1d$$3 $$1/2d/DIRPLACEHOLDER_$$2_2d$$3 $$1/3d/DIRPLACEHOLDER_$$2_3d$$3:!g;'\
		| $(PERL) -pe "s/DIRPLACEHOLDER/$$dir/g;" \
		; done ) > $@ \
	 || (rm -f $@ ; false)


# include all the dependencies
include Makefile.dep




################### some clean-up rules

clean:
	-rm -f *~ */*~ */*/*~ Makefile.dep


###################

.PHONY: all deal.II TAGS
.PHONY: clean
.PHONY: default 1d 2d 3d all

# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# deal.II top level directory
D=..

# get options which hold for all files of the project
include $D/common/Make.global_options


# create lists of file names
cc-files     = $(shell echo source/*.cc)
o-files-base = $(sort $(addprefix $(LIBDIR)/lac/, $(patsubst source/%,%,$(cc-files:.cc=.$(OBJEXT)))))
o-files      = $(o-files-base)
go-files     = $(o-files-base:.$(OBJEXT)=.g.$(OBJEXT))
h-files      = $(sort $(shell echo include/lac/*.h))

inst-in-files= $(shell echo source/*.inst.in)
inst-files   = $(inst-in-files:%.in=%)

# if umfpack was configured, simply link all the umfpack files into liblac.
# these object files should have been generated before we got to the present
# directory
ifeq ($(USE_CONTRIB_UMFPACK),yes)
  umfpack-files = $(sort $(shell echo $(LIBDIR)/contrib/umfpack/*.$(OBJEXT)))
  o-files += $(umfpack-files)
  go-files += $(umfpack-files)
endif


# production rules
$(LIBDIR)/lac/%.g.$(OBJEXT) :
	@echo "=====lac==============debug======$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
$(LIBDIR)/lac/%.$(OBJEXT) :
	@echo "=====lac==============optimized==$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@

source/%.inst : source/%.inst.in $D/common/template-arguments
	@echo "=====lac============================= $(@F)"
	@$D/common/scripts/expand_instantiations $D/common/template-arguments < $< > $@


# rules
lib: $(LIBDIR)/liblac.g$(lib-suffix) $(LIBDIR)/liblac$(lib-suffix)
libg: $(LIBDIR)/liblac.g$(lib-suffix)
libo: $(LIBDIR)/liblac$(lib-suffix)

# in general, when linking shared libraries, we will want to link with
# the underlying libs as well. AIX even requires this. on the other
# hand, on DEC OSF with CXX this will fail with duplicate symbols, so
# make sure that it doesn't happen on that system and that compiler
ifneq ($(GXX-VERSION),compaq_cxx)
  deplibs.g = $(lib-base.g)
  deplibs.o = $(lib-base.o)

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


$(LIBDIR)/liblac$(static-lib-suffix): $(o-files)
	@echo "=====lac==============optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files)
	@$(RANLIB) $@

$(LIBDIR)/liblac.g$(static-lib-suffix): $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files)
	@$(RANLIB) $@


$(LIBDIR)/liblac$(shared-lib-suffix): $(o-files)
	@echo "=====lac==============optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,lac) $(call DEAL_II_ADD_SONAME,lac) $(o-files) $(deplibs.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,lac) $@

$(LIBDIR)/liblac.g$(shared-lib-suffix): $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,lac.g) $(call DEAL_II_ADD_SONAME,lac.g) $(go-files) $(deplibs.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,lac.g) $@


clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep source/*.inst


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
	@echo ===================================== Remaking lac/Makefile
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B$(LIBDIR)/lac" $(cc-files) \
		> $@ \
	 || (rm -f $@ ; false)


# include all the dependencies
include Makefile.dep

.PHONY: clean
.PHONY: lib lib$(static-lib-suffix) lib.g$(static-lib-suffix)

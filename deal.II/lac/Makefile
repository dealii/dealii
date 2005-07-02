# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# deal.II top level directory
D=..


# create lists of file names
cc-files     = $(shell echo source/*.cc)
o-files-base = $(addprefix $(LIBDIR)/lac/, $(patsubst source/%,%,$(cc-files:.cc=.$(OBJEXT))))
o-files      = $(sort $(o-files-base))
go-files     = $(o-files:.$(OBJEXT)=.g.$(OBJEXT))
h-files      = $(sort $(shell echo include/lac/*.h))

# get options which hold for all files of the project
include $D/common/Make.global_options


# production rules
ifneq ($(with-multithreading),no)
  MT = MT
else
  MT = ==
endif

$(LIBDIR)/lac/%.g.$(OBJEXT) :
	@echo "=====lac==============debug======$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
$(LIBDIR)/lac/%.$(OBJEXT) :
	@echo "=====lac==============optimized==$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@


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

  # and UMFPACK
  ifeq ($(USE_CONTRIB_UMFPACK),yes)
    deplibs.g += $(LIBDIR)/liblac_umfpack$(lib-suffix)
    deplibs.o += $(LIBDIR)/liblac_umfpack$(lib-suffix)
  endif

  # for cygwin, also llapack and such are necessary
  # we better link all libraries mentioned in $(libs)
  ifeq ($(findstring cygwin,$(TARGET)),cygwin)
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
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files) $(deplibs.o)

$(LIBDIR)/liblac.g$(shared-lib-suffix): $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files) $(deplibs.g)



clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep


# Rule to generate the dependency file. This file is
# automagically remade whenever needed, i.e. whenever
# one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom
# of this file.
#
# If the file should turn out to be empty, then blow it
# away to let make issue a proper error message rather than
# some obscure follow-up problems
Makefile.dep: $(cc-files) $(h-files) Makefile $D/common/Make.global_options
	@echo ============================ Remaking lac/Makefile
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B$(LIBDIR)/lac" $(cc-files) \
		> Makefile.dep
	@if test -s $@ ; then : else rm $@ ; fi


# include all the dependencies
include Makefile.dep

.PHONY: clean
.PHONY: lib lib$(static-lib-suffix) lib.g$(static-lib-suffix)

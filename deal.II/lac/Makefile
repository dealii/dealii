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
ifeq ($(enable-shared),yes)
  lib: $(LIBDIR)/liblac.g.so $(LIBDIR)/liblac.so
  libg: $(LIBDIR)/liblac.g.so
  libo: $(LIBDIR)/liblac.so
else
  lib: $(LIBDIR)/liblac.g.a $(LIBDIR)/liblac.a
  libg: $(LIBDIR)/liblac.g.a
  libo: $(LIBDIR)/liblac.a
endif

# if we use HSL, then we need to link against libhsl.so. if we don't, then
# don't set this variable
ifeq ($(USE_CONTRIB_HSL),yes)
  LIBHSL=$(lib-contrib-hsl)
else
  LIBHSL=
endif

# in general, when linking shared libraries, we will want to link with
# the underlying libs as well. AIX even requires this. on the other
# hand, on DEC OSF with CXX this will fail with duplicate symbols, so
# make sure that it doesn't happen on that system and that compiler
ifneq ($(GXX-VERSION),compaq_cxx)
  deplibs.g = $(lib-base.g) $(LIBHSL)
  deplibs.o = $(lib-base.o) $(LIBHSL)

  ifeq ($(USE_CONTRIB_PETSC),yes)
    deplibs.g += $(lib-contrib-petsc.g)
    deplibs.o += $(lib-contrib-petsc.o)
  endif

else
  deplibs.g =
  deplibs.o =
endif


$(LIBDIR)/liblac.a: $(o-files)
	@echo "=====lac==============optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files)
	@$(RANLIB) $@

$(LIBDIR)/liblac.g.a: $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files)
	@$(RANLIB) $@


$(LIBDIR)/liblac.so: $(o-files)
	@echo "=====lac==============optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files) $(deplibs.o)

$(LIBDIR)/liblac.g.so: $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files) $(deplibs.g)



clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(h-files) Makefile $D/common/Make.global_options
	@echo ============================ Remaking lac/Makefile
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B$(LIBDIR)/lac" $(cc-files) \
		> Makefile.dep

# include all the dependencies
include Makefile.dep

.PHONY: clean
.PHONY: lib lib.a lib.g.a

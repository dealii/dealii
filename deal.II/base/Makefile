# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002



# deal.II top level directory
D=..

# get options which hold for all files of the project
ifneq ($(CLEAN),yes)
include $D/common/Make.global_options

# create lists of file names
cc-files     = $(shell echo $D/base/source/*.cc)
o-files-base = $(addprefix $(LIBDIR)/base/, $(patsubst $D/base/source/%,%,$(cc-files:.cc=.$(OBJEXT))))
o-files      = $(sort $(o-files-base))
go-files     = $(o-files:.$(OBJEXT)=.g.$(OBJEXT))
h-files      = $(sort $(shell echo include/base/*.h))

# there's an extra .o file for the function parser that we need to link into
# libbase. do similarly for the threading building block things if threading
# is enabled
ifeq ($(enable-parser),yes)
  extra-o-files   = $(LIBDIR)/contrib/functionparser/fparser.$(OBJEXT)
  extra-g.o-files = $(LIBDIR)/contrib/functionparser/fparser.$(OBJEXT)
endif


# production rules
$(LIBDIR)/base/%.g.$(OBJEXT) :
	@echo "=====base=============debug======$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
$(LIBDIR)/base/%.$(OBJEXT) :
	@echo "=====base=============optimized==$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@


# in general, when linking shared libraries, we will want to link with
# the underlying libs as well. AIX even requires this. on the other
# hand, on DEC OSF with CXX this will fail with duplicate symbols, so
# make sure that it doesn't happen on that system and that compiler
#
# for the current directory, this is only the PETSc libraries since we
# reference them in base/source/utilities
ifneq ($(GXX-VERSION),compaq_cxx)
  deplibs.g =
  deplibs.o =

  ifeq ($(USE_CONTRIB_PETSC),yes)
    deplibs.g += $(lib-contrib-petsc.g) $(DEAL_II_PETSC_MPIUNI_LIB)
    deplibs.o += $(lib-contrib-petsc.o) $(DEAL_II_PETSC_MPIUNI_LIB)
  endif

  ifeq ($(USE_CONTRIB_TRILINOS),yes)
    deplibs.g += $(lib-contrib-trilinos)
    deplibs.o += $(lib-contrib-trilinos)
  endif

  ifeq ($(enable-threads),yes)
    deplibs.g += $D/lib/libtbb_debug$(shared-lib-suffix)
    deplibs.o += $D/lib/libtbb$(shared-lib-suffix)
  endif

  ifeq ($(USE_CONTRIB_P4EST),yes)
    deplibs.g += $(DEAL_II_P4EST_DIR)/DEBUG/lib/libp4est.so \
                 $(DEAL_II_P4EST_DIR)/DEBUG/lib/libsc.so \
		 -Wl,-rpath,$(DEAL_II_P4EST_DIR)/DEBUG/lib
    deplibs.o += $(DEAL_II_P4EST_DIR)/FAST/lib/libp4est.so \
                 $(DEAL_II_P4EST_DIR)/FAST/lib/libsc.so \
		 -Wl,-rpath,$(DEAL_II_P4EST_DIR)/FAST/lib
  endif

else
  deplibs.g =
  deplibs.o =
endif

# targets for object code without linking into library
objg: $(go-files)
objo: $(o-files)

# rules for libraries
lib: $(LIBDIR)/libbase.g$(lib-suffix) $(LIBDIR)/libbase$(lib-suffix)
libg: $(LIBDIR)/libbase.g$(lib-suffix)
libo: $(LIBDIR)/libbase$(lib-suffix)

$(LIBDIR)/libbase$(static-lib-suffix): $(o-files) $(extra-o-files)
	@echo "=====base=============optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files) $(extra-o-files)
	@$(RANLIB) $@

$(LIBDIR)/libbase.g$(static-lib-suffix): $(go-files) $(extra-g.o-files)
	@echo "=====base=============debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files) $(extra-g.o-files)
	@$(RANLIB) $@


$(LIBDIR)/libbase$(shared-lib-suffix): $(o-files) $(extra-o-files)
	@echo "=====base=============optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,base) $(call DEAL_II_ADD_SONAME,base) $(o-files) $(extra-o-files)  $(deplibs.o)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,base) $@

$(LIBDIR)/libbase.g$(shared-lib-suffix): $(go-files) $(extra-g.o-files)
	@echo "=====base=============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) $(SHLIBFLAGS) -o $(LIBDIR)/$(call DEAL_II_SHLIB_NAME,base.g) $(call DEAL_II_ADD_SONAME,base.g) $(go-files) $(extra-g.o-files)  $(deplibs.g)
	@ln -f -s $(call DEAL_II_SHLIB_NAME,base.g) $@


# Rule to generate the dependency file. This file is
# automagically remade whenever needed, i.e. whenever
# one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom
# of this file.
#
# If the command fails, then remove Makefile.dep again and fail
Makefile.dep: $(cc-files) $(h-files) Makefile $D/common/Make.global_options
	@echo ============================ Remaking base/Makefile
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B$(LIBDIR)/base" $(cc-files) \
		> $@ \
	 || (rm -f $@ ; false)



# include all the dependencies
include Makefile.dep
endif

clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep source/*.inst


.PHONY: clean
.PHONY: lib lib$(static-lib-suffix) lib.g$(static-lib-suffix)

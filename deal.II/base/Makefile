# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002



# deal.II top level directory
D=..


# create lists of file names
cc-files     = $(shell echo source/*.cc)
o-files-base = $(addprefix $(LIBDIR)/base/, $(patsubst source/%,%,$(cc-files:.cc=.$(OBJEXT))))
o-files      = $(sort $(o-files-base))
go-files     = $(o-files:.$(OBJEXT)=.g.$(OBJEXT))
h-files      = $(sort $(shell echo include/base/*.h))

# get options which hold for all files of the project
include $D/common/Make.global_options


# production rules
ifneq ($(with-multithreading),no)
  MT = MT
else
  MT = ==
endif

$(LIBDIR)/base/%.g.$(OBJEXT) :
	@echo "=====base=============debug======$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
$(LIBDIR)/base/%.$(OBJEXT) :
	@echo "=====base=============optimized==$(MT)== $(<F)"
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@


# rules
ifeq ($(enable-shared),yes)
  lib: $(LIBDIR)/libbase.g.so $(LIBDIR)/libbase.so
  libg: $(LIBDIR)/libbase.g.so
  libo: $(LIBDIR)/libbase.so
else
  lib: $(LIBDIR)/libbase.g.a $(LIBDIR)/libbase.a
  libg: $(LIBDIR)/libbase.g.a
  libo: $(LIBDIR)/libbase.a
endif

$(LIBDIR)/libbase.a: $(o-files)
	@echo "=====base=============optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files)
	@$(RANLIB) $@

$(LIBDIR)/libbase.g.a: $(go-files)
	@echo "=====base=============debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files)
	@$(RANLIB) $@


$(LIBDIR)/libbase.so: $(o-files)
	@echo "=====base=============optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files)

$(LIBDIR)/libbase.g.so: $(go-files)
	@echo "=====base=============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files)




clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(h-files) Makefile $D/common/Make.global_options
	@echo ============================ Remaking base/Makefile
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B$(LIBDIR)/base" $(cc-files) \
		> Makefile.dep


# include all the dependencies
include Makefile.dep

.PHONY: clean
.PHONY: lib lib.a lib.g.a

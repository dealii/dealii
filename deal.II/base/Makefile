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
lib: $(LIBDIR)/libbase.g$(lib-suffix) $(LIBDIR)/libbase$(lib-suffix)
libg: $(LIBDIR)/libbase.g$(lib-suffix)
libo: $(LIBDIR)/libbase$(lib-suffix)

$(LIBDIR)/libbase$(static-lib-suffix): $(o-files)
	@echo "=====base=============optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files)
	@$(RANLIB) $@

$(LIBDIR)/libbase.g$(static-lib-suffix): $(go-files)
	@echo "=====base=============debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files)
	@$(RANLIB) $@


$(LIBDIR)/libbase$(shared-lib-suffix): $(o-files)
	@echo "=====base=============optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files)

$(LIBDIR)/libbase.g$(shared-lib-suffix): $(go-files)
	@echo "=====base=============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files)




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
	@echo ============================ Remaking base/Makefile
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B$(LIBDIR)/base" $(cc-files) \
		> Makefile.dep
	@if test -s $@ ; then : else rm $@ ; fi



# include all the dependencies
include Makefile.dep

.PHONY: clean
.PHONY: lib lib$(static-lib-suffix) lib.g$(static-lib-suffix)

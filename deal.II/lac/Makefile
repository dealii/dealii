# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# deal.II top level directory
D=..


# create lists of file names. note that we need several output
# files for data_out_base.cc, which is why we use the several
# intermediate o-files-base* rules
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

$(LIBDIR)/liblac.a: $(o-files)
	@echo "=====lac==============optimized==$(MT)== Linking library:   $(@F)"
	@ar ru $@ $(o-files)

$(LIBDIR)/liblac.g.a: $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@ar ru $@ $(go-files)


$(LIBDIR)/liblac.so: $(o-files)
	@echo "=====lac==============optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files)

$(LIBDIR)/liblac.g.so: $(go-files)
	@echo "=====lac==============debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files)



clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(h-files) Makefile
	@echo ============================ Remaking lac/Makefile
	@$(PERL) $D/common/scripts/make_dependencies.pl $(INCLUDE) "-B\$$(LIBDIR)/lac" $(cc-files) \
		> Makefile.dep

# include all the dependencies
include Makefile.dep

.PHONY: clean
.PHONY: lib lib.a lib.g.a

# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# deal.II top level directory
D=..


# create lists of file names. note that we need several output
# files for data_out_base.cc, which is why we use the several
# intermediate o-files-base* rules
cc-files     = $(shell echo source/*.cc)
o-files-base1= $(addprefix $(LIBDIR)/lac/, $(patsubst source/%,%,$(cc-files:.cc=.o)))
o-files-base2= $(filter-out %data_out_base.o, $(o-files-base1))
o-files-base3= $(sort $(o-files-base2))
o-files      = $(o-files-base3)
go-files     = $(o-files-base3:.o=.g.o)
h-files      = $(sort $(shell echo include/lac/*.h))

forward-declarations = $D/lac/include/lac/forward_declarations.h

# get options which hold for all files of the project
include $D/common/Make.global_options


# production rules
ifneq ($(with-multithreading),no)
  MT = MT
else
  MT = ==
endif

$(LIBDIR)/lac/%.g.o :
	@echo =====lac==============debug======$(MT)== $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
$(LIBDIR)/lac/%.o :
	@echo =====lac==============optimized==$(MT)== $(<F)
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



$(LIBDIR)/liblac.a: $(forward-declarations) $(o-files)
	@echo =====lac==============optimized==$(MT)== Linking library:   $(@F)
	@ar ru $@ $(o-files)

$(LIBDIR)/liblac.g.a: $(forward-declarations) $(go-files)
	@echo =====lac==============debug======$(MT)== Linking library:   $(@F)
	@ar ru $@ $(go-files)


$(LIBDIR)/liblac.so: $(forward-declarations) $(o-files)
	@echo =====lac==============optimized==$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(o-files)

$(LIBDIR)/liblac.g.so: $(forward-declarations) $(go-files)
	@echo =====lac==============debug======$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(go-files)



clean:
	-rm -f *~ include/*~ include/*/*~ Makefile.dep


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
	@echo ============================ Remaking Makefile
	@$(PERL) $D/common/scripts/make_dependencies.pl $(INCLUDE) "-B\$$(LIBDIR)/lac" $(cc-files) \
		> Makefile.dep

# include all the dependencies
include Makefile.dep

# include a rule which tells how to remake the forward declarations file
include $D/common/scripts/forward_declarations


.PHONY: clean
.PHONY: lib lib.a lib.g.a

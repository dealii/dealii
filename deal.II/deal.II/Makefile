# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999


# deal.II top level directory
D=..


# get options which hold for all files of the project
include $D/common/Make.global_options



# top level targets
default: 2dg 2d
all: 2dg 2d 3dg 3d 1dg 1d


ifeq ($(enable-shared),yes)
  1dg: $(LIBDIR)/libdeal_II_1d.g.so
  1d : $(LIBDIR)/libdeal_II_1d.so
  2dg: $(LIBDIR)/libdeal_II_2d.g.so 
  2d : $(LIBDIR)/libdeal_II_2d.so
  3dg: $(LIBDIR)/libdeal_II_3d.g.so 
  3d : $(LIBDIR)/libdeal_II_3d.so
else
  1dg: $(LIBDIR)/libdeal_II_1d.g.a
  1d : $(LIBDIR)/libdeal_II_1d.a
  2dg: $(LIBDIR)/libdeal_II_2d.g.a
  2d : $(LIBDIR)/libdeal_II_2d.a
  3dg: $(LIBDIR)/libdeal_II_3d.g.a
  3d : $(LIBDIR)/libdeal_II_3d.a
endif


# create lists of file names
cc-files     = $(shell echo source/*/*.cc)
h-files      = $(sort $(shell echo include/*/*.h))

forward-declarations = $D/deal.II/include/grid/forward_declarations.h


# strip subdirectories from cc file names, preprend a "lib/[g]o"
# and change the suffix. We place the object files with or without
# debug info in different directories, since we want to use the
# -frepo switch to gcc; this, however, creates a .rpo file which
# would be the same for both debug and non-debug versions, if we
# don't care about different directories
tmp1        = $(notdir $(cc-files))
o-files-1d  = $(addprefix $(LIBDIR)/1d/, $(tmp1:.cc=_1d.o) )
go-files-1d = $(addprefix $(LIBDIR)/1d/, $(tmp1:.cc=_1d.g.o))
o-files-2d  = $(addprefix $(LIBDIR)/2d/, $(tmp1:.cc=_2d.o) )
go-files-2d = $(addprefix $(LIBDIR)/2d/, $(tmp1:.cc=_2d.g.o))
o-files-3d  = $(addprefix $(LIBDIR)/3d/, $(tmp1:.cc=_3d.o) )
go-files-3d = $(addprefix $(LIBDIR)/3d/, $(tmp1:.cc=_3d.g.o))

o-files  = $(o-files-1d)  $(o-files-2d)  $(o-files-3d)
go-files = $(go-files-1d) $(go-files-2d) $(go-files-3d)


# rules how to generate object files from source files. note that
# there are some files which needs exceptional rules; these are
# listed immediately below
ifneq ($(with-multithreading),no)
  MT = MT
else
  MT = ==
endif

$(LIBDIR)/1d/%.g.o :
	@echo =====deal.II====1d====debug======$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=1 -c $< -o $@
$(LIBDIR)/1d/%.o :
	@echo =====deal.II====1d====optimized==$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=1 -c $< -o $@

$(LIBDIR)/2d/%.g.o :
	@echo =====deal.II====2d====debug======$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=2 -c $< -o $@
$(LIBDIR)/2d/%.o :
	@echo =====deal.II====2d====optimized==$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=2 -c $< -o $@

$(LIBDIR)/3d/%.g.o :
	@echo =====deal.II====3d====debug======$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=3 -c $< -o $@
$(LIBDIR)/3d/%.o :
	@echo =====deal.II====3d====optimized==$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=3 -c $< -o $@




# rules how to make the libraries themselves
$(LIBDIR)/libdeal_II_1d.g.a: $(go-files-1d)
	@echo =====deal.II====1d====debug======$(MT)== Linking library:   $(@F)
	@ar ru $@ $(go-files-1d)
$(LIBDIR)/libdeal_II_1d.a: $(o-files-1d)
	@echo =====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)
	@ar ru $@ $(o-files-1d)


$(LIBDIR)/libdeal_II_2d.g.a: $(go-files-2d)
	@echo =====deal.II====2d====debug======$(MT)== Linking library:   $(@F)
	@ar ru $@ $(go-files-2d)
$(LIBDIR)/libdeal_II_2d.a: $(o-files-2d)
	@echo =====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)
	@ar ru $@ $(o-files-2d)


$(LIBDIR)/libdeal_II_3d.g.a: $(go-files-3d)
	@echo =====deal.II====3d====debug======$(MT)== Linking library:   $(@F)
	@ar ru $@ $(go-files-3d)
$(LIBDIR)/libdeal_II_3d.a: $(o-files-3d)
	@echo =====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)
	@ar ru $@ $(o-files-3d)


$(LIBDIR)/libdeal_II_1d.g.so: $(go-files-1d)
	@echo =====deal.II====1d====debug======$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(go-files-1d)

$(LIBDIR)/libdeal_II_1d.so: $(o-files-1d)
	@echo =====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(o-files-1d)


$(LIBDIR)/libdeal_II_2d.g.so: $(go-files-2d)
	@echo =====deal.II====2d====debug======$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(go-files-2d)

$(LIBDIR)/libdeal_II_2d.so: $(o-files-2d)
	@echo =====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(o-files-2d)


$(LIBDIR)/libdeal_II_3d.g.so: $(go-files-3d)
	@echo =====deal.II====3d====debug======$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(go-files-3d)

$(LIBDIR)/libdeal_II_3d.so: $(o-files-3d)
	@echo =====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)
	@$(CC) $(LDFLAGS) -shared -o $@ $(o-files-3d)





#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
	@echo ============================ Remaking Makefile
	@$(PERL) $D/common/scripts/make_dependencies.pl $(INCLUDE) "-B\$$(LIBDIR)" $(cc-files) \
		| $(PERL) -p -e 's!LIBDIR\)/(.*):!LIBDIR)/DIM_PLACEHOLDER/$$1:!g;' \
		| $(PERL) -pe 's!((\.g)?.o):!_DIM_PLACEHOLDER$$1:!g;' \
	        | $(PERL) -pe 's!^(.*)/DIM_PLACEHOLDER/(.*)_DIM_PLACEHOLDER(\..*):!$$1/1d/$$2_1d$$3 $$1/2d/$$2_2d$$3 $$1/3d/$$2_3d$$3:!g;'\
		> Makefile.dep

# include all the dependencies
include Makefile.dep

# include a rule which tells how to remake the forward declarations file
include $D/common/scripts/forward_declarations




################### some clean-up rules

clean: 
	-rm -f *~ */*~ */*/*~ Makefile.dep


###################

.PHONY: all deal.II TAGS
.PHONY: clean
.PHONY: default 1d 2d 3d all

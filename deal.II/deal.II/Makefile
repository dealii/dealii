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
o-files-1d  = $(addprefix lib/1d/o/,  $(tmp1:.cc=_1d.o) )
go-files-1d = $(addprefix lib/1d/go/, $(tmp1:.cc=_1d.go))
o-files-2d  = $(addprefix lib/2d/o/,  $(tmp1:.cc=_2d.o) )
go-files-2d = $(addprefix lib/2d/go/, $(tmp1:.cc=_2d.go))
o-files-3d  = $(addprefix lib/3d/o/,  $(tmp1:.cc=_3d.o) )
go-files-3d = $(addprefix lib/3d/go/, $(tmp1:.cc=_3d.go))

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

lib/1d/go/%.go :
	@echo =====deal.II====1d====debug======$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=1 -c $< -o $@
lib/1d/o/%.o :
	@echo =====deal.II====1d====optimized==$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=1 -c $< -o $@

lib/2d/go/%.go :
	@echo =====deal.II====2d====debug======$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=2 -c $< -o $@
lib/2d/o/%.o :
	@echo =====deal.II====2d====optimized==$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=2 -c $< -o $@

lib/3d/go/%.go :
	@echo =====deal.II====3d====debug======$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=3 -c $< -o $@
lib/3d/o/%.o :
	@echo =====deal.II====3d====optimized==$(MT)== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=3 -c $< -o $@




# rules how to make the libraries themselves
$(LIBDIR)/libdeal_II_1d.g.a: $(go-files-1d)
	@echo =====deal.II====1d====debug======$(MT)== Linking library:   $(@F)
	@ar ruv $@ $(go-files-1d)
$(LIBDIR)/libdeal_II_1d.a: $(o-files-1d)
	@echo =====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)
	@ar ruv $@ $(o-files-1d)


$(LIBDIR)/libdeal_II_2d.g.a: $(go-files-2d)
	@echo =====deal.II====2d====debug======$(MT)== Linking library:   $(@F)
	@ar ruv $@ $(go-files-2d)
$(LIBDIR)/libdeal_II_2d.a: $(o-files-2d)
	@echo =====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)
	@ar ruv $@ $(o-files-2d)


$(LIBDIR)/libdeal_II_3d.g.a: $(go-files-3d)
	@echo =====deal.II====3d====debug======$(MT)== Linking library:   $(@F)
	@ar ruv $@ $(go-files-3d)
$(LIBDIR)/libdeal_II_3d.a: $(o-files-3d)
	@echo =====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)
	@ar ruv $@ $(o-files-3d)


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
lib/Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
	@echo ============================ Remaking Makefile
	@$(PERL) $D/common/scripts/make_dependencies.pl $(INCLUDE) $(cc-files) \
		| $(PERL) -p -e 's!^lib/(g?o/.*):!lib/DIM_PLACEHOLDER/$$1:!g;' \
		| $(PERL) -pe 's!(\.g?o):!_DIM_PLACEHOLDER$$1:!g;' \
	        | $(PERL) -pe 's!^lib/DIM_PLACEHOLDER/(.*):!lib/1d/$$1 lib/2d/$$1 lib/3d/$$1:!g;'\
	        | $(PERL) -pe 's!(lib/1d/[^ ]+)_DIM_PLACEHOLDER!$$1_1d!g;' \
	        | $(PERL) -pe 's!(lib/2d/[^ ]+)_DIM_PLACEHOLDER!$$1_2d!g;' \
	        | $(PERL) -pe 's!(lib/3d/[^ ]+)_DIM_PLACEHOLDER!$$1_3d!g;' \
		> lib/Makefile.dep

# include all the dependencies
include lib/Makefile.dep

# include a rule which tells how to remake the forward declarations file
include $D/common/scripts/forward_declarations




################### some clean-up rules

clean: 
	-rm -f *~ */*~ */*/*~
	cd lib ; $(MAKE) clean


###################

.PHONY: all deal.II TAGS
.PHONY: clean lib-clean source-clean
.PHONY: default 1d 2d 3d all clean

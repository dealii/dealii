# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# deal.II top level directory
D=..


# get options which hold for all files of the project
include $D/common/Make.global_options



# top level targets
default: 2dg 2d
all: 2dg 2d 3dg 3d 1dg 1d


1dg: $(LIBDIR)/libdeal_II_1d.g$(lib-suffix)
1d : $(LIBDIR)/libdeal_II_1d$(lib-suffix)
2dg: $(LIBDIR)/libdeal_II_2d.g$(lib-suffix) 
2d : $(LIBDIR)/libdeal_II_2d$(lib-suffix)
3dg: $(LIBDIR)/libdeal_II_3d.g$(lib-suffix) 
3d : $(LIBDIR)/libdeal_II_3d$(lib-suffix)

# create lists of file names
cc-files     = $(shell echo source/*/*.cc)
h-files      = $(sort $(shell echo include/*/*.h))


# strip subdirectories from cc file names, preprend a "lib/[g]o"
# and change the suffix. We place the object files with or without
# debug info in different directories, since we want to use the
# -frepo switch to gcc; this, however, creates a .rpo file which
# would be the same for both debug and non-debug versions, if we
# don't care about different directories
tmp1        = $(notdir $(cc-files))
o-files-1d  = $(addprefix $(LIBDIR)/1d/, $(tmp1:.cc=_1d.$(OBJEXT)) )
go-files-1d = $(addprefix $(LIBDIR)/1d/, $(tmp1:.cc=_1d.g.$(OBJEXT)))
o-files-2d  = $(addprefix $(LIBDIR)/2d/, $(tmp1:.cc=_2d.$(OBJEXT)) )
go-files-2d = $(addprefix $(LIBDIR)/2d/, $(tmp1:.cc=_2d.g.$(OBJEXT)))
o-files-3d  = $(addprefix $(LIBDIR)/3d/, $(tmp1:.cc=_3d.$(OBJEXT)) )
go-files-3d = $(addprefix $(LIBDIR)/3d/, $(tmp1:.cc=_3d.g.$(OBJEXT)))

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

$(LIBDIR)/1d/%.g.$(OBJEXT) :
	@echo "=====deal.II====1d====debug======$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=1 -c $< -o $@
$(LIBDIR)/1d/%.$(OBJEXT) :
	@echo "=====deal.II====1d====optimized==$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=1 -c $< -o $@

$(LIBDIR)/2d/%.g.$(OBJEXT) :
	@echo "=====deal.II====2d====debug======$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=2 -c $< -o $@
$(LIBDIR)/2d/%.$(OBJEXT) :
	@echo "=====deal.II====2d====optimized==$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=2 -c $< -o $@

$(LIBDIR)/3d/%.g.$(OBJEXT) :
	@echo "=====deal.II====3d====debug======$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=3 -c $< -o $@
$(LIBDIR)/3d/%.$(OBJEXT) :
	@echo "=====deal.II====3d====optimized==$(MT)== $(<:source/%=%)"
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=3 -c $< -o $@



# in general, when linking shared libraries, we will want to link with
# the underlying libs as well. AIX even requires this. on the other
# hand, on DEC OSF with CXX this will fail with duplicate symbols, so
# make sure that it doesn't happen on that system and that compiler
ifneq ($(GXX-VERSION),compaq_cxx)
  deplibs.g = $(lib-lac.g) $(lib-base.g)
  deplibs.o = $(lib-lac.o) $(lib-base.o)
else
  deplibs.g =
  deplibs.o =
endif



# rules how to make the libraries themselves
$(LIBDIR)/libdeal_II_1d.g$(static-lib-suffix): $(go-files-1d)
	@echo "=====deal.II====1d====debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files-1d)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II_1d$(static-lib-suffix): $(o-files-1d)
	@echo "=====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files-1d)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II_2d.g$(static-lib-suffix): $(go-files-2d)
	@echo "=====deal.II====2d====debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files-2d)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II_2d$(static-lib-suffix): $(o-files-2d)
	@echo "=====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files-2d)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II_3d.g$(static-lib-suffix): $(go-files-3d)
	@echo "=====deal.II====3d====debug======$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(go-files-3d)
	@$(RANLIB) $@
$(LIBDIR)/libdeal_II_3d$(static-lib-suffix): $(o-files-3d)
	@echo "=====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)"
	@$(AR) ru $@ $(o-files-3d)
	@$(RANLIB) $@


$(LIBDIR)/libdeal_II_1d.g$(shared-lib-suffix): $(go-files-1d)
	@echo "=====deal.II====1d====debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files-1d) $(deplibs.g)

$(LIBDIR)/libdeal_II_1d$(shared-lib-suffix): $(o-files-1d)
	@echo "=====deal.II====1d====optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files-1d) $(deplibs.o)


$(LIBDIR)/libdeal_II_2d.g$(shared-lib-suffix): $(go-files-2d)
	@echo "=====deal.II====2d====debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files-2d) $(deplibs.g)

$(LIBDIR)/libdeal_II_2d$(shared-lib-suffix): $(o-files-2d)
	@echo "=====deal.II====2d====optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files-2d) $(deplibs.o)


$(LIBDIR)/libdeal_II_3d.g$(shared-lib-suffix): $(go-files-3d)
	@echo "=====deal.II====3d====debug======$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(go-files-3d) $(deplibs.g)

$(LIBDIR)/libdeal_II_3d$(shared-lib-suffix): $(o-files-3d)
	@echo "=====deal.II====3d====optimized==$(MT)== Linking library:   $(@F)"
	@$(SHLIBLD) $(LDFLAGS) -shared -o $@ $(o-files-3d) $(deplibs.o)





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
	@echo "============================ Remaking deal.II/Makefile"
	@$D/common/scripts/make_dependencies $(INCLUDE) "-B\$$(LIBDIR)" $(cc-files) \
		| $(PERL) -p -e 's!LIBDIR\)/(.*):!LIBDIR)/DIM_PLACEHOLDER/$$1:!g;' \
		| $(PERL) -pe 's!((\.g)?.$(OBJEXT)):!_DIM_PLACEHOLDER$$1:!g;' \
	        | $(PERL) -pe 's!^(.*)/DIM_PLACEHOLDER/(.*)_DIM_PLACEHOLDER(\..*):!$$1/1d/$$2_1d$$3 $$1/2d/$$2_2d$$3 $$1/3d/$$2_3d$$3:!g;'\
		> Makefile.dep
	@if test -s $@ ; then : else rm $@ ; fi


# include all the dependencies
include Makefile.dep




################### some clean-up rules

clean: 
	-rm -f *~ */*~ */*/*~ Makefile.dep


###################

.PHONY: all deal.II TAGS
.PHONY: clean
.PHONY: default 1d 2d 3d all

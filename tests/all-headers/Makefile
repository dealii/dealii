# collect all header files, and for each try to do the following:
# generate a file tmp.cc that contains nothing but "#include <xxx>"
# with this particular header file. then try to compile it. this 
# sometimes fails since this header does not include all other headers
# it actually needs. we don't usually see this since we may be 
# including these missing headers into our files before, so that nobody
# notices the problem. expose it this way.

include ../Makefile.paths
include $D/common/Make.global_options

############################################################
# First how to create executables, including all necessary
# flags:
############################################################

flags     = $(CXXFLAGS.g)

ifeq ($(findstring gcc,$(GXX_VERSION)),gcc)
flags += -Wno-missing-noreturn
endif

%.g.$(OBJEXT) : %.cc Makefile
	@echo =====debug========= $<
	@$(CXX) $(flags) -c $< -o $@

%.$(OBJEXT) : %.cc Makefile
	@echo =====optimized===== $<
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@


default:all

# collect all headers
HEADERS = $(shell cd $D; echo */include/*/*.h)
TESTS   = $(HEADERS:%=test-%)


# we need two auxiliary functions. the first one takes a filename and replaces
# all occurrences of slashes by hyphens, so that we can create a .cc file that
# corresponds to a header file in a different directory
subst_slashes = $(shell echo $(1) | $(PERL) -pi -e 's/\//-/g;')

# the second one takes a path to a header file and makes sure that we can
# #include it, i.e. it drops the top-level xxx/include/ part of the filename,
# since we only find include directives like <grid/tria.h>
normalize = $(shell echo $(1) | $(PERL) -pi -e 's/.*include.//g;')


# next the target that creates the target rules. First, we list all the header
# files as the head of the rule, followed by a colon. Then the actual rule
# that includes writing the .cc file that only includes a single header file,
# compiles it, and then deletes everything again
Makefile.dep: $(HEADERS:%=$D/%) Makefile
	@for i in $(HEADERS) ; do \
		echo test-$$i \\ ; \
	 done > Makefile.dep
	@echo : >> Makefile.dep
	@echo '	@echo "#include <$$(call normalize,$$@)>" > $$(call subst_slashes,$$@).cc' >> Makefile.dep
	@echo '	@$$(MAKE) $$(call subst_slashes,$$@).g.o' >> Makefile.dep
	@echo '	@rm -f $$(call subst_slashes,$$@).g.o' >> Makefile.dep
	@echo '	@rm -f $$(call subst_slashes,$$@).cc' >> Makefile.dep

include Makefile.dep


# here's what's to be done
all: $(TESTS)


# create a report, where we get exactly one line per test
report:
	@for test in $(TESTS) ; do \
		testname=all-headers/`echo $$test | perl -pi -e 's/.*include.//g; s/\//-/g;'` ; \
		if (($(MAKE) $$test 2>&1) > /dev/null); then \
			echo `date -u +"%Y-%m-%d %H:%M"` + $$testname ; \
		else \
			echo `date -u +"%Y-%m-%d %H:%M"` - $$testname ; \
		fi ; \
	done


clean:
	-rm -f Makefile.dep *.$(OBJEXT) \
               *.g.$(OBJEXT) *.output *~


.PHONY: report clean distclean

# $Id$
# Copyright W. Bangerth, 2003


# deal.II top level directory
D=../..

include $D/common/Make.global_options

# if compaq cxx and we have shared, then remove -shared again from
# LDFLAGS, since this prevents generation of a working executable. this
# is clearly a hack, but I'm too lazy for now to work things out
ifeq ($(GXX-VERSION),compaq_cxx)
  XLDFLAGS = $(LDFLAGS:-shared=)
endif

# rules to generate the program. let is be dependent on the global flags for
# the following reason: if you go to a different system, but keep with the
# same directory, you get into trouble because make_dependencies is not
# deleted by "make clean"; following this, we generate empty Makefile.dep's
# because we can't run the executable. So better be careful and regenerate
# make_dependencies whenever we think it may be beneficial
make_dependencies : make_dependencies.cc $D/common/Make.global_options
	@echo ============================ Compiling $@
	@$(CXX) $(CXXFLAGS.o) $(LDFLAGS) $< -o $@

clean:
	-rm -f make_dependencies

.PHONY: clean

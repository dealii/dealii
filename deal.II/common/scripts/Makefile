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

make_dependencies : make_dependencies.cc
	@echo ============================ Compiling $@
	@$(CXX) $(CXXFLAGS.o) -c $<
	@$(CXX) $(XLDFLAGS) -o $@ make_dependencies.o

clean:
	-rm -f make_dependencies

.PHONY: clean

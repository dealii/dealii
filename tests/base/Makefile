# $Id$
#
# Set root directory of DEAL distribution.
# This is only one letter to make the rest MORE readable.

D=../..

include $D/deal.II/Make.global_options

%.go : %.cc
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@
%.o : %.cc
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@


vpath %.a: $D/base/lib
vpath %.a: $D/lac/lib
vpath %.a: $D/deal.II/lib

test: run-reference run-logtest run-tensor


reference: reference.go -lbase.g
	$(CXX) $(CXXFLAGS) -o $@ $^

run-reference: reference
	-./reference

logtest: logtest.go -lbase.g
	$(CXX) $(CXXFLAGS) -o $@ $^

run-logtest: logtest
	-./logtest

tensor: tensor.go -lbase.g -llac.g
	$(CXX) $(CXXFLAGS) -o $@ $^

run-tensor: tensor
	-./tensor


clean:
	-rm *.o *.go reference logtest tensor



.PHONY: test run-reference run-logtest run-tensor clean

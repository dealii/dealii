# $Id$
#
# Set root directory of DEAL distribution.
# This is only one letter to make the rest MORE readable.

D=../..

CXX       = c++
INCLUDE   = -I$D/deal.II/include -I$D/lac/include \
            -I$D/base/include

CXXFLAGS  = -DDEBUG -g -Wall -W -pedantic -Wconversion \
            -Winline -Woverloaded-virtual -fno-rtti  \
            $(INCLUDE) -Ddeal_II_dimension=$(deal_II_dimension)


%.go : %.cc
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@
%.o : %.cc
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@


vpath %.a: $D/base/lib
vpath %.a: $D/lac/lib
vpath %.a: $D/deal.II/lib

ifeq ($(shell uname),Linux)
CXX       = /home/wolf/bin/gcc/bin/c++
endif


test: run-reference run-logtest run-tensor


reference: reference.go
	$(CXX) $(CXXFLAGS) -o $@ $^ $D/base/lib/libbase.g.a

run-reference: reference
	./reference

logtest: logtest.go
	$(CXX) $(CXXFLAGS) -o $@ $^ $D/base/lib/libbase.g.a

run-logtest: logtest
	./logtest

tensor: tensor.go
	$(CXX) $(CXXFLAGS) -o $@ $^ $D/base/lib/libbase.g.a

run-tensor: tensor
	./tensor


clean:
	-rm *.o *.go reference logtest tensor



.PHONY: test run-reference run-logtest run-tensor clean

# $Id$

#this file uses many small rules to allow
#   make -j
#on multiprocessor machines. Feel free to use this feature,
#it gives you quicker compile time and annoys those FEATFLOW
#people ;-)

include ../Make.global_options


cc-files = $(wildcard */*.cc)
h-files  = $(wildcard ../include/*/*.h)

# strip subdirectories from cc file names, preprend a "../lib/"
# and change the suffix:
tmp1     = $(notdir $(cc-files))
tmp2     = $(addprefix ../lib/, $(tmp1))
o-files  = $(tmp2:.cc=.o)
go-files = $(tmp2:.cc=.go)


deal_II_dimension=2

CXX       = c++
INCLUDE   = -I../include -I../../lac/include \
            -I../../base/include -I../../mia/include
CXXFLAGS.g= -DDEBUG -g -Wall -W -pedantic -Wconversion \
            -Winline -Woverloaded-virtual -fno-rtti  \
            $(INCLUDE) -Ddeal_II_dimension=$(deal_II_dimension)
CXXFLAGS  =-O3 -Wuninitialized -finline-functions -ffast-math \
           -funroll-loops -felide-constructors -fnonnull-objects \
           -fno-rtti $(INCLUDE) \
           -Ddeal_II_dimension=$(deal_II_dimension)

ifeq ($(shell uname),Linux)
CXX       = /home/wolf/bin/gcc/bin/c++
endif

ifeq ($(shell uname),SunOS)
CXX       = /usr/local/gcc/egcs-19980803/bin/c++
endif


%.go :
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o :
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@






all: ../lib/libdeal_II.g.a ../lib/libdeal_II.a

../lib/libdeal_II.g.a: $(go-files)
	@echo ============================ Updating library:   $@
	@ar ruv $@ $(go-files)

../lib/libdeal_II.a: $(o-files)
	@echo ============================ Updating library:   $@
	@ar ruv $@ $(o-files)


clean:
	@cd basic   ; rm -f *.o *.go *~
	@cd grid    ; rm -f *.o *.go *~
	@cd fe      ; rm -f *.o *.go *~
	@cd numerics; rm -f *.o *.go *~



.PHONY: all clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#use perl to generate rules for the .go files as well
#as to make rules not for tria.o and the like, but
#rather for libgrid.a(tria.o)
#
#Later addition: the perl script was made to compile
#right into a library. Since this is not what we want
#anymore, but since 1) I am too lazy to change the
#script and 2) others may use the script as it was,
#we don't change the script but rather pipe its
#result through another perl script to get what we
#want.
../lib/Makefile.dep: $(cc-files) $(h-files) Makefile ../Make.global_options
	@echo ============================ Remaking Makefile
	@perl ../Make_dep.pl ../lib/libgrid $(INCLUDE) $(cc-files) \
		| perl -p -e 's!^.*\(.*/(.*)\):!../lib/$$1:!g;' \
		> ../lib/Makefile.dep


include ../lib/Makefile.dep

# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1999

# Template for makefiles for the examples subdirectory. In principle,
# everything should be done automatically if you set the target file
# here correctly. We get deduce it from the files in the present
# directory:
target   = $(basename $(shell echo step-*.cc))

# All dependencies between files should be updated by the included
# file Makefile.dep if necessary. Object files are compiled into
# the archives ./Obj.a and ./Obj.g.a. By default, the debug version
# is used to link. It you don't like that, change the following
# variable to "off"
debug-mode = on



###############################################################################
# Internals

#deal include base path
D = ../../../..

include $D/common/Make.global_options

# get lists of files we need


# list of libraries needed to link with
libs.g   = $(lib-deal2-2d.g) \
	   $(lib-lac.g)      \
           $(lib-base.g)
libs     = $(lib-deal2-2d.o) \
	   $(lib-lac.o)      \
           $(lib-base.o)


# check whether we use debug mode or not
ifeq ($(debug-mode),on)
  libraries = $(target).go $(libs.g)
  flags     = $(CXXFLAGS.g)
else
  libraries = $(target).go $(libs)
  flags     = $(CXXFLAGS.o)
endif



# make rule for the target. $^ is the object file $(target).g?o
$(target) : $(libraries)
	@echo ============================ Linking $@
	@$(CXX) $(flags) -o $@ $^

# rule how to run the program
run: $(target)
	@echo ============================ Running $<
	@./$(target)


# rule to make object files
%.go : %.cc
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	-rm -f *.o *.go *~ Makefile.dep $(target) *gmv *gnuplot *gpl *eps



.PHONY: clean


# Rule to generate the dependency file. This file is
# automagically remade whenever needed, i.e. whenever
# one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom
# of this file.
#
# Since the script prefixes the output names by lib/g?o, we have to
# strip that again (the script was written for the main libraries and
# large projects where object files are put into subdirs)
Makefile.dep: $(target).cc Makefile \
              $(shell echo $D/base/include/base/*.h \
                           $D/lac/include/lac/*.h  \
                           $D/deal.II/include/*/*.h)
	@echo ============================ Remaking Makefile
	@perl $D/common/scripts/make_dependencies.pl  $(INCLUDE) $(target).cc \
		| perl -pi -e 's!lib/g?o/!!g;' \
		> Makefile.dep


include Makefile.dep


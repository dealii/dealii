############################################################
# Directory with the DEAL distribution
############################################################

D = ../../deal.II


############################################################
# Include general settings for including DEAL libraries
############################################################

include $D/common/Make.global_options


############################################################
# Set debug-mode as a default
############################################################

debug-mode = on


############################################################
# Define library names
############################################################

libs   = -ldeal_II_2d -llac -lbase
libs.g = -ldeal_II_2d.g -llac.g -lbase.g


############################################################
# Select compiler flags according to debug-mode
############################################################

ifeq ($(debug-mode),on)
libraries = $(libs.g)
flags     = $(CXXFLAGS.g) $(CXXFLAGS)
endif

ifeq ($(debug-mode),off)
libraries = $(libs)
flags     = $(CXXFLAGS.o) $(CXXFLAGS)
endif

%.go : %.cc Makefile
	@echo =====debug========= $<
	@$(CXX) $(flags) -c $< -o $@
%.o : %.cc Makefile
	@echo =====optimized===== $<
	@$(CXX) $(flags) -c $< -o $@


all: check

############################################################
# Typical block for building a running program
#
# 1. provide a list of source files in ...-cc-files
#
# 2. generate the list of object files according to debug-mode
#
# 3. make executable
#
# 4. Explicit dependencies of object files (will be automatic soon)
#
############################################################

show_shapes-cc-files = show_shapes.cc

ifeq ($(debug-mode),on)
show_shapes-o-files = $(show_shapes-cc-files:.cc=.go)
else
show_shapes-o-files = $(show_shapes-cc-files:.cc=.o)
endif

show_shapes: $(show_shapes-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^

############################################################

show_transform-cc-files = show_transform.cc

ifeq ($(debug-mode),on)
show_transform-o-files = $(show_transform-cc-files:.cc=.go)
else
show_transform-o-files = $(show_transform-cc-files:.cc=.o)
endif

show_transform: $(show_transform-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^

############################################################
# Continue with other targets if needed
############################################################

check: show_transform show_shapes
	@echo Running programs
	@show_transform
	@show_shapes
	@echo Checking results
	$(MAKE) dodiff

%.check: %.dat
	diff $@ $<

dat = $(wildcard *.dat)
check = $(dat:.dat=.check)

dodiff: $(check)
	@echo OK

############################################################
# Cleanup targets
############################################################

clean:
	rm -f *.o *.go

veryclean: clean
	rm -f show_shapes show_transform *.inp *.gpl *.eps *.gnuplot *.dat

############################################################
# Automatic generation of dependencies
############################################################

all-cc-files = $(show_shapes-cc-files) $(show_transform-cc-files)

Make.depend: $(all-cc-files)
	@echo =====Dependecies=== Make.depend
	@$(CXX) $(flags) $^ -M > $@
	@perl -pi~ -e 's/(^[^.]+)\.o:/\1.o \1.go:/;' $@

include Make.depend

.PHONY: clean veryclean check

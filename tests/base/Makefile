############################################################
# Directory with the DEAL distribution
############################################################

D = ../..


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

libs   = -llac -lbase
libs.g = -llac.g -lbase.g


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


all: logtest.check reference.check tensor.check
exe: $(all:.check=.testcase)
output: $(all:.check=.output)
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

logtest-cc-files = logtest.cc

ifeq ($(debug-mode),on)
logtest-o-files = $(logtest-cc-files:.cc=.go)
else
logtest-o-files = $(logtest-cc-files:.cc=.o)
endif

logtest.testcase: $(logtest-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^


############################################################


tensor-cc-files = tensor.cc

ifeq ($(debug-mode),on)
tensor-o-files = $(tensor-cc-files:.cc=.go)
else
tensor-o-files = $(tensor-cc-files:.cc=.o)
endif

tensor.testcase: $(tensor-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^



############################################################


reference-cc-files = reference.cc

ifeq ($(debug-mode),on)
reference-o-files = $(reference-cc-files:.cc=.go)
else
reference-o-files = $(reference-cc-files:.cc=.o)
endif

reference.testcase: $(reference-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^



############################################################
# Continue with other targets if needed
############################################################


############################################################
# Continue with other targets if needed
############################################################


target1-cc-files = t1.cc t2.cc t3.cc

ifeq ($(debug-mode),on)
target1-o-files = $(target1-cc-files:.cc=.go)
else
target1-o-files = $(target1-cc-files:.cc=.o)
endif

target1: $(target1-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^


############################################################
# Postprocessing
############################################################

%.output:%.testcase
	@echo =====Running======= $<
	@$<
	@perl -pi.bak -e 's/JobId.*//;s/value.*//;' $@
	@rm $@.bak

%.check:%.output
	@-diff $< $(patsubst %.output,%.checked, $<) && echo '=====OK============'
############################################################
# Cleanup targets
############################################################

clean:
	rm -f *.o *.go *.output

veryclean: clean
	rm -f *.testcase *.inp *.gpl *.eps *.gnuplot

############################################################
# Automatic generation of dependencies
############################################################

all-cc-files = $(wildcard *.cc)

Make.depend: $(all-cc-files)
	@echo =====Dependecies=== Make.depend
	@$(CXX) $(flags) $^ -M > $@
	@perl -pi~ -e 's/(^[^.]+)\.o:/\1.o \1.go:/;' $@

include Make.depend

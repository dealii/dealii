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


all: solver.check mgbase.check mg.check
exe: $(all:.check=.testcase) benchmark
run: $(all:.check=.output)
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

solver-cc-files = solver.cc testmatrix.cc

ifeq ($(debug-mode),on)
solver-o-files = $(solver-cc-files:.cc=.go)
else
solver-o-files = $(solver-cc-files:.cc=.o)
endif

solver.testcase: $(solver-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^


############################################################


mgbase-cc-files = mgbase.cc

ifeq ($(debug-mode),on)
mgbase-o-files = $(mgbase-cc-files:.cc=.go)
else
mgbase-o-files = $(mgbase-cc-files:.cc=.o)
endif

mgbase.testcase: $(mgbase-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^



############################################################


mg-cc-files = mg.cc testmatrix.cc

ifeq ($(debug-mode),on)
mg-o-files = $(mg-cc-files:.cc=.go)
else
mg-o-files = $(mg-cc-files:.cc=.o)
endif

mg.testcase: $(mg-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^



############################################################
# Continue with other targets if needed
############################################################


benchmark-cc-files = benchmark.cc

ifeq ($(debug-mode),on)
benchmark-o-files = $(benchmark-cc-files:.cc=.go)
else
benchmark-o-files = $(benchmark-cc-files:.cc=.o)
endif

benchmark.testcase: $(benchmark-o-files) $(libraries)
	$(CXX) $(flags) -o $@ $^

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
	@diff $< $(patsubst %.output,%.checked, $<) && echo '=====OK============'
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

all-cc-files = $(shell echo *.cc)

Make.depend: $(all-cc-files)
	@echo =====Dependencies=== Make.depend
	@$(CXX) $(flags) $^ -M > $@
	@perl -pi~ -e 's/(^[^.]+)\.o:/\1.o \1.go:/;' $@

include Make.depend

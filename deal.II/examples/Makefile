# $Id$
#
# This Makefile only recurses into the subdirs
include ../common/Make.global_options

all: default

# existing examples. take dirnames and strip 'step'
steps = $(shell echo step-? step-??)


# default is: build all examples. for each example, there is a target
# build-step-N, where N in [1...]
default: $(addprefix build-,$(steps)) doxygen

# run example programs; make a target run-step-N for each N
run: $(addprefix run-,$(steps)) run-doxygen

# clean subdirs; make a target clean-step-N for each N
clean:
	-rm -f *~ */*~ */*.$(OBJEXT) */*gmv */*gnuplot */*gpl */*eps */*pov */*tex \
	$(addsuffix $(EXEEXT),$(basename $(shell echo step-*/step-*.cc))) \
	doxygen/*.exe */Makefile.dep


# for each build/run/clean target: strip the build- prefix of the
# target and build in that directory
build-step-%:
	cd $(@:build-%=%) ; $(MAKE)
run-step-%:
	cd $(@:run-%=%) ; $(MAKE) run
clean-step-%:
	cd $(@:clean-%=%) ; $(MAKE) clean
doxygen:
	cd doxygen ; $(MAKE)
run-doxygen:
	cd doxygen ; $(MAKE) run
clean-doxygen:
	cd doxygen ; $(MAKE) clean

# all targets in this directory do not produce files, so they are
# .PHONY:
.PHONY: $(addprefix build-step-,$(steps)) \
        $(addprefix run-step-,$(steps))   \
        $(addprefix clean-step-,$(steps)) \
	doxygen doxygen-run doxygen-clean
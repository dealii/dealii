# $Id$
# Copyright W. Bangerth, University of Heidelberg, 2002

D = ..

# get options which hold for all files of the project
include $D/common/Make.global_options


default:
	@echo There is no default target in this directory


# clean everything
clean: clean-objects

# only clean object files, not libraries and executables
clean-objects:  clean-base clean-lac \
		clean-1d clean-2d clean-3d \
		clean-contrib


# clean targets for the individual subdirectories. remove
# object files and template instantiation files.
# extract dir name from target
clean-base clean-lac clean-1d clean-2d clean-3d:
	-cd $(@:clean-%=%) ; rm *.$(OBJEXT) *.ti *.ii

clean-contrib:
	-rm contrib/*/*.$(OBJEXT)

distclean: clean
	-rm -f lib* bin/*

.PHONY: clean clean-objects 
.PHONY: clean-base clean-lac clean-1d clean-2d clean-3d
.PHONY: clean-contrib

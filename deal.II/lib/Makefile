# $Id$
# Copyright W. Bangerth, University of Heidelberg, 2002


default:
	@echo There is no default target in this directory


# clean everything
clean: clean-objects

# only clean object files, not libraries
clean-objects:  clean-base clean-lac \
		clean-1d clean-2d clean-3d \
		clean-contrib


# clean targets for the individual subdirectories. remove
# object files and template instantiation files.
# extract dir name from target
clean-base clean-lac clean-1d clean-2d clean-3d:
	-cd $(@:clean-%=%) ; rm *.o *.ti *.ii

clean-contrib:
	-rm contrib/*/*.o

veryclean: clean
	-rm -f lib*

.PHONY: clean clean-objects 
.PHONY: clean-base clean-lac clean-1d clean-2d clean-3d
.PHONY: clean-contrib

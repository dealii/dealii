# $Id$
# Copyright W. Bangerth, University of Heidelberg, 2002


default:
	@echo There is no default target in this directory


# clean everything
clean: clean-objects
	-rm lib*

# only clean object files, not libraries
clean-objects: clean-base clean-lac clean-1d clean-2d clean-3d


# clean targets for the individual subdirectories.
# extract dir name from target
clean-%:
	-cd $(@:clean-%=%) ; rm *.o


.PHONY: clean clean-objects clean-base clean-lac clean-1d clean-2d clean-3d

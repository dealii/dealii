# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998

# note that this file does not declare anything apart from
# clean-up rules. to compile the library, use the make file
# in the source subdir to re-make the libraries


default:
	@echo "------------------------------------------"
	@echo "There is no default rule in this Makefile."
	@echo "------------------------------------------"

# remove all object files, libraries and temporary (st*) files
clean: clean_1d clean_2d clean_3d
	@echo "Cleaning present directory..."
	-@rm -f st* lib*

clean_1d:
	@echo "Cleaning 1d directory..."
	-@cd 1d ; rm -f o/*.o go/*.go

clean_2d:
	@echo "Cleaning 2d directory..."
	-@cd 2d ; rm -f o/*.o go/*.go

clean_3d:
	@echo "Cleaning 3d directory..."
	-@cd 3d ; rm -f o/*.o go/*.go



.PHONY: clean_1d clean_2d clean_3d

# $Id$

structure.png: structure.dot
	@echo ================== Making $@
	@dot -Tpng -Tcmapx -O $<
	@mv $<.png $(basename $<).png
	@mv $<.cmapx $(basename $<).cmapx

###########################################################################
# $Id$
###########################################################################

all: lapack_templates.h
	cp $< ../../lac/include/lac/

lapack_templates.h: templates.pl lapack_templates.h.in
	perl  $^ > $@

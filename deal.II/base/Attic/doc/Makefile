# $Id$


KDOCFLAGS   = -I../../../deal.II/doc/kdoc ../../../deal.II/doc/kdoc/kdoc -a -p
kdoc.inc  = $(wildcard ../include/base/*.h)



doc-html: kdoc cvslog


# make kdoc doc; make sure that the *.kdoc files exist by
# using the dependancies and the following rules.
kdoc: $(kdoc.inc)
	cd doc.kdoc ; perl $(KDOCFLAGS) -ubase   -dbase   \
		base   $(kdoc.inc:..%=../..%)

cvslog:
	@cd .. ; ../deal.II/doc/cvslog/cvs2html -o doc/cvslog/base #-a -k

clean:
	rm -f 	doc.kdoc/base/* \
		cvslog/* \
		*~


.PHONY: doc-html kdoc cvslog clean

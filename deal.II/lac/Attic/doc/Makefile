# $Id$


KDOCFLAGS   = -I../kdoc ../../../deal.II/doc/kdoc/kdoc -a -p
kdoc.inc  = $(wildcard ../include/lac/*.h)



doc-html: kdoc cvslog


# make kdoc doc; make sure that the *.kdoc files exist by
# using the dependancies and the following rules.
kdoc: $(kdoc.inc)
	cd doc.kdoc ; perl $(KDOCFLAGS) -ulac   -dlac   \
		lac   $(kdoc.inc:..%=../..%)

cvslog:
	@cd .. ; ../deal.II/doc/cvslog/cvs2html -o doc/cvslog/lac -a -k

clean:
	rm -f 	doc.kdoc/lac/* \
		cvslog/* \
		*~


.PHONY: doc-html kdoc cvslog clean

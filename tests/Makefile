# $Id$
# Common Makefile for all test directories

all:
	cd base ; $(MAKE)
	cd lac ; $(MAKE)
	cd fe ; $(MAKE)
	cd deal.II ; $(MAKE)


clean: clean-base clean-lac clean-deal.II clean-fe

clean-base:
	cd base ; $(MAKE) clean

clean-lac:
	cd lac ; $(MAKE) clean

clean-fe:
	cd fe ; $(MAKE) clean

clean-deal.II:
	cd deal.II ; $(MAKE) clean


.PHONY: all clean clean-base clean-lac clean-fe clean-deal.II

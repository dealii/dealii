# $Id$
# Copyright (C) 2005 by Wolfgang Bangerth


D = ../..

include $D/common/Make.global_options

ifneq ($(with-multithreading),no)
  MT = MT
else
  MT = ==
endif



# rules
$(LIBDIR)/contrib/functionparser/fparser.$(OBJEXT): fparser.h fparser.cc
	@echo =====function parser==optimized==$(MT)== $(<F)
	@$(CXX) $(CXXFLAGS.o) -c fparser.cc -o $@



clean:
	-rm -f *~ 


.PHONY: clean

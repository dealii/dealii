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
$(LIBDIR)/contrib/functionparser/fparser.$(OBJEXT): fparser.cc fparser.h
	@echo =====function parser==optimized==$(MT)== $(<F)
	@$(CXX) $(CXXFLAGS.o) -c fparser.cc -o $@



clean:
	-rm -f *~ 


.PHONY: clean

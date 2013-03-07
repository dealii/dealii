#------------------------------------------------------------------------------
# GNUmakefile for the AMD MATLAB mexFunction
#------------------------------------------------------------------------------

all: amd

include ../Make/Make.include

MX = $(MEX) -I../Include

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid

INC = ../Include/amd.h ../Source/amd_internal.h

MEXAMD = $(addsuffix .o, $(subst amd_,amd_m_,$(AMD)))

amd_m_%.o: ../Source/amd_%.c $(INC)
	$(MX) -DDINT -c $<
	- $(MV) amd_$*.o $@

amd: amd_mex.c $(MEXAMD) $(INC)
	$(MX) -output amd amd_mex.c $(MEXAMD)

#------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#------------------------------------------------------------------------------

purge: clean
	- $(RM) amd.mex* amd.dll

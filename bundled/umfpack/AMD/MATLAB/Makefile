#------------------------------------------------------------------------------
# compile the AMD mexFunction for MATLAB (original make only)
#------------------------------------------------------------------------------

# This is a very ugly Makefile, and is only provided for those who do not
# have GNU make.  Note that it is not used if you have GNU make.  It ignores
# dependency checking and just compiles everything.  It was created
# automatically, via make -n using the GNUmakefile.  That way, I don't have
# maintain two Makefiles.

all: amd

include ../Make/Make.include

MX = $(MEX) -I../Include

amd:
	$(MX) -DDINT -c ../Source/amd_aat.c
	$(MV) amd_aat.o amd_m_aat.o
	$(MX) -DDINT -c ../Source/amd_1.c
	$(MV) amd_1.o amd_m_1.o
	$(MX) -DDINT -c ../Source/amd_2.c
	$(MV) amd_2.o amd_m_2.o
	$(MX) -DDINT -c ../Source/amd_dump.c
	$(MV) amd_dump.o amd_m_dump.o
	$(MX) -DDINT -c ../Source/amd_postorder.c
	$(MV) amd_postorder.o amd_m_postorder.o
	$(MX) -DDINT -c ../Source/amd_post_tree.c
	$(MV) amd_post_tree.o amd_m_post_tree.o
	$(MX) -DDINT -c ../Source/amd_defaults.c
	$(MV) amd_defaults.o amd_m_defaults.o
	$(MX) -DDINT -c ../Source/amd_order.c
	$(MV) amd_order.o amd_m_order.o
	$(MX) -DDINT -c ../Source/amd_control.c
	$(MV) amd_control.o amd_m_control.o
	$(MX) -DDINT -c ../Source/amd_info.c
	$(MV) amd_info.o amd_m_info.o
	$(MX) -DDINT -c ../Source/amd_valid.c
	$(MV) amd_valid.o amd_m_valid.o
	$(MX) -output amd amd_mex.c amd_m_aat.o \
	    amd_m_1.o amd_m_2.o amd_m_dump.o amd_m_postorder.o \
	    amd_m_post_tree.o amd_m_defaults.o amd_m_order.o amd_m_control.o \
	    amd_m_info.o amd_m_valid.o

#------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#------------------------------------------------------------------------------

purge: clean
	- $(RM) amd.mex* amd.dll

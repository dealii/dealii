//----------------------------  multithread_info.cc  ----------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  multithread_info.cc  ----------------


#include <base/multithread_info.h>

#if defined(__linux__)
#  include <fstream>
#  include <string>
#endif

#if defined(__sun__) || defined(__osf__) || defined(_AIX)
#  include <unistd.h>
#endif

#if defined(__sgi__)
#  include <unistd.h>
#endif

#if defined(__MACH__) && defined(__ppc__) && defined(__APPLE__)
#  include <sys/types.h>
#  include <sys/sysctl.h>
#endif

#if DEAL_II_USE_MT == 1

/* Detecting how many processors a given machine has is something that
   varies greatly between operating systems. For a few operating
   systems, we have figured out how to do that below, but some others
   are still missing. If you find a way to do this on your favorite
   system, please let us know.
 */


#  if defined(__linux__)

unsigned int MultithreadInfo::get_n_cpus()
{
				   // on linux, we simply count the
				   // number of lines listing
				   // individual processors when
				   // reading from /proc/cpuinfo
  std::ifstream cpuinfo;
  std::string search;
  unsigned int nCPU = 0;
  
  cpuinfo.open("/proc/cpuinfo");

  AssertThrow(cpuinfo,ExcProcNotPresent());
  
  while(cpuinfo)
    {
      cpuinfo >> search;
      if (search.find("processor")!=std::string::npos)
	nCPU++;	  
    }
  cpuinfo.close();
  
  return nCPU;
}

#  elif defined(__sun__) || defined(__osf__) || defined(_AIX)

unsigned int MultithreadInfo::get_n_cpus()
{
  return sysconf(_SC_NPROCESSORS_ONLN);
}


#  elif defined(__sgi__)

unsigned int MultithreadInfo::get_n_cpus()
{
  return sysconf(_SC_NPROC_ONLN);
}

#  elif defined(__MACH__) && defined(__ppc__) && defined(__APPLE__)
// This is only tested on a dual G5 2.5GHz running MacOSX 10.3.6 
// and gcc version 3.3 20030304 (Apple Computer, Inc. build 1666)
// If it doesnt work please contact the mailinglist.
unsigned int MultithreadInfo::get_n_cpus()
{
	int mib[2];
	int n_cpus;
	size_t len;
	
	mib[0] = CTL_HW;
	mib[1] = HW_NCPU;
	len = sizeof(n_cpus);
	sysctl(mib, 2, &n_cpus, &len, NULL, 0);
	
	return n_cpus;
}


#  else

// If you get n_cpus=1 although you are on a multi-processor machine,
// then this may have two reasons: either because the system macros,
// e.g.__linux__, __sgi__, etc. weren't defined by the compiler or the
// detection of processors is really not implemented for your specific
// system. In the first case you can add e.g. -D__sgi__ to your
// compiling flags, in the latter case you need to implement the
// get_n_cpus() function for your system.
//
// In both cases, this #else case is compiled, a fact that you can
// easily verify by uncommenting the following #error directive,
// recompiling and getting a compilation error right at that line.
// After definition of the system macro or the implementation of the
// new detection this #error message during compilation shouldn't
// occur any more.
//
// Please send all new implementations of detection of processors to
// the deal.II mailing list, such that it can be included into the
// next deal.II release.

//#error Detection of Processors not supported on this OS. Setting n_cpus=1 by default.

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#  endif


#else				 // not in MT mode

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#endif


MultithreadInfo::MultithreadInfo ()
                :
                n_cpus (get_n_cpus()),
                n_default_threads (n_cpus)
{}



unsigned int
MultithreadInfo::memory_consumption ()
{
				   // only simple data elements, so
				   // use sizeof operator
  return sizeof (MultithreadInfo);
}



// definition of the variable which is declared `extern' in the .h file
MultithreadInfo multithread_info;

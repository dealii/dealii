//----------------------------  multithread_info.cc  ----------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
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

#if defined(__sun__)
#  include <unistd.h>
#endif




#ifdef DEAL_II_USE_MT

#if defined(__linux__)

unsigned int MultithreadInfo::get_n_cpus()
{
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


#elif defined(__sun__) || defined(__osf__)


unsigned int MultithreadInfo::get_n_cpus()
{
  return sysconf(_SC_NPROCESSORS_ONLN);
}


#else

// if you get to see this warning, then this can have two reasons: either
// because you fell through all of the above #if clauses and on your system
// the detection of processors is really not implemented, or you are using
// a compiler that does not understand the #warning directive, in which case
// you will see a report on the screen even if the detection of processors
// is implemented for your system. In the latter case, you need not worry
// about the warning (although it is acknowledged that it is annoying),
// otherwise you might want to consider implementing the missing feature
// and submitting it back to the authors of the library.
#warning Detection of Processors not supported on this OS

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#endif


#else				 // not in multithreadmode

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#endif


MultithreadInfo::MultithreadInfo () :
                n_cpus (get_n_cpus()),
                n_default_threads (n_cpus)
{};



unsigned int
MultithreadInfo::memory_consumption ()
{
				   // only simple data elements, so
				   // use sizeof operator
  return sizeof (MultithreadInfo);
};



// definition of the variable which is declared `extern' in the .h file
MultithreadInfo multithread_info;

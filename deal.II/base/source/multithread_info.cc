//----------------------------  multithread_info.cc  ----------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  multithread_info.cc  ----------------

#include <base/multithread_info.h>

#if defined(__linux__)

#include <fstream.h>
#include <string>
#endif

#if defined(__sun__)
#include <unistd.h>
#endif

#ifdef DEAL_II_USE_MT

#if defined(__linux__)

unsigned int MultithreadInfo::get_n_cpus()
{
  ifstream cpuinfo;
  string search;
  unsigned int nCPU = 0;
  
  cpuinfo.open("/proc/cpuinfo");

  AssertThrow(cpuinfo,ExcProcNotPresent());
  
  while((cpuinfo))
    {
      cpuinfo >> search;
      if (search.find("processor")!=string::npos)
	nCPU++;	  
    }
  cpuinfo.close();
  
  return nCPU;
}

#elif defined(__sun__)

unsigned int MultithreadInfo::get_n_cpus()
{
  return sysconf(_SC_NPROCESSORS_ONLN);
}

#else

#warning Detection of Processors not supported on this OS

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#endif

				 // not in multithreadmode
#else

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#endif

MultithreadInfo multithread_info;

MultithreadInfo::MultithreadInfo () :
                n_cpus (get_n_cpus()),
                n_default_threads (n_cpus)
{};






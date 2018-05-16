// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

/*
  Test that OpenMP is not messing with thread affinity, which will stop TBB
  from creating threads.
 */


#include <deal.II/grid/tria.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>

#if defined(__linux__)
#include <sched.h>
#include <sys/sysinfo.h>
#endif

bool
getaffinity(unsigned int &bits_set,unsigned int &mask)
{
  bits_set = 0;
  mask = 0x00;

#if defined(__linux__)
  cpu_set_t my_set;
  CPU_ZERO(&my_set);

  unsigned int len = sizeof(my_set);
  int   ret = sched_getaffinity(0, len, &my_set);

  if (ret!=0)
    {
      printf("sched_getaffinity() failed, return value: %d\n", ret);
      return false;
    }
  for (int i=0; i<CPU_SETSIZE; ++i)
    bits_set += CPU_ISSET(i,&my_set);

  mask = *(int *)(&my_set);
#else
  // sadly we don't have an implementation
  // for mac/windows
#endif
  return true;
}

int
get_num_thread_env()
{
  const char *penv = getenv ("DEAL_II_NUM_THREADS");
  if (penv!=nullptr)
    {
      int max_threads_env = -1;
      try
        {
          max_threads_env = dealii::Utilities::string_to_int(std::string(penv));
        }
      catch (...)
        {
          return -1;
        }
      return max_threads_env;
    }

  return -1;
}


int
main ()
{
  // we need this, otherwise gcc will not link against deal.II
  dealii::Triangulation<2> test;

  unsigned int bits_set, mask;
  if (!getaffinity(bits_set, mask))
    return 1;

  unsigned int nprocs = dealii::MultithreadInfo::n_cores();
  unsigned int tbbprocs = dealii::MultithreadInfo::n_threads();
  int env = get_num_thread_env();
  printf("aff_ncpus=%d, mask=%08X, nprocs=%d, tbb_threads=%d, DEAL_II_NUM_THREADS=%d\n",
         bits_set, mask, nprocs, tbbprocs, env );

  if (bits_set !=0  && bits_set!=nprocs)
    {
      printf("Warning: sched_getaffinity() returns that we can only use %d out of %d CPUs.\n",bits_set, nprocs);
      return 2;
    }
  if (env != -1 && nprocs != tbbprocs)
    {
      printf("Warning: number of threads is set to %d in environment using DEAL_II_NUM_THREADS.\n", env);
      return 0; // do not return an error!
    }
  if (nprocs != tbbprocs)
    {
      printf("Warning: for some reason TBB only wants to use %d out of %d CPUs.\n",
             tbbprocs, nprocs);
      return 3;
    }

  return 0;
}

//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/thread_management.h>
#include <iostream>
#include <cstdlib>

#ifdef DEAL_II_USE_MT_POSIX
#  include <list>
#endif

#ifndef DEAL_II_USE_DIRECT_ERRNO_H
#  include <errno.h>
#else
#  include </usr/include/errno.h>
#endif
#include <sys/errno.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_SYS_SYSCALL_H
#  include <sys/syscall.h>
#endif

DEAL_II_NAMESPACE_OPEN


namespace Threads 
{
  namespace internal
  {
                                     // counter and access mutex for the
                                     // number of threads
    volatile unsigned int n_existing_threads_counter = 1;
    ThreadMutex  n_existing_threads_mutex;

  
    void register_thread ()
    {
      ThreadMutex::ScopedLock lock (n_existing_threads_mutex);
      ++n_existing_threads_counter;
    }


  
    void deregister_thread () 
    {
      ThreadMutex::ScopedLock lock (n_existing_threads_mutex);
      --n_existing_threads_counter;
      Assert (n_existing_threads_counter >= 1,
              ExcInternalError());
    }



    void handle_std_exception (const std::exception &exc) 
    {
      std::cerr << std::endl << std::endl
                << "---------------------------------------------------------"
                << std::endl
                << "In one of the sub-threads of this program, an exception\n"
                << "was thrown and not caught. Since exceptions do not\n"
                << "propagate to the main thread, the library has caught it.\n"
                << "The information carried by this exception is given below.\n"
                << std::endl
                << "---------------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "---------------------------------------------------------"
                << std::endl;
      std::abort ();
    }



    void handle_unknown_exception ()
    {
      std::cerr << std::endl << std::endl
                << "---------------------------------------------------------"
                << std::endl
                << "In one of the sub-threads of this program, an exception\n"
                << "was thrown and not caught. Since exceptions do not\n"
                << "propagate to the main thread, the library has caught it.\n"
                << "The information carried by this exception is given below.\n"
                << std::endl
                << "---------------------------------------------------------"
                << std::endl;
      std::cerr << "Type of exception is unknown, but not std::exception.\n"
                << "No additional information is available.\n"
                << "---------------------------------------------------------"
                << std::endl;
      std::abort ();
    }
  }
  

  
  unsigned int n_existing_threads () 
  {
    ThreadMutex::ScopedLock lock (internal::n_existing_threads_mutex);
    return internal::n_existing_threads_counter;
  }
  

  unsigned int this_thread_id ()
  {
#ifdef SYS_gettid
    const int this_id = syscall(SYS_gettid);
#elif HAVE_GETPID
    const pid_t this_id = getpid();
#else
    const pid_t this_id = 0;
#endif
    return static_cast<unsigned int>(this_id);
  }


  
#if DEAL_II_USE_MT != 1
  DummyBarrier::DummyBarrier (const unsigned int  count,
			      const char         *,
			      void               *)
  {
    Assert (count == 1, ExcBarrierSizeNotUseful(count));
  }


#else
#  ifdef DEAL_II_USE_MT_POSIX


#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS    
  PosixThreadBarrier::PosixThreadBarrier (const unsigned int  count,
					  const char         *,
					  void               *)
  {
    pthread_barrier_init (&barrier, 0, count);
  }

#else

  PosixThreadBarrier::PosixThreadBarrier (const unsigned int  count,
					  const char         *,
					  void               *)
                  : count (count)
  {
                                     // throw an exception unless we
                                     // have the special case that a
                                     // count of 1 is given, since
                                     // then waiting for a barrier is
                                     // a no-op, and we don't need the
                                     // POSIX functionality
    AssertThrow (count == 1,
		 ExcMessage ("Your local POSIX installation does not support\n"
			     "POSIX barriers. You will not be able to use\n"
			     "this class, but the rest of the threading\n"
			     "functionality is available."));
  }
#endif



  PosixThreadBarrier::~PosixThreadBarrier ()
  {
#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS
    pthread_barrier_destroy (&barrier);
#else
                                     // unless the barrier is a no-op,
                                     // complain again (how did we get
                                     // here then?)
    if (count != 1)
      std::abort ();
#endif
  }



  int
  PosixThreadBarrier::wait ()
  {
#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS    
    return pthread_barrier_wait (&barrier);
#else
                                     // in the special case, this
                                     // function is a no-op. otherwise
                                     // complain about the missing
                                     // POSIX functions
    if (count == 1)
      return 0;
    else
      {
        std::abort ();
        return 1;
      };
#endif
  }
  


  
#  endif
#endif  
  


  std::vector<std::pair<unsigned int,unsigned int> >
  split_interval (const unsigned int begin,
		  const unsigned int end,
		  const unsigned int n_intervals)
  {
    Assert (end >= begin, ExcInternalError());
    
    const unsigned int n_elements              = end-begin;
    const unsigned int n_elements_per_interval = n_elements / n_intervals;
    const unsigned int residual                = n_elements % n_intervals;
    
    std::vector<std::pair<unsigned int,unsigned int> > return_values (n_intervals);

    return_values[0].first = begin;
    for (unsigned int i=0; i<n_intervals; ++i)
      {
	if (i != n_intervals-1) 
	  {
	    return_values[i].second = (return_values[i].first
				       + n_elements_per_interval);
					     // distribute residual in
					     // division equally among
					     // the first few
					     // subintervals
	    if (i < residual)
	      ++return_values[i].second;
	    return_values[i+1].first = return_values[i].second;
	  }
	else
	  return_values[i].second = end;
      };
    return return_values;
  }
  
}   // end namespace Thread

DEAL_II_NAMESPACE_CLOSE

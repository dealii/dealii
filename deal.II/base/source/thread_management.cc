//----------------------------  thread_management.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  thread_management.cc  ---------------------------


#include <base/thread_management.h>
#include <iostream>
#ifdef DEAL_II_USE_MT_POSIX
#  include <list>
#endif

#ifndef DEAL_II_USE_DIRECT_ERRNO_H
#  include <errno.h>
#else
#  include </usr/include/errno.h>
#endif
#include <sys/errno.h>



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
    const unsigned int n = internal::n_existing_threads_counter;
    return n;
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

  PosixThreadMutex::PosixThreadMutex ()
  {
    pthread_mutex_init (&mutex, 0);
  }



  PosixThreadMutex::~PosixThreadMutex ()
  {
    pthread_mutex_destroy (&mutex);
  }


  PosixThreadCondition::PosixThreadCondition ()
  {
    pthread_cond_init (&cond, 0);
  }



  PosixThreadCondition::~PosixThreadCondition ()
  {
    pthread_cond_destroy (&cond);
  }
  

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


#if (DEAL_II_USE_MT == 1) && defined(DEAL_II_USE_MT_POSIX)

  namespace internal
  {
    thread_description_base::~thread_description_base ()
    {
                                       // if we are here, then the
                                       // last listener is just about
                                       // to cease to exists. there
                                       // are two possibilities (since
                                       // there are two types of
                                       // possible listeners, and
                                       // exactly one listener is
                                       // around):
                                       //
                                       // 1. there are no Thread<>
                                       // objects around any more, but
                                       // we are at the end of the
                                       // entry_point function of the
                                       // new thread; the thread is
                                       // still running, though, and
                                       // we are called from this new
                                       // thread.
                                       //
                                       // if this is the case,
                                       // apparently none of the
                                       // orginal listeners has called
                                       // join(), otherwise the thread
                                       // would not be running any
                                       // more. since there are no
                                       // listeners any more, nobody
                                       // will ever call join() on
                                       // this thread any
                                       // more. however, to avoid a
                                       // resource leak, somebody has
                                       // to do something like that --
                                       // on the other hand, we can't
                                       // since calling pthread_join()
                                       // on ourself will create a
                                       // deadlock, of course. so
                                       // detach the thread.
                                       //
                                       // note that in this case,
                                       // since nobody has ever called
                                       // join(), was_joined must be
                                       // false
                                       //
                                       // 2. the thread has already
                                       // ended. in this case, the
                                       // last Thread<> object
                                       // refering to it just went out
                                       // of scope. If this is the
                                       // case, to prevent a memory
                                       // leak, we have to call join
                                       // on the thread, if this has
                                       // not yet happened, i.e. if
                                       // was_joined is false.
                                       //
                                       // if we are in case 2, then
                                       // the destructor is called
                                       // from another than the
                                       // presently running thread
                                       //
                                       //
                                       // now: how do we find out
                                       // which case we're in? we
                                       // could either call
                                       // pthread_join and see whether
                                       // we get an error back, or if
                                       // we get back an error on
                                       // pthread_detach. from the man
                                       // pages, the second way seems
                                       // a little bit safer, so go it:
      if (was_joined == false)
        {
                                           // assume case 1:
          int error = pthread_detach (thread);
          if (error == 0)
            return;

                                           // ouch, could not
                                           // detach. see if thread
                                           // could not be found any
                                           // more:
          if (error == ESRCH)
                                             // ok, this is the
                                             // case. then we are in
                                             // branch 2, and need to
                                             // join the thread
            join ();
          else
                                             // something went
                                             // terribly wrong
            AssertThrow (false, ExcInternalError());
        }
    }

    
    void
    thread_description_base::create (void * (*p) (void *), void *d)
    {
                                       // start new thread. retry until
                                       // we either succeed or get an
                                       // error other than EAGAIN
      int error = 0;
      do
        {
          error = pthread_create (&thread, 0, p, d);
        }
      while (error == EAGAIN);

      AssertThrow (error == 0, ExcInternalError());
    }
      
    

    void
    thread_description_base::join () const
    {
                                       // use Schmidt's double
                                       // checking pattern: if thread
                                       // was already joined, then
                                       // return immediately
      if (was_joined)
        return;

                                       // otherwise make sure that
                                       // only one thread can enter
                                       // the following section at a
                                       // time
      ThreadMutex::ScopedLock lock(join_mutex);
                                       // while getting the lock,
                                       // another thread may have set
                                       // was_joined to true, so check
                                       // again (this is the double
                                       // checking pattern) and only
                                       // if we are really sure that
                                       // this has not happened, call
                                       // pthread_join
      if (!was_joined)
        {
          const int error = pthread_join (thread, 0);
          AssertThrow (error == 0, ExcInternalError());
        }

                                       // set the flag
      was_joined = true;
    }
    
  } // end namespace internal

#endif  // (DEAL_II_USE_MT == 1) && defined(DEAL_II_USE_MT_POSIX)
  
}   // end namespace Thread

//----------------------------  thread_management.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  thread_management.cc  ---------------------------


#include <base/thread_management.h>
#ifdef DEAL_II_USE_MT_POSIX
#  include <list>
#endif


namespace Threads 
{
#ifndef DEAL_II_USE_MT
  void DummyThreadManager::spawn (const FunPtr fun_ptr,
				  void *       fun_data,
				  int          /*flags*/) const
  {
    (*fun_ptr) (fun_data);
  };
  

  
  DummyBarrier::DummyBarrier (const unsigned int  count,
			      const char         *,
			      void               *)
  {
    Assert (count == 1, ExcBarrierSizeNotUseful(count));
  };


#else
#  ifdef DEAL_II_USE_MT_POSIX

  PosixThreadMutex::PosixThreadMutex ()
  {
    pthread_mutex_init (&mutex, 0);
  };



  PosixThreadMutex::~PosixThreadMutex ()
  {
    pthread_mutex_destroy (&mutex);
  };
  

  PosixThreadBarrier::PosixThreadBarrier (const unsigned int  count,
					  const char         *,
					  void               *)
  {
#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS    
    pthread_barrier_init (&barrier, 0, count);
#else
    AssertThrow (false,
		 ExcMessage ("Your local POSIX installation does not support\n"
			     "POSIX barriers. You will not be able to use\n"
			     "this class, but the rest of the threading\n"
			     "functionality is available."));
#endif
  };



  PosixThreadBarrier::~PosixThreadBarrier ()
  {
#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS
    pthread_barrier_destroy (&barrier);
#else
    std::abort ();
#endif
  };



  int
  PosixThreadBarrier::wait ()
  {
#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS    
    return pthread_barrier_wait (&barrier);
#else
    std::abort ();
    return 1;
#endif
  };
  


  PosixThreadManager::PosixThreadManager ()
		  :
		  thread_id_list (new std::list<pthread_t>())
  {};


  PosixThreadManager::~PosixThreadManager ()
  {
				     // wait for all threads, and
				     // release memory
    wait ();
    delete reinterpret_cast<std::list<pthread_t>*>(thread_id_list);
  };



  void
  PosixThreadManager::spawn (const FunPtr fun_ptr,
			     void *       fun_data,
			     int)
  {
    std::list<pthread_t> &tid_list
      = *reinterpret_cast<std::list<pthread_t>*>(thread_id_list);

    tid_list.push_back (pthread_t());

    const int error
      = pthread_create (&tid_list.back(), 0, fun_ptr, fun_data);

    Assert (error == 0, ExcInternalError());
  };
  


  void
  PosixThreadManager::wait () const
  {
    std::list<pthread_t> &tid_list
      = *reinterpret_cast<std::list<pthread_t>*>(thread_id_list);

    for (std::list<pthread_t>::iterator i=tid_list.begin();
	 i != tid_list.end(); ++i)
      pthread_join (*i, 0);
  };
  
#  endif
#endif  
  
  FunDataCounter::FunDataCounter () :
		  n_fun_encapsulation_objects (0),
		  n_fun_data_base_objects (0)
  {};
  

  
  FunDataCounter::~FunDataCounter () 
  {   
    AssertThrow (n_fun_encapsulation_objects == 0,
		 ExcObjectsExist("FunEncapsulation", n_fun_encapsulation_objects));
    AssertThrow (n_fun_data_base_objects == 0,
		 ExcObjectsExist("FunDataBase", n_fun_data_base_objects));
  };
      


/**
 * This is the global object which we will use to count the number of
 * threads generated, and which is used to complain when there is a
 * memory leak.
 */
  FunDataCounter fun_data_counter;


  FunEncapsulation::FunEncapsulation () :
		  fun_data_base (0)
  {
				     // keep some statistics on the
				     // number of variables around
    ++fun_data_counter.n_fun_encapsulation_objects;
  };



  FunEncapsulation::FunEncapsulation (FunDataBase *fun_data_base) :
		  fun_data_base (fun_data_base)
  {
				     // keep some statistics on the
				     // number of variables around
    ++fun_data_counter.n_fun_encapsulation_objects;
  };



  FunEncapsulation::FunEncapsulation (const FunEncapsulation &fun_data) :
		  fun_data_base (fun_data.fun_data_base->clone ())
  {
				     // keep some statistics on the
				     // number of variables around
    ++fun_data_counter.n_fun_encapsulation_objects;
  };


  FunEncapsulation::~FunEncapsulation ()
  {
				     // wait until the data in
				     // fun_data_base has been copied
				     // and is no more needed
    fun_data_base->lock.acquire();
				     // now that we have got the lock,
				     // we can delete the data
    delete fun_data_base;
    fun_data_base = 0;

				   // keep some statistics on the
				   // number of variables around
    --fun_data_counter.n_fun_encapsulation_objects;
  };
    

  const FunEncapsulation &
  FunEncapsulation::operator = (const FunEncapsulation &/*fun_data*/)
  {
				     // this is not implemented at
				     // present. return dummy value
				     // instead
    Assert (false, ExcNotImplemented());
    const FunEncapsulation * const p = 0;
    return *p;
  };



  FunDataBase::FunDataBase (const ThreadEntryPoint thread_entry_point) :
		  thread_entry_point (thread_entry_point)
  {
				     // keep some statistics on the
				     // number of variables around
    ++fun_data_counter.n_fun_data_base_objects;
  };



  FunDataBase::FunDataBase (const FunDataBase &fun_data_base) :
		  thread_entry_point (fun_data_base.thread_entry_point)
  {
				     // keep some statistics on the
				     // number of variables around
    ++fun_data_counter.n_fun_data_base_objects;
  };



  FunDataBase::~FunDataBase ()
  {
				     // invalidate pointer for security
				     // reasons. accesses to this
				     // pointer after lifetime of this
				     // object will then fail
    thread_entry_point = 0;

				     // keep some statistics on the
				     // number of variables around
    --fun_data_counter.n_fun_data_base_objects;
  };


    
  void spawn (ThreadManager          &thread_manager,
	      const FunEncapsulation &fun_data)
  {
				     // lock the @p{fun_data_base} object
				     // to avoid destruction while its
				     // data is still accessed. the lock
				     // is released by the new thread
				     // once it has copied all data
    fun_data.fun_data_base->lock.acquire ();
				     // now start the new thread
#ifdef DEAL_II_USE_MT
#  if defined(DEAL_II_USE_MT_ACE)
    thread_manager.spawn (*fun_data.fun_data_base->thread_entry_point,
			  (void*)&fun_data,
			  THR_SCOPE_SYSTEM | THR_DETACHED);
#  elif defined(DEAL_II_USE_MT_POSIX)
    thread_manager.spawn (*fun_data.fun_data_base->thread_entry_point,
			  (void*)&fun_data,
			  0);    
#  endif    
#else
    thread_manager.spawn (*fun_data.fun_data_base->thread_entry_point,
			  (void*)&fun_data,
			  0);
#endif
  };


  
  void spawn_n (ThreadManager          &thread_manager,
		const FunEncapsulation &fun_data,
		const unsigned int      n_threads)
  {
    for (unsigned int i=0; i<n_threads; ++i)
      spawn (thread_manager, fun_data);
  };
  



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
  };  
  
  
 
};   // end namespace Threads

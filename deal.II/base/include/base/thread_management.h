//----------------------------  thread_management.h  ---------------------------
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
//----------------------------  thread_management.h  ---------------------------
#ifndef __deal2__thread_management_h
#define __deal2__thread_management_h


#include <base/config.h>
#include <base/exceptions.h>

#include <iterator>
#include <vector>
#include <utility>

#include <list>
#include <boost_local/type_traits.hpp>
#include <boost_local/tuple/tuple.hpp>
#include <boost_local/shared_ptr.hpp>

#if DEAL_II_USE_MT == 1
#  if defined(DEAL_II_USE_MT_POSIX)
#    include <pthread.h>
#  endif
#endif




namespace Threads 
{
/**
 * This class is used instead of a true lock class when not using
 * multithreading. It allows to write programs such that they start
 * new threads and/or lock objects in multithreading mode, and use
 * dummy thread management and synchronisation classes instead when
 * running in single-thread mode. Specifically, the @p{spawn} functions
 * only call the function but wait for it to return instead of running
 * in on another thread, and the mutices do nothing really. The only
 * reason to provide such a function is that the program can be
 * compiled both in MT and non-MT mode without difference.
 *
 * @author Wolfgang Bangerth, 2000, 2003
 */
  class DummyThreadMutex
  {
    public:
                                       /**
                                        * Scoped lock class. When you
                                        * declare an object of this
                                        * type, you have to pass it a
                                        * mutex, which is locked in
                                        * the constructor of this
                                        * class and unlocked in the
                                        * destructor. The lock is thus
                                        * held during the entire
                                        * lifetime of this object,
                                        * i.e. until the end of the
                                        * present scope, which
                                        * explains where the name
                                        * comes from. This pattern of
                                        * using locks with mutexes
                                        * follows the
                                        * resource-acquisition-is-initialization
                                        * pattern, and was used first
                                        * for mutexes by Doug
                                        * Schmidt. It has the
                                        * advantage that locking a
                                        * mutex this way is
                                        * thread-safe, i.e. when an
                                        * exception is thrown between
                                        * the locking and unlocking
                                        * point, the destructor makes
                                        * sure that the mutex is
                                        * unlocked; this would not
                                        * automatically be the case
                                        * when you lock and unlock the
                                        * mutex "by hand", i.e. using
                                        * @p{acquire} and @p{release}.
                                        */
      class ScopedLock
      {
        public:
                                           /**
                                            * Constructor. Lock the
                                            * mutex. Since this is a
                                            * dummy mutex class, this
                                            * of course does nothing.
                                            */
          ScopedLock (DummyThreadMutex &) {};

                                           /**
                                            * Destructor. Unlock the
                                            * mutex. Since this is a
                                            * dummy mutex class, this
                                            * of course does nothing.
                                            */
          ~ScopedLock () {};
      };
      
				       /**
					* Simulate acquisition of the
					* mutex. As this class does
					* nothing really, this
					* function does nothing as
					* well.
					*/
      inline void acquire () const {};

				       /**
					* Simulate release of the
					* mutex. As this class does
					* nothing really, this
					* function does nothing as
					* well.
					*/
      inline void release () const {};
  };



/**
 * This class is used in single threaded mode instead of a class
 * implementing real condition variable semantics. It allows to write
 * programs such that they start new threads and/or lock objects in
 * multithreading mode, and use dummy thread management and
 * synchronisation classes instead when running in single-thread
 * mode. Specifically, the @p{spawn} functions only call the function
 * but wait for it to return instead of running in on another thread,
 * and the mutices do nothing really. The only reason to provide such
 * a function is that the program can be compiled both in MT and
 * non-MT mode without difference.
 *
 * In this particular case, just as with mutexes, the functions do
 * nothing, and by this provide the same semantics of condition
 * variables as in multi-threaded mode.
 *
 * @author Wolfgang Bangerth, 2003
 */
  class DummyThreadCondition
  {
    public:
				       /**
					* Signal to a single listener
					* that a condition has been
					* met, i.e. that some data
					* will now be available. Since
					* in single threaded mode,
					* this function of course does
					* nothing.
					*/
      inline void signal () {};

				       /**
					* Signal to multiple listener
					* that a condition has been
					* met, i.e. that some data
					* will now be available. Since
					* in single threaded mode,
					* this function of course does
					* nothing.
					*/
      inline void broadcast () {};

				       /**
					* Wait for the condition to be
					* signalled. Signal variables
					* need to be guarded by a
					* mutex which needs to be
					* given to this function as an
					* argument, see the man page
					* of @p{posix_cond_wait} for a
					* description of the
					* mechanisms. Since in single
					* threaded mode, this function
					* of course does nothing, but
					* returns immediately.
					*/
      inline void wait (DummyThreadMutex &) {};
  };


  
/**
 * This class is used instead of a true barrier class when not using
 * multithreading. It allows to write programs such that they use the
 * same class names in multithreading and non-MT mode and thus may be
 * compiled with or without thread-support without the need to use
 * conditional compilation. Since a barrier class only makes sense in
 * non-multithread mode if only one thread is to be synchronised
 * (otherwise, the barrier could not be left, since the one thread is
 * waiting for some other part of the program to reach a certain point
 * of execution), the constructor of this class throws an exception if
 * the @p{count} argument denoting the number of threads that need to
 * be synchronised is not equal to one.
 *
 * @author Wolfgang Bangerth, 2001
 */
  class DummyBarrier
  {
    public:
				       /**
					* Constructor. Since barriers
					* are only useful in
					* single-threaded mode if the
					* number of threads to be
					* synchronised is one, this
					* constructor raises an
					* exception if the @p{count}
					* argument is one.
					*/
      DummyBarrier (const unsigned int  count,
		    const char         *name = 0,
		    void               *arg  = 0);

				       /**
					* Wait for all threads to
					* reach this point. Since
					* there may only be one
					* thread, return immediately,
					* i.e. this function is a
					* no-op.
					*/
      int wait () { return 0; };

				       /**
					* Dump the state of this
					* object. Here: do nothing.
					*/
      void dump () {};

				       /**
					* Exception.
					*/
      DeclException1 (ExcBarrierSizeNotUseful,
		      int,
		      << "In single-thread mode, other barrier sizes than 1 are not "
		      << "useful. You gave " << arg1);
  };
  
  
#if DEAL_II_USE_MT == 1
#  if defined(DEAL_II_USE_MT_POSIX)

				   /**
				    * Class implementing a Mutex with
				    * the help of POSIX functions.
				    *
				    * @author Wolfgang Bangerth, 2002, 2003
				    */
  class PosixThreadMutex 
  {
    public:
                                       /**
                                        * Scoped lock class. When you
                                        * declare an object of this
                                        * type, you have to pass it a
                                        * mutex, which is locked in
                                        * the constructor of this
                                        * class and unlocked in the
                                        * destructor. The lock is thus
                                        * held during the entire
                                        * lifetime of this object,
                                        * i.e. until the end of the
                                        * present scope, which
                                        * explains where the name
                                        * comes from. This pattern of
                                        * using locks with mutexes
                                        * follows the
                                        * resource-acquisition-is-initialization
                                        * pattern, and was used first
                                        * for mutexes by Doug
                                        * Schmidt. It has the
                                        * advantage that locking a
                                        * mutex this way is
                                        * thread-safe, i.e. when an
                                        * exception is thrown between
                                        * the locking and unlocking
                                        * point, the destructor makes
                                        * sure that the mutex is
                                        * unlocked; this would not
                                        * automatically be the case
                                        * when you lock and unlock the
                                        * mutex "by hand", i.e. using
                                        * @p{acquire} and @p{release}.
                                        */
      class ScopedLock
      {
        public:
                                           /**
                                            * Constructor. Lock the
                                            * mutex.
                                            */
          ScopedLock (PosixThreadMutex &m) : mutex(m) { mutex.acquire(); };

                                           /**
                                            * Destructor. Unlock the
                                            * mutex. Since this is a
                                            * dummy mutex class, this
                                            * of course does nothing.
                                            */
          ~ScopedLock () { mutex.release (); };
          
        private:
                                           /**
                                            * Store the address of the
                                            * mutex object.
                                            */
          PosixThreadMutex &mutex;
      };
      
				       /**
					* Constructor. Initialize the
					* underlying POSIX mutex data
					* structure.
					*/
      PosixThreadMutex ();

				       /**
					* Destructor. Release all
					* resources.
					*/
      ~PosixThreadMutex ();
      
				       /**
					* Acquire a mutex.
					*/
      inline void acquire () { pthread_mutex_lock(&mutex); };

				       /**
					* Release the mutex again.
					*/
      inline void release () { pthread_mutex_unlock(&mutex); };

    private:
				       /**
					* Data object storing the
					* POSIX data which we need to
					* call the POSIX functions.
					*/
      pthread_mutex_t mutex;

                                       /**
                                        * Make the class implementing
                                        * condition variables a
                                        * friend, since it needs
                                        * access to the
                                        * @p{pthread_mutex_t}
                                        * structure.
                                        */
      friend class PosixThreadCondition;
  };


				   /**
				    * Class implementing a condition
				    * variable with the help of POSIX
				    * functions. The semantics of this
				    * class and its member functions
				    * are the same as those of the
				    * POSIX functions.
				    *
				    * @author Wolfgang Bangerth, 2003
				    */
  class PosixThreadCondition
  {
    public:
				       /**
					* Constructor. Initialize the
					* underlying POSIX condition
					* variable data structure.
					*/
      PosixThreadCondition ();

				       /**
					* Destructor. Release all
					* resources.
					*/
      ~PosixThreadCondition ();
      
				       /**
					* Signal to a single listener
					* that a condition has been
					* met, i.e. that some data
					* will now be available.
					*/
      inline void signal () { pthread_cond_signal(&cond); };

				       /**
					* Signal to multiple listener
					* that a condition has been
					* met, i.e. that some data
					* will now be available.
					*/
      inline void broadcast () { pthread_cond_broadcast(&cond); };

				       /**
					* Wait for the condition to be
					* signalled. Signal variables
					* need to be guarded by a
					* mutex which needs to be
					* given to this function as an
					* argument, see the man page
					* of @p{posix_cond_wait} for a
					* description of the
					* mechanisms.
					*/
      inline void wait (PosixThreadMutex &mutex)
        { pthread_cond_wait(&cond, &mutex.mutex); };

    private:
				       /**
					* Data object storing the
					* POSIX data which we need to
					* call the POSIX functions.
					*/
      pthread_cond_t cond;
  };


				   /**
				    * Implementation of a thread
				    * barrier class, based on the
				    * POSIX thread functions. POSIX
				    * barriers are a relatively new
				    * feature and are not supported on
				    * all systems.
				    *
				    * If the configuration detected
				    * the absence of these functions,
				    * then barriers will not be
				    * available, and creating objects
				    * of this class will result in an
				    * exception been thrown unless the
				    * count given for the parties
				    * waiting for the barrier is equal
				    * to one (as in this case waiting
				    * for the barrier is a
				    * no-operation, and we can
				    * dispense with the POSIX
				    * functions at all). The rest of
				    * the threading functionality will
				    * be available in its full extent,
				    * though, even if POSIX barriers
				    * are not available.
				    *
				    * @author Wolfgang Bangerth, 2002
				    */
  class PosixThreadBarrier 
  {
    public:
				       /**
					* Constructor. Initialize the
					* underlying POSIX barrier data
					* structure.
					*/
      PosixThreadBarrier (const unsigned int  count,
			  const char         *name = 0,
			  void               *arg  = 0);

				       /**
					* Destructor. Release all
					* resources.
					*/
      ~PosixThreadBarrier ();      

				       /**
					* Wait for all threads to
					* reach this point. The return
					* value is zero for all
					* participating threads except
					* for one, for which the
					* return value is some
					* non-zero value. The
					* operating system picks the
					* special thread by some not
					* further known method.
					*/
      int wait ();

    private:
				       /**
					* Data object storing the
					* POSIX data which we need to
					* call the POSIX functions.
					*/
#ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS
      pthread_barrier_t barrier;
#else
      unsigned int count;
#endif
  };


                                   /**
                                    * If using POSIX functions, then
                                    * alias the POSIX wrapper classes
                                    * to the names we use throughout
                                    * the library.
                                    */
  typedef PosixThreadMutex     ThreadMutex;
  typedef PosixThreadCondition ThreadCondition;  
  typedef PosixThreadBarrier   Barrier;

#  else
#    error Not Implemented
#  endif
#else
				   /**
				    * In non-multithread mode, the
				    * mutex and thread management
				    * classes are aliased to dummy
				    * classes that actually do
				    * nothing, in particular not lock
				    * objects. Likewise for the
				    * barrier class.
				    */
  typedef DummyThreadMutex     ThreadMutex;
  typedef DummyThreadCondition ThreadCondition;  
  typedef DummyBarrier         Barrier;
#endif
  
}


namespace Threads
{

                                   /**
                                    * Return the number of presently
                                    * existing threads. This function
                                    * may be useful in a situation
                                    * where a large number of threads
                                    * are concurrently, and you want
                                    * to decide whether creation of
                                    * another thread is reasonable or
                                    * whether running the respective
                                    * operation sequentially is more
                                    * useful since already many more
                                    * threads than processors are
                                    * running.
                                    *
                                    * Note that the function returns
                                    * the total number of threads, not
                                    * those actually running. Some of
                                    * the threads may be waiting for
                                    * locks and mutexes, or may be
                                    * sleeping until they are
                                    * delivered with data to work on.
                                    * 
                                    * Upon program start, this number
                                    * is one. It is increased each
                                    * time a thread is created using
                                    * the @ref{Threads::spawn} or
                                    * @ref{Threads::spawn_n}
                                    * functions. It is decreased once
                                    * a thread terminates by returning
                                    * from the function that was
                                    * spawned.
                                    *
                                    * Note that this means that only
                                    * threads created and terminated
                                    * through the interfaces of this
                                    * namespace are taken care of. If
                                    * threads are created by directly
                                    * calling the respective functions
                                    * of the operating system
                                    * (e.g. @p{pthread_create} for the
                                    * POSIX thread interface), or if
                                    * they are killed (e.g. either
                                    * through @p{pthread_exit} from
                                    * the spawned thread, or
                                    * @p{pthread_kill} from another
                                    * thread), then these events are
                                    * not registered and counted for
                                    * the result of this function.
                                    */
  unsigned int n_existing_threads ();

				   /**
				    * Split the range @p{[begin,end)}
				    * into @p{n_intervals} subintervals
				    * of equal size. The last interval
				    * will be a little bit larger, if
				    * the number of elements in the
				    * whole range is not exactly
				    * divisible by @p{n_intervals}. The
				    * type of the iterators has to
				    * fulfill the requirements of a
				    * forward iterator,
				    * i.e. @p{operator++} must be
				    * available, and of course it must
				    * be assignable.
				    *
				    * A list of subintervals is
				    * returned as a vector of pairs of
				    * iterators, where each pair
				    * denotes the range
				    * @p{[begin[i],end[i])}.
				    */
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator,ForwardIterator> >
  split_range (const ForwardIterator &begin,
	       const ForwardIterator &end,
	       const unsigned int n_intervals);

				   /**
				    * Split the interval @p{[begin,end)}
				    * into subintervals of (almost)
				    * equal size. This function works
				    * mostly as the one before, with
				    * the difference that instead of
				    * iterators, now values are taken
				    * that define the whole interval.
				    */
  std::vector<std::pair<unsigned int,unsigned int> >
  split_interval (const unsigned int begin,
		  const unsigned int end,
		  const unsigned int n_intervals);

                                   /**
                                    * A namespace in which helper
                                    * functions and the like for the
                                    * threading subsystem are
                                    * implemented. The members of this
                                    * namespace are not meant for
                                    * public use.
                                    *
                                    * The classes inside this
                                    * namespace are suppressed by the
                                    * documentation script in the
                                    * class overview table, to keep it
                                    * short.
                                    * 
                                    * @author Wolfgang Bangerth, 2003
                                    */
  namespace internal
  {
                                     /**
                                      * If in a sub-thread an
                                      * exception is thrown, it is not
                                      * propagated to the main
                                      * thread. Therefore, the
                                      * exception handler that is
                                      * provided by the applications
                                      * main function or some of its
                                      * other parts will not be able
                                      * to catch these
                                      * exceptions. Therefore, we have
                                      * to provide an exception
                                      * handler in the top function of
                                      * each sub-thread that at least
                                      * catches the exception and
                                      * prints some information,
                                      * rather than letting the
                                      * operating system to just kill
                                      * the program without a
                                      * message. In each of the
                                      * functions we use as entry
                                      * points to new threads, we
                                      * therefore install a try-catch
                                      * block, and if an exception of
                                      * type @p{std::exception} is
                                      * caught, it passes over control
                                      * to this function, which will
                                      * then provide some output.
                                      */
    void handle_std_exception (const std::exception &exc);

                                     /**
                                      * Same as above, but the type of
                                      * the exception is not derived
                                      * from @p{std::exception}, so
                                      * there is little way to provide
                                      * something more useful.
                                      */
    void handle_unknown_exception ();

                                     /**
                                      * The following function is used
                                      * for internal bookkeeping of the
                                      * number of existing threads. It
                                      * is not thought for use in
                                      * application programs, but only
                                      * for use in the template
                                      * functions below.
                                      */
    void register_thread ();
  
                                     /**
                                      * The following function is used
                                      * for internal bookkeeping of the
                                      * number of existing threads. It
                                      * is not thought for use in
                                      * application programs, but only
                                      * for use in the template
                                      * functions below.
                                      */
    void deregister_thread ();  
  }
  
}   // end declarations of namespace Threads


/* ----------- implementation of functions in namespace Threads ---------- */
namespace Threads 
{  
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator,ForwardIterator> >
  split_range (const ForwardIterator &begin,
	       const ForwardIterator &end,
	       const unsigned int     n_intervals)
  {
    typedef std::pair<ForwardIterator,ForwardIterator> IteratorPair;

				     // in non-multithreaded mode, we
				     // often have the case that this
				     // function is called with
				     // n_intervals==1, so have a
				     // shortcut here to handle that
				     // case efficiently
    
    if (n_intervals==1)
      return (std::vector<IteratorPair>
	      (1, IteratorPair(begin, end)));
    
				     // if more than one interval
				     // requested, do the full work
    const unsigned int n_elements              = distance (begin, end);
    const unsigned int n_elements_per_interval = n_elements / n_intervals;
    const unsigned int residual                = n_elements % n_intervals;
    
    std::vector<IteratorPair> return_values (n_intervals);

    return_values[0].first = begin;
    for (unsigned int i=0; i<n_intervals; ++i)
      {
	if (i != n_intervals-1) 
	  {
	    return_values[i].second = return_values[i].first;
					     // note: the cast is
					     // performed to avoid a
					     // warning of gcc that in
					     // the library `dist>=0'
					     // is checked (dist has a
					     // template type, which
					     // here is unsigned if no
					     // cast is performed)
	    std::advance (return_values[i].second,
			  static_cast<signed int>(n_elements_per_interval));
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
      }
    return return_values;
  }  


  namespace internal 
  {
                                     /**
                                      * A type that is used to
                                      * distinguish argument lists of
                                      * functions by enumeration.
                                      */
    template <int> struct int2type
    {
    };
  } 


  namespace internal 
  {
                                     /**
                                      * Given an arbitrary type RT,
                                      * store an element of it and grant
                                      * access to it through functions
                                      * get() and set(). There are
                                      * specializations for reference
                                      * types (which cannot be set), and
                                      * for type void.
                                      */
    template <typename RT> struct return_value
    {
      private:
        RT value;
      public:
        inline return_value () : value() {}
      
        inline RT get () const { return value; }
        inline void set (RT v) { value = v; }
    };

  
                                     /**
                                      * Given an arbitrary type RT,
                                      * store an element of it and grant
                                      * access to it through functions
                                      * get() and set(). This is the
                                      * specialization for reference
                                      * types: since they cannot be set
                                      * after construction time, we
                                      * store a pointer instead, that
                                      * holds the address of the object
                                      * being referenced.
                                      */
    template <typename RT> struct return_value<RT &>
    {
      private:
        RT * value;
      public:
        inline return_value () : value(0) {}

        inline RT & get () const { return *value; }
        inline void set (RT & v) { value = &v; }
    };

  
                                     /**
                                      * Given an arbitrary type RT,
                                      * store an element of it and grant
                                      * access to it through functions
                                      * get() and set(). This is the
                                      * specialization for type void:
                                      * there is obviously nothing to
                                      * store, so no function set(), and
                                      * a function get() that returns
                                      * void.
                                      */
    template <> struct return_value<void> {
        static inline void get () {}
    };
  }

  

  namespace internal
  {
                                     /**
                                      * Call arbitrary functions with
                                      * return type RT. For each number
                                      * of arguments to these functions,
                                      * there is an instance of the
                                      * do_call function in this class
                                      * that unpacks the argument list
                                      * (which is passed by reference)
                                      * and calls the function. The
                                      * number of arguments is
                                      * distinguished by the last
                                      * argument. The return value
                                      * object is the second last. A
                                      * second version of the do_call
                                      * function is used to call member
                                      * function pointers, in which case
                                      * there is an additional argument
                                      * at second position holding a
                                      * reference to the object with
                                      * which the member function
                                      * pointer is to be called.
                                      *
                                      * There is a specialization of
                                      * this class for the case that the
                                      * return type is void, in which
                                      * case there is no return value to
                                      * be set.
                                      */
    template <typename RT>
    struct Caller 
    {
                                         /**
                                          * Call a function with 0
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<0> &)
          {
            ret_val.set ((*fun_ptr) ());
          }

                                         /**
                                          * Call a member function with
                                          * 0 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<0> &)
          {
            ret_val.set ((obj.*fun_ptr) ());
          }


                                         /**
                                          * Call a function with 1
                                          * argument, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<1> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>()));
          }

                                         /**
                                          * Call a member function with
                                          * 1 argument, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<1> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>()));
          }




                                         /**
                                          * Call a function with 2
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<2> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>()));
          }

                                         /**
                                          * Call a member function with
                                          * 2 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<2> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>()));
          }


                                         /**
                                          * Call a function with 3
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<3> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>()));
          }

                                         /**
                                          * Call a member function with
                                          * 3 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<3> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>()));
          }


                                         /**
                                          * Call a function with 4
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<4> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>()));
          }

                                         /**
                                          * Call a member function with
                                          * 4 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<4> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>()));
          }


                                         /**
                                          * Call a function with 5
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<5> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>(),
                                     arg_list.template get<4>()));
          }

                                         /**
                                          * Call a member function with
                                          * 5 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<5> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>(),
                                         arg_list.template get<4>()));
          }


                                         /**
                                          * Call a function with 6
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<6> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>(),
                                     arg_list.template get<4>(),
                                     arg_list.template get<5>()));
          }

                                         /**
                                          * Call a member function with
                                          * 6 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<6> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>(),
                                         arg_list.template get<4>(),
                                         arg_list.template get<5>()));
          }


                                         /**
                                          * Call a function with 7
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<7> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>(),
                                     arg_list.template get<4>(),
                                     arg_list.template get<5>(),
                                     arg_list.template get<6>()));
          }

                                         /**
                                          * Call a member function with
                                          * 7 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<7> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>(),
                                         arg_list.template get<4>(),
                                         arg_list.template get<5>(),
                                         arg_list.template get<6>()));
          }


                                         /**
                                          * Call a function with 8
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<8> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>(),
                                     arg_list.template get<4>(),
                                     arg_list.template get<5>(),
                                     arg_list.template get<6>(),
                                     arg_list.template get<7>()));
          }

                                         /**
                                          * Call a member function with
                                          * 8 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<8> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>(),
                                         arg_list.template get<4>(),
                                         arg_list.template get<5>(),
                                         arg_list.template get<6>(),
                                         arg_list.template get<7>()));
          }


                                         /**
                                          * Call a function with 9
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<9> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>(),
                                     arg_list.template get<4>(),
                                     arg_list.template get<5>(),
                                     arg_list.template get<6>(),
                                     arg_list.template get<7>(),
                                     arg_list.template get<8>()));
          }

                                         /**
                                          * Call a member function with
                                          * 9 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<9> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>(),
                                         arg_list.template get<4>(),
                                         arg_list.template get<5>(),
                                         arg_list.template get<6>(),
                                         arg_list.template get<7>(),
                                         arg_list.template get<8>()));
          }


                                         /**
                                          * Call a function with 10
                                          * arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<10> &)
          {
            ret_val.set ((*fun_ptr) (arg_list.template get<0>(),
                                     arg_list.template get<1>(),
                                     arg_list.template get<2>(),
                                     arg_list.template get<3>(),
                                     arg_list.template get<4>(),
                                     arg_list.template get<5>(),
                                     arg_list.template get<6>(),
                                     arg_list.template get<7>(),
                                     arg_list.template get<8>(),
                                     arg_list.template get<9>()));
          }

                                         /**
                                          * Call a member function with
                                          * 10 arguments, and set the
                                          * return value.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<RT> &ret_val,
                                    const int2type<10> &)
          {
            ret_val.set ((obj.*fun_ptr) (arg_list.template get<0>(),
                                         arg_list.template get<1>(),
                                         arg_list.template get<2>(),
                                         arg_list.template get<3>(),
                                         arg_list.template get<4>(),
                                         arg_list.template get<5>(),
                                         arg_list.template get<6>(),
                                         arg_list.template get<7>(),
                                         arg_list.template get<8>(),
                                         arg_list.template get<9>()));
          }
    };



  
                                     /**
                                      * Call arbitrary functions with
                                      * void return type. For each
                                      * number of arguments to these
                                      * functions, there is an instance
                                      * of the do_call function in this
                                      * class that unpacks the argument
                                      * list (which is passed by
                                      * reference) and calls the
                                      * function. The number of
                                      * arguments is distinguished by
                                      * the last argument. The return
                                      * value object is the second last,
                                      * but since the return type is
                                      * void, this is of course simply
                                      * ignored. A second version of the
                                      * do_call function is used to call
                                      * member function pointers, in
                                      * which case there is an
                                      * additional argument at second
                                      * position holding a reference to
                                      * the object with which the member
                                      * function pointer is to be
                                      * called.
                                      */
    template <>
    struct Caller<void>
    {
                                         /**
                                          * Call a void function with 0
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &,
                                    internal::return_value<void> &,
                                    const int2type<0> &)
          {
            (*fun_ptr) ();
          }

                                         /**
                                          * Call a void member function
                                          * with 0 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &,
                                    internal::return_value<void> &,
                                    const int2type<0> &)
          {
            (obj.*fun_ptr) ();
          }


                                         /**
                                          * Call a void function with 1
                                          * argument.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<1> &)
          {
            (*fun_ptr) (arg_list.template get<0>());
          }

                                         /**
                                          * Call a void member function
                                          * with 1 argument.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<1> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>());
          }




                                         /**
                                          * Call a void function with 2
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<2> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>());
          }

                                         /**
                                          * Call a void member function
                                          * with 2 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<2> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>());
          }


                                         /**
                                          * Call a void function with 3
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<3> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>());
          }

                                         /**
                                          * Call a void member function
                                          * with 3 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<3> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>());
          }


                                         /**
                                          * Call a void function with 4
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<4> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>());
          }

                                         /**
                                          * Call a void member function
                                          * with 4 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<4> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>());
          }


                                         /**
                                          * Call a void function with 5
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<5> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>(),
                        arg_list.template get<4>());
          }

                                         /**
                                          * Call a void member function
                                          * with 5 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<5> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>(),
                            arg_list.template get<4>());
          }


                                         /**
                                          * Call a void function with 6
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<6> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>(),
                        arg_list.template get<4>(),
                        arg_list.template get<5>());
          }

                                         /**
                                          * Call a void member function
                                          * with 6 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<6> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>(),
                            arg_list.template get<4>(),
                            arg_list.template get<5>());
          }


                                         /**
                                          * Call a void function with 7
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<7> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>(),
                        arg_list.template get<4>(),
                        arg_list.template get<5>(),
                        arg_list.template get<6>());
          }

                                         /**
                                          * Call a void member function
                                          * with 7 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<7> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>(),
                            arg_list.template get<4>(),
                            arg_list.template get<5>(),
                            arg_list.template get<6>());
          }


                                         /**
                                          * Call a void function with 8
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<8> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>(),
                        arg_list.template get<4>(),
                        arg_list.template get<5>(),
                        arg_list.template get<6>(),
                        arg_list.template get<7>());
          }

                                         /**
                                          * Call a void member function
                                          * with 8 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<8> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>(),
                            arg_list.template get<4>(),
                            arg_list.template get<5>(),
                            arg_list.template get<6>(),
                            arg_list.template get<7>());
          }


                                         /**
                                          * Call a void function with 9
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<9> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>(),
                        arg_list.template get<4>(),
                        arg_list.template get<5>(),
                        arg_list.template get<6>(),
                        arg_list.template get<7>(),
                        arg_list.template get<8>());
          }

                                         /**
                                          * Call a void member function
                                          * with 9 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<9> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>(),
                            arg_list.template get<4>(),
                            arg_list.template get<5>(),
                            arg_list.template get<6>(),
                            arg_list.template get<7>(),
                            arg_list.template get<8>());
          }


                                         /**
                                          * Call a void function with 10
                                          * arguments.
                                          */
        template <typename PFun, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<10> &)
          {
            (*fun_ptr) (arg_list.template get<0>(),
                        arg_list.template get<1>(),
                        arg_list.template get<2>(),
                        arg_list.template get<3>(),
                        arg_list.template get<4>(),
                        arg_list.template get<5>(),
                        arg_list.template get<6>(),
                        arg_list.template get<7>(),
                        arg_list.template get<8>(),
                        arg_list.template get<9>());
          }

                                         /**
                                          * Call a void member function
                                          * with 10 arguments.
                                          */
        template <typename PFun, typename C, typename ArgList>
        static inline void do_call (PFun     fun_ptr,
                                    C       &obj,
                                    ArgList &arg_list,
                                    internal::return_value<void> &,
                                    const int2type<10> &)
          {
            (obj.*fun_ptr) (arg_list.template get<0>(),
                            arg_list.template get<1>(),
                            arg_list.template get<2>(),
                            arg_list.template get<3>(),
                            arg_list.template get<4>(),
                            arg_list.template get<5>(),
                            arg_list.template get<6>(),
                            arg_list.template get<7>(),
                            arg_list.template get<8>(),
                            arg_list.template get<9>());
          }
    };
  
  

                                     /**
                                      * Call an arbitrary function by
                                      * dispatching to the functions in
                                      * the Caller class based on the
                                      * number of elements in the
                                      * argument list and the return
                                      * type.
                                      */
    template <typename RT, typename PFun, typename ArgList>
    static inline void call (PFun     fun_ptr,
                             ArgList &arg_list,
                             internal::return_value<RT> &ret_val)
    {
      Caller<RT>::do_call (fun_ptr, arg_list, ret_val,
                           int2type<boost::tuples::length<ArgList>::value>());
    }


  
                                     /**
                                      * Call an arbitrary member
                                      * function by dispatching to the
                                      * functions in the Caller class
                                      * based on the number of elements
                                      * in the argument list and the
                                      * return type.
                                      */
    template <typename RT, typename PFun, typename C, typename ArgList>
    static inline void call (PFun     fun_ptr,
                             C       &obj,
                             ArgList &arg_list,
                             internal::return_value<RT> &ret_val)
    {
      Caller<RT>::do_call (fun_ptr, obj, arg_list, ret_val,
                           int2type<boost::tuples::length<ArgList>::value>());
    }  
  }



  namespace internal
  {
                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments, and whether the
                                      * second argument is a const or
                                      * non-const class, dependening on
                                      * which the member function will
                                      * also me const or
                                      * non-const. There are
                                      * specializations of this class
                                      * for each number of arguments,
                                      * and the const and non-const
                                      * versions.
                                      */
    template <typename RT, class C, typename ArgList,
              int length = boost::tuples::length<ArgList>::value>
    struct mem_fun_ptr_helper;
  

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 0 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 0>
    {
        typedef RT (C::*type) ();
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 0 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 0>
    {
        typedef RT (C::*type) () const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 1 argument
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 1>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 1 argument
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 1>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 2 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 2>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 2 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 2>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 3 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 3>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 3 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 3>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 4 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 4>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 4 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 4>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 5 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 5>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 5 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 5>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 6 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 6>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 6 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 6>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 7 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 7>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 7 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 7>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 8 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 8>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type,
                               typename boost::tuples::element<7,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 8 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 8>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type,
                               typename boost::tuples::element<7,ArgList>::type) const;
    };


                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 9 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 9>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type,
                               typename boost::tuples::element<7,ArgList>::type,
                               typename boost::tuples::element<8,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 9 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 9>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type,
                               typename boost::tuples::element<7,ArgList>::type,
                               typename boost::tuples::element<8,ArgList>::type) const;
    };



                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 10 arguments
                                      * and non-const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, C, ArgList, 10>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type,
                               typename boost::tuples::element<7,ArgList>::type,
                               typename boost::tuples::element<8,ArgList>::type,
                               typename boost::tuples::element<9,ArgList>::type);
    };

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 10 arguments
                                      * and const member functions.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr_helper<RT, const C, ArgList, 10>
    {
        typedef RT (C::*type) (typename boost::tuples::element<0,ArgList>::type,
                               typename boost::tuples::element<1,ArgList>::type,
                               typename boost::tuples::element<2,ArgList>::type,
                               typename boost::tuples::element<3,ArgList>::type,
                               typename boost::tuples::element<4,ArgList>::type,
                               typename boost::tuples::element<5,ArgList>::type,
                               typename boost::tuples::element<6,ArgList>::type,
                               typename boost::tuples::element<7,ArgList>::type,
                               typename boost::tuples::element<8,ArgList>::type,
                               typename boost::tuples::element<9,ArgList>::type) const;
    };

  

                                     /**
                                      * Construct a pointer to member
                                      * function based on the template
                                      * arguments, and whether the
                                      * second argument is a const or
                                      * non-const class, dependening on
                                      * which the member function will
                                      * also me const or non-const. We
                                      * do this by dispatching to the
                                      * mem_fun_ptr_helper classes that
                                      * are overloaded on the number of
                                      * elements and the const/non-const
                                      * decision.
                                      *
                                      * Note that the last template
                                      * argument for the
                                      * mem_fun_ptr_helper class is
                                      * automatically computed in the
                                      * default argument to the general
                                      * template.
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_ptr
    {
        typedef typename mem_fun_ptr_helper<RT,C,ArgList>::type type;
    };  
  }


  namespace internal
  {
                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments, and whether the
                                      * second argument is a const or
                                      * non-const class, dependening on
                                      * which the member function will
                                      * also me const or
                                      * non-const. There are
                                      * specializations of this class
                                      * for each number of arguments.
                                      */
    template <typename RT, typename ArgList,
              int length = boost::tuples::length<ArgList>::value>
    struct fun_ptr_helper;
  

                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 0 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 0>
    {
        typedef RT (*type) ();
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 1 argument.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 1>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 2 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 2>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 3 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 3>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 4 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 4>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 5 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 5>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type,
                            typename boost::tuples::element<4,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 6 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 6>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type,
                            typename boost::tuples::element<4,ArgList>::type,
                            typename boost::tuples::element<5,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 7 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 7>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type,
                            typename boost::tuples::element<4,ArgList>::type,
                            typename boost::tuples::element<5,ArgList>::type,
                            typename boost::tuples::element<6,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 8 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 8>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type,
                            typename boost::tuples::element<4,ArgList>::type,
                            typename boost::tuples::element<5,ArgList>::type,
                            typename boost::tuples::element<6,ArgList>::type,
                            typename boost::tuples::element<7,ArgList>::type);
    };


                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 9 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 9>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type,
                            typename boost::tuples::element<4,ArgList>::type,
                            typename boost::tuples::element<5,ArgList>::type,
                            typename boost::tuples::element<6,ArgList>::type,
                            typename boost::tuples::element<7,ArgList>::type,
                            typename boost::tuples::element<8,ArgList>::type);
    };



                                     /**
                                      * Construct a pointer to non-member
                                      * function based on the template
                                      * arguments. This is the
                                      * specialization for 10 arguments.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 10>
    {
        typedef RT (*type) (typename boost::tuples::element<0,ArgList>::type,
                            typename boost::tuples::element<1,ArgList>::type,
                            typename boost::tuples::element<2,ArgList>::type,
                            typename boost::tuples::element<3,ArgList>::type,
                            typename boost::tuples::element<4,ArgList>::type,
                            typename boost::tuples::element<5,ArgList>::type,
                            typename boost::tuples::element<6,ArgList>::type,
                            typename boost::tuples::element<7,ArgList>::type,
                            typename boost::tuples::element<8,ArgList>::type,
                            typename boost::tuples::element<9,ArgList>::type);
    };

  

                                     /**
                                      * Construct a pointer to
                                      * non-member function based on the
                                      * template arguments. We do this
                                      * by dispatching to the
                                      * fun_ptr_helper classes that are
                                      * overloaded on the number of
                                      * elements.
                                      *
                                      * Note that the last template
                                      * argument for the
                                      * fun_ptr_helper class is
                                      * automatically computed in the
                                      * default argument to the general
                                      * template.
                                      */
    template <typename RT, typename ArgList>
    struct fun_ptr
    {
        typedef typename fun_ptr_helper<RT,ArgList>::type type;
    };  
  }


  namespace internal 
  {
                                     /**
                                      * Extract the Nth element of the
                                      * type list and make it a
                                      * reference.
                                      */
    template <int N, typename Tuple>
    struct add_reference_to_Nth
    {
        typedef typename boost::tuples::element<N,Tuple>::type ArgType;
        typedef typename boost::add_reference<ArgType>::type type;
    };
  
                                     /**
                                      * Specializations of this template
                                      * declare a typedef to a tuple
                                      * type that has the same basic
                                      * types as the first template
                                      * argument, but all references
                                      * instead of values. The second
                                      * argument is used to distinguish
                                      * between the lengths of argument
                                      * lists. The default argument
                                      * makes it possible to omit this
                                      * length argument.
                                      */
    template <typename Tuple, int = boost::tuples::length<Tuple>::value>
    struct tie_args_helper;


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 0.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,0>
    {
        typedef Tuple type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 1.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,1>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 2.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,2>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 3.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,3>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 4.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,4>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 5.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,5>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type,
                     typename add_reference_to_Nth<4,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 6.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,6>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type,
                     typename add_reference_to_Nth<4,Tuple>::type,
                     typename add_reference_to_Nth<5,Tuple>::type>
        type;
    };



                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 7.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,7>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type,
                     typename add_reference_to_Nth<4,Tuple>::type,
                     typename add_reference_to_Nth<5,Tuple>::type,
                     typename add_reference_to_Nth<6,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 8.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,8>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type,
                     typename add_reference_to_Nth<4,Tuple>::type,
                     typename add_reference_to_Nth<5,Tuple>::type,
                     typename add_reference_to_Nth<6,Tuple>::type,
                     typename add_reference_to_Nth<7,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 9.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,9>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type,
                     typename add_reference_to_Nth<4,Tuple>::type,
                     typename add_reference_to_Nth<5,Tuple>::type,
                     typename add_reference_to_Nth<6,Tuple>::type,
                     typename add_reference_to_Nth<7,Tuple>::type,
                     typename add_reference_to_Nth<8,Tuple>::type>
        type;
    };


                                     /**
                                      * Make a tuple type of all
                                      * references out of the given
                                      * tuple. Specialization for tuple
                                      * of length 10.
                                      */
    template <typename Tuple>
    struct tie_args_helper<Tuple,10>
    {
        typedef 
        boost::tuple<typename add_reference_to_Nth<0,Tuple>::type,
                     typename add_reference_to_Nth<1,Tuple>::type,
                     typename add_reference_to_Nth<2,Tuple>::type,
                     typename add_reference_to_Nth<3,Tuple>::type,
                     typename add_reference_to_Nth<4,Tuple>::type,
                     typename add_reference_to_Nth<5,Tuple>::type,
                     typename add_reference_to_Nth<6,Tuple>::type,
                     typename add_reference_to_Nth<7,Tuple>::type,
                     typename add_reference_to_Nth<8,Tuple>::type,
                     typename add_reference_to_Nth<9,Tuple>::type>
        type;
    };


  
                                     /**
                                      * Declare a typedef to a tuple
                                      * type that has the same basic
                                      * types as the template
                                      * argument, but all references
                                      * instead of values.
                                      *
                                      * Do so by redirecting to the
                                      * tie_args_helper specializations;
                                      * note that the second argument of
                                      * these templates is computed
                                      * automatically by the default
                                      * argument specification.
                                      */
    template <typename Tuple>
    struct tie_args 
    {
        typedef typename tie_args_helper<Tuple>::type type;
    };
  }

#if (DEAL_II_USE_MT == 1)
#  if defined(DEAL_II_USE_MT_POSIX)
  
  namespace internal 
  {
                                     /**
                                      * Base class describing a
                                      * thread. This is the basic
                                      * class abstracting the
                                      * operating system's POSIX
                                      * implementation into a C++
                                      * class. It provides a mechanism
                                      * to start a new thread, as well
                                      * as for joining it.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    struct thread_description_base
    {
      private:
                                         /**
                                          * Variable holding the data
                                          * the operating system needs
                                          * to work with a thread.
                                          */
        pthread_t             thread;

                                         /**
                                          * Store whether the
                                          * @p{join()} member function
                                          * as already been called. If
                                          * @p{true}, then @p{join}
                                          * will return immediately,
                                          * otherwise it needs to go
                                          * through a call to the
                                          * operating system.
                                          *
                                          * This class is generated
                                          * exactly once per thread,
                                          * but is kept in the
                                          * background: user's should
                                          * not directly access this
                                          * class. Access to it is
                                          * performed through counted
                                          * pointers, both from
                                          * @p{Thread<>} objects as
                                          * well as from the thread
                                          * entry point function on
                                          * the new thread. It is only
                                          * deleted when all users
                                          * have freed it, which means
                                          * that also the new thread
                                          * has ended.
                                          */
        mutable volatile bool was_joined;

                                         /**
                                          * Mutex used to synchronise
                                          * calls to the @p{join()}
                                          * function.
                                          */
        mutable ThreadMutex   join_mutex;

      public:

                                         /**
                                          * Constructor.
                                          */
        thread_description_base () : was_joined (false) {};

                                         /**
                                          * Destructor.
                                          */
        virtual ~thread_description_base ();

                                         /**
                                          * Create a new thread with
                                          * the given thread entry
                                          * point and arguments. Store
                                          * the result of the
                                          * operation in the
                                          * @p{thread} member variable
                                          * for further use.
                                          */
        void create (void * (*p) (void *), void *d);

                                         /**
                                          * Join a thread, i.e. wait
                                          * for it to finish. This
                                          * function can safely be
                                          * called from different
                                          * threads at the same time,
                                          * and can also be called
                                          * more than once.
                                          */
        void join () const;
    };

#  else       // some other threading model
#    error Not Implemented
#  endif     // defined(DEAL_II_USE_MT_POSIX)

                                     /**
                                      * Class derived from
                                      * @ref{thread_description_base}
                                      * that also provides the
                                      * possibility to store a return
                                      * value.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    template <typename RT>
    struct thread_description : public thread_description_base
    {
        return_value<RT> ret_val;
    };

                                     // forward declare another class
    template <typename, typename> struct wrapper_base;
  }



                                   /**
                                    * User visible class describing a
                                    * thread. Relays all real calls to
                                    * the internal thread object
                                    * abstracting the operating
                                    * system's functions, to which it
                                    * keeps a shared pointer. This
                                    * object can be freely copied
                                    * around in user space.
                                    *
                                    * The default value of the
                                    * template argument is @p{void},
                                    * so if the function you are
                                    * calling on a new thread has no
                                    * return value, you can omit the
                                    * template argument.
                                    * 
                                    * @author Wolfgang Bangerth, 2003
                                    */
  template <typename RT = void>
  class Thread
  {
                                       /**
                                        * Construct a thread object
                                        * with a pointer to an
                                        * internal thread object. This
                                        * is the constructor used to
                                        * the @p{spawn} family of
                                        * functions.
                                        *
                                        * We would like to make this
                                        * constructor private and only
                                        * grant the
                                        * @p{wrapper_base::fire_up}
                                        * function friendship, but
                                        * granting friendship to
                                        * functions in other
                                        * namespaces doesn't work with
                                        * some compilers, so only do
                                        * so if the configure script
                                        * decided that this is safe.
                                        */
#if defined(DEAL_II_NAMESP_TEMPL_FRIEND_BUG2) || defined(DEAL_II_NAMESP_TEMPL_FRIEND_BUG)
    public:
#endif
      Thread (const boost::shared_ptr<internal::thread_description<RT> > &td)
                      : thread_descriptor (td) {};

    public:
      
                                       /**
                                        * Default constructor. You
                                        * can't do much with a thread
                                        * object constructed this way,
                                        * except for assigning it a
                                        * thread object that holds
                                        * data created by the
                                        * @p{spawn} functions.
                                        */
      Thread () {};

                                       /**
                                        * Join the thread represented
                                        * by this object, i.e. wait
                                        * for it to finish. You can't
                                        * call this function if you
                                        * have used the default
                                        * constructor of this class
                                        * and have not assigned a
                                        * thread object to it.
                                        */
      void join () const {
        AssertThrow (thread_descriptor, ExcNoThread());
        thread_descriptor->join ();
      };

                                       /**
                                        * Get the return value of the
                                        * function of the
                                        * thread. Since this is only
                                        * available once the thread
                                        * finishes, this implicitely
                                        * also calls @p{join()}.
                                        */
      RT return_value () {
        join ();
        return thread_descriptor->ret_val.get();
      };


                                       /**
                                        * Check for equality of thread
                                        * objects. Since objects of
                                        * this class store an implicit
                                        * pointer to an object that
                                        * exists exactly once for each
                                        * thread, the check is simply
                                        * to compare these pointers.
                                        */
      bool operator == (const Thread &t) {
        return thread_descriptor == t.thread_descriptor;
      };

                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcNoThread);
      
    private:
                                       /**
                                        * Shared pointer to the object
                                        * representing the thread, and
                                        * abstracting operating system
                                        * functions to work on
                                        * it. Boost's shared pointer
                                        * implementation will make
                                        * sure that that object lives
                                        * as long as there is at least
                                        * one subscriber to it.
                                        */
      boost::shared_ptr<internal::thread_description<RT> > thread_descriptor;

#if !defined(DEAL_II_NAMESP_TEMPL_FRIEND_BUG2) && !defined(DEAL_II_NAMESP_TEMPL_FRIEND_BUG)
      template <typename, typename> friend struct internal::wrapper_base;
#endif
  };



  namespace internal
  {

                                     /**
                                      * Base class for the classes
                                      * wrapping function pointers and
                                      * arguments for non-member and
                                      * member functions. The first
                                      * template class denotes the
                                      * return type of the function
                                      * being called. The second one
                                      * is a class that provides a
                                      * function @p{entry_point},
                                      * which will be used as an entry
                                      * point for the new thread. In
                                      * the classes derived from this
                                      * one, the second template
                                      * argument is actually the
                                      * derived class itself, using
                                      * something like the
                                      * Barton-Nackman trick.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    template <typename RT, typename EntryPointClass>
    struct wrapper_base
    {

                                         /**
                                          * Start a new thread, wait
                                          * until it has copied the
                                          * data out of this object,
                                          * and return the thread
                                          * descriptor.
                                          */
        Thread<RT> fire_up () {
          thread_descriptor =
            DescriptionPointer(new internal::thread_description<RT>());

          ThreadMutex::ScopedLock lock (mutex);        
          thread_descriptor->create (&EntryPointClass::entry_point,
				     (void *)this);
          condition.wait (mutex);

          return thread_descriptor;
        }

      protected:
                                         /**
                                          * Typedef for shared
                                          * pointers to the objects
                                          * describing threads on the
                                          * OS level.
                                          */
        typedef
        boost::shared_ptr<internal::thread_description<RT> >
        DescriptionPointer;

                                         /**
                                          * Shared pointer to the
                                          * unique object describing a
                                          * thread on the OS level.
                                          */
        DescriptionPointer thread_descriptor;

                                         /**
                                          * Mutex and condition
                                          * variable used to
                                          * synchronise calling and
                                          * called thread.
                                          */
        mutable ThreadMutex     mutex;    
        mutable ThreadCondition condition;
    };
  

                                     /**
                                      * Wrap the arguments to a
                                      * non-member or static member
                                      * function and provide an entry
                                      * point for a new thread that
                                      * unwraps this data and calls
                                      * the function with them.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    template <typename RT, typename ArgList>
    struct fun_wrapper : public wrapper_base<RT, fun_wrapper<RT,ArgList> >
    {
                                         /**
                                          * Typedef for the type of
                                          * the function to be called
                                          * on the new thread.
                                          */
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;

                                         /**
                                          * Typedef for a typelist of
                                          * reference arguments to the
                                          * function to be called.
                                          */
        typedef typename internal::tie_args<ArgList>::type ArgReferences;

                                         /**
                                          * Constructor. Store the
                                          * necessary information
                                          * about the function to be
                                          * called and with which
                                          * arguments.
                                          *
                                          * Pass down the address of
                                          * this class's thread entry
                                          * point function. This way,
                                          * we can ensure that object
                                          * and thread entry point
                                          * function always are in
                                          * synch with respect to
                                          * their knowledge of the
                                          * types involved.
                                          */
        fun_wrapper (FunPtr               fun_ptr,
                     const ArgReferences &args)
                        : fun_ptr (fun_ptr),
                          args (args)  {};
      private:
                                         /**
                                          * Default constructor. Made
                                          * private and not
                                          * implemented to prevent
                                          * calling.
                                          */
        fun_wrapper ();

                                         /**
                                          * Copy constructor. Made
                                          * private and not
                                          * implemented to prevent
                                          * calling.
                                          */
        fun_wrapper (const fun_wrapper &);

                                         /**
                                          * Pointer to the function to
                                          * be called on the new
                                          * thread.
                                          */
        FunPtr        fun_ptr;

                                         /**
                                          * References to the
                                          * arguments with which the
                                          * function is to be called.
                                          */
        ArgReferences args;
      

                                         /**
                                          * Entry point for the new
                                          * thread.
                                          */
        static void * entry_point (void *arg)
          {
            const wrapper_base<RT, fun_wrapper> *w
              = reinterpret_cast<const wrapper_base<RT, fun_wrapper>*> (arg);
            const fun_wrapper *wrapper
              = static_cast<const fun_wrapper*> (w);

                                             // copy information from
                                             // the stack of the
                                             // calling thread
            FunPtr    fun_ptr = wrapper->fun_ptr;
            ArgList   args    = wrapper->args;

            boost::shared_ptr<internal::thread_description<RT> >
              thread_descriptor  = wrapper->thread_descriptor;
          
                                             // signal the fact that
                                             // we have copied all the
                                             // information that is
                                             // needed
            {
              ThreadMutex::ScopedLock lock (wrapper->mutex);
              wrapper->condition.signal ();
            }

                                             // call the
                                             // function. since an
                                             // exception that is
                                             // thrown from one of the
                                             // called functions will
                                             // not propagate to the
                                             // main thread, it will
                                             // kill the program if
                                             // not treated here
                                             // before we return to
                                             // the operating system's
                                             // thread library
            internal::register_thread ();
            try 
              {
                internal::call (fun_ptr, args,
                                thread_descriptor->ret_val);
              }
            catch (const std::exception &exc)
              {
                internal::handle_std_exception (exc);
              }
            catch (...)
              {
                internal::handle_unknown_exception ();
              }
            internal::deregister_thread ();
          
            return 0;
          };


					 /**
					  * Make the base class a
					  * friend, so that it can
					  * access the thread entry
					  * point function.
					  */
	template <typename, typename> friend class wrapper_base;
    };

  
  
                                     /**
                                      * Wrap the arguments to a member
                                      * function and provide an entry
                                      * point for a new thread that
                                      * unwraps this data and calls
                                      * the function with them.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    template <typename RT, class C, typename ArgList>
    struct mem_fun_wrapper : public wrapper_base<RT, mem_fun_wrapper<RT,C,ArgList> >
    {
                                         /**
                                          * Typedef for the type of
                                          * the function to be called
                                          * on the new thread.
                                          */
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;

                                         /**
                                          * Typedef for a typelist of
                                          * reference arguments to the
                                          * function to be called.
                                          */
        typedef typename internal::tie_args<ArgList>::type ArgReferences;
      
                                         /**
                                          * Constructor. Store the
                                          * necessary information
                                          * about the function to be
                                          * called and with which
                                          * arguments.
                                          *
                                          * Pass down the address of
                                          * this class's thread entry
                                          * point function. This way,
                                          * we can ensure that object
                                          * and thread entry point
                                          * function always are in
                                          * synch with respect to
                                          * their knowledge of the
                                          * types involved.
                                          */
        mem_fun_wrapper (MemFunPtr            mem_fun_ptr,
                         C                   &c,
                         const ArgReferences &args)
                        :
			c (c),
			mem_fun_ptr (mem_fun_ptr),
			args (args)  {};
      private:
                                         /**
                                          * Default constructor. Made
                                          * private and not
                                          * implemented to prevent
                                          * calling.
                                          */
        mem_fun_wrapper ();

                                         /**
                                          * Copy constructor. Made
                                          * private and not
                                          * implemented to prevent
                                          * calling.
                                          */
        mem_fun_wrapper (const mem_fun_wrapper &);
      
                                         /**
                                          * Pointer to the function to
                                          * be called on the new
                                          * thread, as well as the
                                          * object with which this has
                                          * to happen.
                                          */
        C            &c;
        MemFunPtr     mem_fun_ptr;

                                         /**
                                          * References to the
                                          * arguments with which the
                                          * function is to be called.
                                          */
        ArgReferences args;

                                         /**
                                          * Entry point for the new
                                          * thread.
                                          */
        static void * entry_point (void *arg)
          {
            const wrapper_base<RT,mem_fun_wrapper> *w
              = reinterpret_cast<const wrapper_base<RT,mem_fun_wrapper>*> (arg);
            const mem_fun_wrapper *wrapper
              = static_cast<const mem_fun_wrapper*> (w);

                                             // copy information from
                                             // the stack of the
                                             // calling thread
            MemFunPtr mem_fun_ptr = wrapper->mem_fun_ptr;
            C        &c           = wrapper->c;
            ArgList   args        = wrapper->args;

            boost::shared_ptr<internal::thread_description<RT> >
              thread_descriptor  = wrapper->thread_descriptor;

                                             // signal the fact that
                                             // we have copied all the
                                             // information that is
                                             // needed
            {
              ThreadMutex::ScopedLock lock (wrapper->mutex);
              wrapper->condition.signal ();
            }
          
                                             // call the
                                             // function. since an
                                             // exception that is
                                             // thrown from one of the
                                             // called functions will
                                             // not propagate to the
                                             // main thread, it will
                                             // kill the program if
                                             // not treated here
                                             // before we return to
                                             // the operating system's
                                             // thread library
            internal::register_thread ();
            try 
              {
                internal::call (mem_fun_ptr, c, args,
                                thread_descriptor->ret_val);
              }
            catch (const std::exception &exc)
              {
                internal::handle_std_exception (exc);
              }
            catch (...)
              {
                internal::handle_unknown_exception ();
              }
            internal::deregister_thread ();
          
            return 0;
          };

					 /**
					  * Make the base class a
					  * friend, so that it can
					  * access the thread entry
					  * point function.
					  */
	template <typename, typename> friend class wrapper_base;
    };
  }


  namespace internal
  {
                                     /**
                                      * General template declaration
                                      * of a class that is used to
                                      * encapsulate arguments to
                                      * global and static member
                                      * functions, make sure a new
                                      * thread is created and that
                                      * function being run on that
                                      * thread.
                                      *
                                      * Although this general template
                                      * is not implemented at all, the
                                      * default template argument
                                      * makes sure that whenever using
                                      * the name of this class, the
                                      * last template argument will be
                                      * computed correctly from the
                                      * previous arguments, and the
                                      * correct specialization for
                                      * this last template argument be
                                      * used, even though we need to
                                      * specify it.
                                      */
    template <typename RT, typename ArgList,
              int length = boost::tuples::length<ArgList>::value>
    class fun_encapsulator;


                                     /**
                                      * General template declaration
                                      * of a class that is used to
                                      * encapsulate arguments to
                                      * non-static member functions,
                                      * make sure a new thread is
                                      * created and that function
                                      * being run on that thread.
                                      *
                                      * Although this general template
                                      * is not implemented at all, the
                                      * default template argument
                                      * makes sure that whenever using
                                      * the name of this class, the
                                      * last template argument will be
                                      * computed correctly from the
                                      * previous arguments, and the
                                      * correct specialization for
                                      * this last template argument be
                                      * used, even though we need to
                                      * specify it.
                                      */
    template <typename RT, typename C, typename ArgList,
              int length = boost::tuples::length<ArgList>::value>
    class mem_fun_encapsulator;
  }


// ----------- encapsulators for member functions not taking any parameters

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with no arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 0>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() () {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                ArgList()).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * no arguments.
                                    */
  template <typename RT, typename C>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<> >
  spawn (C &c, RT (C::*fun_ptr)()) {
    return internal::mem_fun_encapsulator<RT, C, boost::tuple<> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with no
                                    * arguments.
                                    */
  template <typename RT, typename C>
  inline
  internal::mem_fun_encapsulator<RT,const C,boost::tuple<> >
  spawn (const C &c, RT (C::*fun_ptr)() const) {
    return internal::mem_fun_encapsulator<RT, const C, boost::tuple<> > (c,fun_ptr);
  }




// ----------- encapsulators for unary member functions

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 1 argument.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 1>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 1 argument.
                                    */
  template <typename RT, typename C, typename Arg1>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1)) {
    return internal::mem_fun_encapsulator<RT, C, boost::tuple<Arg1> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 1
                                    * argument.
                                    */
  template <typename RT, typename C, typename Arg1>
  inline
  internal::mem_fun_encapsulator<RT,const C,boost::tuple<Arg1> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1) const) {
    return internal::mem_fun_encapsulator<RT, const C, boost::tuple<Arg1> > (c,fun_ptr);
  }




// ----------- encapsulators for binary member functions

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 2 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 2>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,
                                                           arg2)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 2 arguments.
                                    */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2)) {
    return internal::mem_fun_encapsulator<RT, C, boost::tuple<Arg1, Arg2> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 2
                                    * arguments.
                                    */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::mem_fun_encapsulator<RT,const C,boost::tuple<Arg1, Arg2> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2) const) {
    return internal::mem_fun_encapsulator<RT, const C, boost::tuple<Arg1, Arg2> > (c,fun_ptr);
  }



// ----------- encapsulators for ternary member functions

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 3 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 3>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,
                                                           arg2,
                                                           arg3)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 3 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 3
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::mem_fun_encapsulator<RT,const C,boost::tuple<Arg1, Arg2, Arg3> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3> > (c,fun_ptr);
  }




// ----------- encapsulators for member functions with 4 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 4 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 4>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 4 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3, Arg4> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3, Arg4> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 4
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::mem_fun_encapsulator<RT,const C,boost::tuple<Arg1, Arg2, Arg3, Arg4> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3, Arg4> > (c,fun_ptr);
  }




// ----------- encapsulators for member functions with 5 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 5 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 5>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4,
                                                           arg5)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 5 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 5
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::mem_fun_encapsulator<RT,const C,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5> > (c,fun_ptr);
  }


// ----------- encapsulators for member functions with 6 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 6 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 6>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4,
                                                           arg5,arg6)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 6 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 6
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::mem_fun_encapsulator<RT,const C,
                               boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6> > (c,fun_ptr);
  }


// ----------- encapsulators for member functions with 7 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 7 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 7>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4,
                                                           arg5,arg6,
                                                           arg7)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 7 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                 Arg4, Arg5, Arg6, Arg7> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 7
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::mem_fun_encapsulator<RT,const C,
                               boost::tuple<Arg1, Arg2, Arg3,
                                            Arg4, Arg5, Arg6, Arg7> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7> > (c,fun_ptr);
  }


// ----------- encapsulators for member functions with 8 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 8 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 8>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4,
                                                           arg5,arg6,
                                                           arg7,arg8)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 8 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                 Arg4, Arg5, Arg6,
                                                 Arg7, Arg8> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 8
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::mem_fun_encapsulator<RT,const C,
                               boost::tuple<Arg1, Arg2, Arg3,
                                            Arg4, Arg5, Arg6,
                                            Arg7, Arg8> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                      Arg6,Arg7,Arg8) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8> > (c,fun_ptr);
  }


// ----------- encapsulators for member functions with 9 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 9 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 9>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4,
                                                           arg5,arg6,
                                                           arg7,arg8,
                                                           arg9)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 9 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                 Arg4, Arg5, Arg6,
                                                 Arg7, Arg8, Arg9> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8,Arg9)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 9
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::mem_fun_encapsulator<RT,const C,
                               boost::tuple<Arg1, Arg2, Arg3,
                                            Arg4, Arg5, Arg6,
                                            Arg7, Arg8, Arg9> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                      Arg6,Arg7,Arg8,Arg9) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9> > (c,fun_ptr);
  }


// ----------- encapsulators for member functions with 10 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for member
                                      * functions with 10 arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_encapsulator<RT, C, ArgList, 10>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_encapsulator (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9,
                    typename boost::tuples::element<9,ArgList>::type arg10) {
          return mem_fun_wrapper<RT,C,ArgList> (mem_fun_ptr, c,
                                                boost::tie(arg1,arg2,
                                                           arg3,arg4,
                                                           arg5,arg6,
                                                           arg7,arg8,
                                                           arg9, arg10)).fire_up ();
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 10 arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9,
            typename Arg10>
  inline
  internal::mem_fun_encapsulator<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                 Arg4, Arg5, Arg6,
                                                 Arg7, Arg8, Arg9, Arg10> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8,Arg9,Arg10)) {
    return internal::mem_fun_encapsulator<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9,
      Arg10> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 10
                                    * arguments.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9,
            typename Arg10>
  inline
  internal::mem_fun_encapsulator<RT,const C,
                               boost::tuple<Arg1, Arg2, Arg3,
                                            Arg4, Arg5, Arg6,
                                            Arg7, Arg8, Arg9, Arg10> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                      Arg6,Arg7,Arg8,Arg9,Arg10) const) {
    return internal::mem_fun_encapsulator<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9,
      Arg10> > (c,fun_ptr);
  }



// ----------- encapsulators for functions not taking any parameters

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with no arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 0>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() () {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          ArgList()).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with no arguments.
                                    */
  template <typename RT>
  inline
  internal::fun_encapsulator<RT,boost::tuple<> >
  spawn (RT (*fun_ptr)()) {
    return internal::fun_encapsulator<RT, boost::tuple<> > (fun_ptr);
  }




// ----------- encapsulators for unary functions

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 1 argument.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 1>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 1 argument.
                                    */
  template <typename RT, typename Arg1>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1> >
  spawn (RT (*fun_ptr)(Arg1)) {
    return internal::fun_encapsulator<RT, boost::tuple<Arg1> > (fun_ptr);
  }




// ----------- encapsulators for binary functions

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 2 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 2>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,
                                                     arg2)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 2 arguments.
                                    */
  template <typename RT, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2> >
  spawn (RT (*fun_ptr)(Arg1,Arg2)) {
    return internal::fun_encapsulator<RT, boost::tuple<Arg1, Arg2> > (fun_ptr);
  }



// ----------- encapsulators for ternary functions

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 3 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 3>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,
                                                     arg2,
                                                     arg3)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 3 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3> > (fun_ptr);
  }




// ----------- encapsulators for functions with 4 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 4 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 4>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 4 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3, Arg4> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3, Arg4> > (fun_ptr);
  }




// ----------- encapsulators for functions with 5 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 5 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 5>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4,
                                                     arg5)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 5 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5> > (fun_ptr);
  }


// ----------- encapsulators for functions with 6 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 6 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 6>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4,
                                                     arg5,arg6)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 6 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6> > (fun_ptr);
  }


// ----------- encapsulators for functions with 7 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 7 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 7>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4,
                                                     arg5,arg6,
                                                     arg7)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 7 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6, Arg7> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7> > (fun_ptr);
  }


// ----------- encapsulators for functions with 8 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 8 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 8>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4,
                                                     arg5,arg6,
                                                     arg7,arg8)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 8 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6,
                                           Arg7, Arg8> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                       Arg6,Arg7,Arg8)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8> > (fun_ptr);
  }


// ----------- encapsulators for functions with 9 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 9 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 9>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4,
                                                     arg5,arg6,
                                                     arg7,arg8,
                                                     arg9)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 9 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6,
                                           Arg7, Arg8, Arg9> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                       Arg6,Arg7,Arg8,Arg9)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9> > (fun_ptr);
  }


// ----------- encapsulators for functions with 10 arguments

  namespace internal
  {
                                     /**
                                      * Encapsulator class for
                                      * functions with 10 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 10>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_encapsulator (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9,
                    typename boost::tuples::element<9,ArgList>::type arg10) {
          return fun_wrapper<RT,ArgList> (fun_ptr,
                                          boost::tie(arg1,arg2,
                                                     arg3,arg4,
                                                     arg5,arg6,
                                                     arg7,arg8,
                                                     arg9, arg10)).fire_up ();
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function for
                                    * non-member or static member
                                    * functions with 10 arguments.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9,
            typename Arg10>
  inline
  internal::fun_encapsulator<RT,boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6,
                                           Arg7, Arg8, Arg9, Arg10> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                       Arg6,Arg7,Arg8,Arg9,Arg10)) {
    return internal::fun_encapsulator<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9,
      Arg10> > (fun_ptr);
  }

#else  // #if (DEAL_II_USE_MT == 1)

  template <typename RT = void>
  class Thread 
  {
    public:
                                       /**
                                        * Default constructor.
                                        */
      Thread ()  {};

                                       /**
                                        * Initialize the return value
                                        * of this object using the
                                        * given member-function
                                        * pointer, object, and
                                        * argument list.
                                        */
      template <typename PFun, typename C, typename ArgRefs>
      Thread (PFun     fun_ptr,
              C       &obj,
              ArgRefs  arg_refs)
        {
          internal::call (fun_ptr, obj, arg_refs, rv);
        }

                                       /**
                                        * Initialize the return value
                                        * of this object using the
                                        * given function pointer, and
                                        * argument list.
                                        */
      template <typename PFun, typename ArgRefs>
      Thread (PFun    fun_ptr,
              ArgRefs arg_refs)
        {
          internal::call (fun_ptr, arg_refs, rv);
        }

                                       /**
                                        * Get the return value of the
                                        * function of the
                                        * thread.
                                        */
      RT return_value () const 
        {
          return rv.get();
        }

                                       /**
                                        * Join this thread. Is of
                                        * course a no-op in this of no
                                        * thread support.
                                        */
      void join () const {}

                                       /**
                                        * Compare for equality of
                                        * threads. Since thrheads are
                                        * not supported, there can
                                        * only be exactly one thread,
                                        * and the result is
                                        * @p{true}. This function
                                        * doesn't make much sense,
                                        * though, when threads are not
                                        * supported.
                                        */
      bool operator == (const Thread &)
        {
          return true;
        }
      
    private:
                                       /**
                                        * Store the return value of
                                        * the thread function.
                                        */
      internal::return_value<RT> rv;
  };
  

  namespace internal
  {
                                     /**
                                      * General template declaration
                                      * of a class that is used to
                                      * forward arguments to
                                      * global and static member
                                      * functions.
                                      *
                                      * Although this general template
                                      * is not implemented at all, the
                                      * default template argument
                                      * makes sure that whenever using
                                      * the name of this class, the
                                      * last template argument will be
                                      * computed correctly from the
                                      * previous arguments, and the
                                      * correct specialization for
                                      * this last template argument be
                                      * used, even though we need to
                                      * specify it.
                                      */
    template <typename RT, typename ArgList,
              int length = boost::tuples::length<ArgList>::value>
    class fun_forwarder;

                                     /**
                                      * General template declaration
                                      * of a class that is used to
                                      * forward arguments to
                                      * non-static member functions.
                                      *
                                      * Although this general template
                                      * is not implemented at all, the
                                      * default template argument
                                      * makes sure that whenever using
                                      * the name of this class, the
                                      * last template argument will be
                                      * computed correctly from the
                                      * previous arguments, and the
                                      * correct specialization for
                                      * this last template argument be
                                      * used, even though we need to
                                      * specify it.
                                      */
    template <typename RT, typename C, typename ArgList,
              int length = boost::tuples::length<ArgList>::value>
    class mem_fun_forwarder;    
  }


  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with no arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 0>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() () {
          return Thread<RT> (mem_fun_ptr, c,
                             ArgList());
        }
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with no arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<> >
  spawn (C &c, RT (C::*fun_ptr)()) {
    return internal::mem_fun_forwarder<RT, C, boost::tuple<> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * no arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C>
  inline
  internal::mem_fun_forwarder<RT,const C,boost::tuple<> >
  spawn (const C &c, RT (C::*fun_ptr)() const) {
    return internal::mem_fun_forwarder<RT, const C, boost::tuple<> > (c,fun_ptr);
  }




// ----------- forwarders for unary member functions

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 1 argument.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 1>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions with
                                    * 1 argument. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C, typename Arg1>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1)) {
    return internal::mem_fun_forwarder<RT, C, boost::tuple<Arg1> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function for
                                    * const member functions with 1
                                    * argument. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C, typename Arg1>
  inline
  internal::mem_fun_forwarder<RT,const C,boost::tuple<Arg1> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1) const) {
    return internal::mem_fun_forwarder<RT, const C, boost::tuple<Arg1> > (c,fun_ptr);
  }




// ----------- forwarders for binary member functions

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 2
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 2>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,
                                        arg2));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 2 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2)) {
    return internal::mem_fun_forwarder<RT, C, boost::tuple<Arg1, Arg2> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 2 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::mem_fun_forwarder<RT,const C,boost::tuple<Arg1, Arg2> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2) const) {
    return internal::mem_fun_forwarder<RT, const C, boost::tuple<Arg1, Arg2> > (c,fun_ptr);
  }



// ----------- forwarders for ternary member functions

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 3
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 3>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,
                                        arg2,
                                        arg3));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 3 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 3 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::mem_fun_forwarder<RT,const C,boost::tuple<Arg1, Arg2, Arg3> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3> > (c,fun_ptr);
  }




// ----------- forwarders for member functions with 4 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 4
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 4>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 4 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3, Arg4> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3, Arg4> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 4 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::mem_fun_forwarder<RT,const C,boost::tuple<Arg1, Arg2, Arg3, Arg4> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3, Arg4> > (c,fun_ptr);
  }




// ----------- forwarders for member functions with 5 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 5
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 5>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 5 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 5 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::mem_fun_forwarder<RT,const C,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5> > (c,fun_ptr);
  }


// ----------- forwarders for member functions with 6 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 6
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 6>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 6 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 6 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::mem_fun_forwarder<RT,const C,
                              boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6> > (c,fun_ptr);
  }


// ----------- forwarders for member functions with 7 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 7
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 7>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 7 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                Arg4, Arg5, Arg6, Arg7> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 7 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::mem_fun_forwarder<RT,const C,
                              boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6, Arg7> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7> > (c,fun_ptr);
  }


// ----------- forwarders for member functions with 8 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 8
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 8>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7,arg8));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 8 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                Arg4, Arg5, Arg6,
                                                Arg7, Arg8> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 8 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::mem_fun_forwarder<RT,const C,
                              boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6,
                                           Arg7, Arg8> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                      Arg6,Arg7,Arg8) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8> > (c,fun_ptr);
  }


// ----------- forwarders for member functions with 9 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 9
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 9>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7,arg8,
                                        arg9));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 9 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                Arg4, Arg5, Arg6,
                                                Arg7, Arg8, Arg9> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8,Arg9)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 9 arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::mem_fun_forwarder<RT,const C,
                              boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6,
                                           Arg7, Arg8, Arg9> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                      Arg6,Arg7,Arg8,Arg9) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9> > (c,fun_ptr);
  }


// ----------- forwarders for member functions with 10 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for member
                                      * functions with 10
                                      * arguments.
                                      */
    template <typename RT, typename C, typename ArgList>
    class mem_fun_forwarder<RT, C, ArgList, 10>
    {
        typedef typename internal::mem_fun_ptr<RT,C,ArgList>::type MemFunPtr;      

      public:
        inline mem_fun_forwarder (C &c, MemFunPtr mem_fun_ptr)
                        : c (c), mem_fun_ptr(mem_fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9,
                    typename boost::tuples::element<9,ArgList>::type arg10) {
          return Thread<RT> (mem_fun_ptr, c,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7,arg8,
                                        arg9, arg10));
        };
    
      private:
        C         &c;
        MemFunPtr  mem_fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the non-const spawn
                                    * function for member functions
                                    * with 10 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9,
            typename Arg10>
  inline
  internal::mem_fun_forwarder<RT,C,boost::tuple<Arg1, Arg2, Arg3,
                                                Arg4, Arg5, Arg6,
                                                Arg7, Arg8, Arg9, Arg10> >
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8,Arg9,Arg10)) {
    return internal::mem_fun_forwarder<RT, C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9,
      Arg10> > (c,fun_ptr);
  }

                                   /**
                                    * Overload of the spawn function
                                    * for const member functions with
                                    * 10 arguments. This is the
                                    * version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9,
            typename Arg10>
  inline
  internal::mem_fun_forwarder<RT,const C,
                              boost::tuple<Arg1, Arg2, Arg3,
                                           Arg4, Arg5, Arg6,
                                           Arg7, Arg8, Arg9, Arg10> >
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                      Arg6,Arg7,Arg8,Arg9,Arg10) const) {
    return internal::mem_fun_forwarder<RT, const C,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9,
      Arg10> > (c,fun_ptr);
  }



// ----------- forwarders for functions not taking any parameters

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with no arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 0>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() () {
          return Thread<RT> (fun_ptr,
                             ArgList());
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with no
                                    * arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT>
  inline
  internal::fun_forwarder<RT,boost::tuple<> >
  spawn (RT (*fun_ptr)()) {
    return internal::fun_forwarder<RT, boost::tuple<> > (fun_ptr);
  }




// ----------- forwarders for unary functions

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 1 argument.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 1>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }

 
                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 1 argument. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename Arg1>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1> >
  spawn (RT (*fun_ptr)(Arg1)) {
    return internal::fun_forwarder<RT, boost::tuple<Arg1> > (fun_ptr);
  }




// ----------- forwarders for binary functions

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 2 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 2>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,
                                        arg2));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 2 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT, typename Arg1, typename Arg2>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2> >
  spawn (RT (*fun_ptr)(Arg1,Arg2)) {
    return internal::fun_forwarder<RT, boost::tuple<Arg1, Arg2> > (fun_ptr);
  }



// ----------- forwarders for ternary functions

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 3 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 3>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,
                                        arg2,
                                        arg3));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 3 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3> > (fun_ptr);
  }




// ----------- forwarders for functions with 4 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 4 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 4>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 4 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3, Arg4> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3, Arg4> > (fun_ptr);
  }




// ----------- forwarders for functions with 5 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 5 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 5>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 5 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5> > (fun_ptr);
  }


// ----------- forwarders for functions with 6 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 6 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 6>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 6 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6> > (fun_ptr);
  }


// ----------- forwarders for functions with 7 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 7 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 7>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 7 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3,
                                          Arg4, Arg5, Arg6, Arg7> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7> > (fun_ptr);
  }


// ----------- forwarders for functions with 8 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 8 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 8>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7,arg8));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 8 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3,
                                          Arg4, Arg5, Arg6,
                                          Arg7, Arg8> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                       Arg6,Arg7,Arg8)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8> > (fun_ptr);
  }


// ----------- forwarders for functions with 9 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 9 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 9>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7,arg8,
                                        arg9));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 9 arguments. This
                                    * is the version of the @p{spawn}
                                    * function for the case that
                                    * threading is not enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3,
                                          Arg4, Arg5, Arg6,
                                          Arg7, Arg8, Arg9> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                       Arg6,Arg7,Arg8,Arg9)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9> > (fun_ptr);
  }


// ----------- forwarders for functions with 10 arguments

  namespace internal
  {
                                     /**
                                      * Forwarder class for functions
                                      * with 10 arguments.
                                      */
    template <typename RT, typename ArgList>
    class fun_forwarder<RT, ArgList, 10>
    {
        typedef typename internal::fun_ptr<RT,ArgList>::type FunPtr;      

      public:
        inline fun_forwarder (FunPtr fun_ptr)
                        : fun_ptr(fun_ptr) {};

        inline
        Thread<RT>
        operator() (typename boost::tuples::element<0,ArgList>::type arg1,
                    typename boost::tuples::element<1,ArgList>::type arg2,
                    typename boost::tuples::element<2,ArgList>::type arg3,
                    typename boost::tuples::element<3,ArgList>::type arg4,
                    typename boost::tuples::element<4,ArgList>::type arg5,
                    typename boost::tuples::element<5,ArgList>::type arg6,
                    typename boost::tuples::element<6,ArgList>::type arg7,
                    typename boost::tuples::element<7,ArgList>::type arg8,
                    typename boost::tuples::element<8,ArgList>::type arg9,
                    typename boost::tuples::element<9,ArgList>::type arg10) {
          return Thread<RT> (fun_ptr,
                             boost::tie(arg1,arg2,
                                        arg3,arg4,
                                        arg5,arg6,
                                        arg7,arg8,
                                        arg9, arg10));
        };
    
      private:
        FunPtr  fun_ptr;
    };
  
  }


                                   /**
                                    * Overload of the spawn function
                                    * for non-member or static member
                                    * functions with 10
                                    * arguments. This is the version
                                    * of the @p{spawn} function for
                                    * the case that threading is not
                                    * enabled.
                                    */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9,
            typename Arg10>
  inline
  internal::fun_forwarder<RT,boost::tuple<Arg1, Arg2, Arg3,
                                          Arg4, Arg5, Arg6,
                                          Arg7, Arg8, Arg9, Arg10> >
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                       Arg6,Arg7,Arg8,Arg9,Arg10)) {
    return internal::fun_forwarder<RT,
      boost::tuple<Arg1, Arg2, Arg3,
      Arg4, Arg5, Arg6,
      Arg7, Arg8, Arg9,
      Arg10> > (fun_ptr);
  }
  

  
#endif // #if (DEAL_II_USE_MT == 1)



                                   /**
                                    * A container for thread
                                    * objects. Allows to add new
                                    * thread objects and wait for them
                                    * all together. The thread objects
                                    * need to have the same return
                                    * value for the called function.
                                    *
                                    * @author Wolfgang Bangerth, 2003
                                    */
  template <typename RT = void>
  class ThreadGroup 
  {
    public:
                                       /**
                                        * Add another thread object to
                                        * the collection.
                                        */
      ThreadGroup & operator += (const Thread<RT> &t) {
        threads.push_back (t);
        return *this;
      };

                                       /**
                                        * Wait for all threads in the
                                        * collection to finish. It is
                                        * not a problem if some of
                                        * them have already been
                                        * waited for, i.e. you may
                                        * call this function more than
                                        * once, and you can also add
                                        * new thread objects between
                                        * subsequent calls to this
                                        * function if you want.
                                        */
      void join_all () const {
        for (typename std::list<Thread<RT> >::const_iterator
               t=threads.begin(); t!=threads.end(); ++t)
          t->join ();
      };
    
    private:
                                       /**
                                        * List of thread objects.
                                        */
      std::list<Thread<RT> > threads;
  };
  
  
}   // end of implementation of namespace Threads




//----------------------------   thread_management.h     ---------------------------
// end of #ifndef __deal2__thread_management_h
#endif
//----------------------------   thread_management.h     ---------------------------

// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#ifndef __deal2__thread_management_h
#define __deal2__thread_management_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/std_cxx11/tuple.h>
#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/std_cxx11/bind.h>

#ifdef DEAL_II_WITH_THREADS
#  include <deal.II/base/std_cxx11/thread.h>
#  include <deal.II/base/std_cxx11/mutex.h>
#  include <deal.II/base/std_cxx11/condition_variable.h>
#endif

#include <iterator>
#include <vector>
#include <list>
#include <utility>


#ifdef DEAL_II_WITH_THREADS
#  ifdef DEAL_II_USE_MT_POSIX
#    include <pthread.h>
#  endif
#  include <tbb/task.h>
#endif



DEAL_II_NAMESPACE_OPEN

/*!@addtogroup threads */
/*@{*/


/**
 * A namespace for the implementation of thread management in
 * deal.II. Most of the content of this namespace is discussed in
 * detail in one of the reports linked to from the documentation page
 * of deal.II.
 *
 * @ingroup threads
 */
namespace Threads
{
  /**
   * This class is used instead of a true lock class when not using
   * multithreading. It allows to write programs such that they start new
   * threads and/or lock objects in multithreading mode, and use dummy thread
   * management and synchronization classes instead when running in
   * single-thread mode. Specifically, the <tt>new_thread</tt> functions only
   * call the function but wait for it to return instead of running in on
   * another thread, and the mutices do nothing really. The only reason to
   * provide such a function is that the program can be compiled both in MT and
   * non-MT mode without difference.
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
     * <tt>acquire</tt> and <tt>release</tt>.
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
      ScopedLock (DummyThreadMutex &) {}

      /**
       * Destructor. Unlock the
       * mutex. Since this is a
       * dummy mutex class, this
       * of course does nothing.
       */
      ~ScopedLock () {}
    };

    /**
     * Simulate acquisition of the
     * mutex. As this class does
     * nothing really, this
     * function does nothing as
     * well.
     */
    inline void acquire () const {}

    /**
     * Simulate release of the
     * mutex. As this class does
     * nothing really, this
     * function does nothing as
     * well.
     */
    inline void release () const {}
  };



  /**
   * This class is used in single threaded mode instead of a class
   * implementing real condition variable semantics. It allows to write
   * programs such that they start new threads and/or lock objects in
   * multithreading mode, and use dummy thread management and
   * synchronisation classes instead when running in single-thread
   * mode. Specifically, the <tt>new_thread</tt> functions only call the function
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
    inline void signal () const {}

    /**
     * Signal to multiple listener
     * that a condition has been
     * met, i.e. that some data
     * will now be available. Since
     * in single threaded mode,
     * this function of course does
     * nothing.
     */
    inline void broadcast () const {}

    /**
     * Wait for the condition to be
     * signalled. Signal variables
     * need to be guarded by a
     * mutex which needs to be
     * given to this function as an
     * argument, see the man page
     * of <tt>posix_cond_wait</tt> for a
     * description of the
     * mechanisms. Since in single
     * threaded mode, this function
     * of course does nothing, but
     * returns immediately.
     */
    inline void wait (DummyThreadMutex &) const {}
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
   * the <tt>count</tt> argument denoting the number of threads that need to
   * be synchronized is not equal to one.
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
     * exception if the <tt>count</tt>
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
    int wait () const
    {
      return 0;
    }

    /**
     * Dump the state of this
     * object. Here: do nothing.
     */
    void dump () const {}

    /** @addtogroup Exceptions
     * @{ */

    /**
     * Exception.
     */
    DeclException1 (ExcBarrierSizeNotUseful,
                    int,
                    << "In single-thread mode, other barrier sizes than 1 are not "
                    << "useful. You gave " << arg1);

    //@}
  };


#ifdef DEAL_II_WITH_THREADS

  /**
   * Class implementing a
   * Mutex. Mutexes are used to lock
   * data structures to ensure that
   * only a single thread of
   * execution can access them at the
   * same time.
   *
   * <h3>Copy semantics</h3>
   *
   * When copied, the receiving
   * object does not receive any
   * state from the object being
   * copied, i.e. an entirely new
   * mutex is created. This is
   * consistent with expectations if
   * a mutex is used as a member
   * variable to lock the other
   * member variables of a class: in
   * that case, the mutex of the
   * copied-to object should only
   * guard the members of the
   * copied-to object, not the
   * members of both the copied-to
   * and copied-from object.
   *
   * @author Wolfgang Bangerth, 2002, 2003, 2009
   */
  class Mutex
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
     * <tt>acquire</tt> and <tt>release</tt>.
     */
    class ScopedLock
    {
    public:
      /**
       * Constructor. Lock the
       * mutex.
       */
      ScopedLock (Mutex &m) : mutex(m)
      {
        mutex.acquire();
      }

      /**
       * Destructor. Unlock the
       * mutex. Since this is a
       * dummy mutex class, this
       * of course does nothing.
       */
      ~ScopedLock ()
      {
        mutex.release ();
      }

    private:
      /**
       * Store the address of the
       * mutex object.
       */
      Mutex &mutex;
    };

    /**
     * Default constructor.
     */
    Mutex ()
    {}

    /**
     * Copy constructor. As
     * discussed in this class's
     * documentation, no state is
     * copied from the object given
     * as argument.
     */
    Mutex (const Mutex &)
      :
      mutex()
    {}


    /**
     * Acquire a mutex.
     */
    inline void acquire ()
    {
      mutex.lock();
    }

    /**
     * Release the mutex again.
     */
    inline void release ()
    {
      mutex.unlock();
    }

  private:
    /**
     * Data object storing the
     * mutex data
     */
    std_cxx11::mutex mutex;

    /**
     * Make the class implementing
     * condition variables a friend, since
     * it needs to access the mutex.
     */
    friend class ConditionVariable;
  };


  /**
   * Class implementing a condition
   * variable. The semantics of this
   * class and its member functions
   * are the same as those of the
   * POSIX functions.
   *
   * @author Wolfgang Bangerth, 2003
   */
  class ConditionVariable
  {
  public:
    /**
     * Signal to a single listener
     * that a condition has been
     * met, i.e. that some data
     * will now be available.
     */
    inline void signal ()
    {
      condition_variable.notify_one();
    }

    /**
     * Signal to multiple listener
     * that a condition has been
     * met, i.e. that some data
     * will now be available.
     */
    inline void broadcast ()
    {
      condition_variable.notify_all();
    }

    /**
     * Wait for the condition to be
     * signalled. Signal variables
     * need to be guarded by a
     * mutex which needs to be
     * given to this function as an
     * argument, see the man page
     * of <tt>pthread_cond_wait</tt> for a
     * description of the
     * mechanisms.
     *
     * The mutex is assumed held at the
     * entry to this function but is
     * released upon exit.
     */
    inline void wait (Mutex &mutex)
    {
      std_cxx11::unique_lock<std_cxx11::mutex> lock(mutex.mutex,
                                                    std_cxx11::adopt_lock);
      condition_variable.wait (lock);
    }

  private:
    /**
     * Data object storing the
     * necessary data.
     */
    std_cxx11::condition_variable condition_variable;
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
   * Provide a backward compatible name (we
   * used ThreadMutex up to release 6.1, but
   * since it is already in a namespace
   * Threads this seems redundant).
   *
   * @deprecated
   */
  typedef Mutex ThreadMutex DEAL_II_DEPRECATED;

  /**
   * Provide a backward compatible name (we
   * used ThreadCondition up to release 6.1,
   * but since it is already in a namespace
   * Threads this seems redundant).
   *
   * @deprecated
   */
  typedef ConditionVariable ThreadCondition DEAL_II_DEPRECATED;

  /**
   * If using POSIX functions, then
   * alias the POSIX wrapper classes
   * to the names we use throughout
   * the library.
   */
  typedef PosixThreadBarrier Barrier;

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
  typedef DummyThreadMutex     Mutex;
  typedef DummyThreadMutex     ThreadMutex;

  /**
   * In non-multithread mode, the
   * mutex and thread management
   * classes are aliased to dummy
   * classes that actually do
   * nothing, in particular not lock
   * objects. Likewise for the
   * barrier class.
   */
  typedef DummyThreadCondition ConditionVariable;
  typedef DummyThreadCondition ThreadCondition;

  /**
   * In non-multithread mode, the
   * mutex and thread management
   * classes are aliased to dummy
   * classes that actually do
   * nothing, in particular not lock
   * objects. Likewise for the
   * barrier class.
   */
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
   * the Threads::new_thread
   * function. It is decreased once
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
   * (e.g. <tt>pthread_create</tt>
   * for the POSIX thread interface),
   * or if they are killed
   * (e.g. either through
   * <tt>pthread_exit</tt> from the
   * spawned thread, or
   * <tt>pthread_kill</tt> from
   * another thread), then these
   * events are not registered and
   * counted for the result of this
   * function. Likewise, threads that
   * the Threading Building Blocks
   * library may have created to work
   * on tasks created using the
   * Threads::new_task functions are
   * not counted.
   *
   * @ingroup threads
   */
  unsigned int n_existing_threads ();

  /**
   * Return a number used as id of
   * this thread. This number is
   * generated using the system call
   * <tt>getpid</tt>, or, if it
   * exists <tt>gettid</tt>. The
   * result of either is converted to
   * an integer and returned by this
   * function.
   *
   * @todo As of now, none of our
   * systems seems to support
   * <tt>gettid</tt>, so that part of
   * the code is untested yet.
   *
   * @ingroup threads
   */
  unsigned int this_thread_id ();

  /**
   * Split the range <tt>[begin,end)</tt>
   * into <tt>n_intervals</tt> subintervals
   * of equal size. The last interval
   * will be a little bit larger, if
   * the number of elements in the
   * whole range is not exactly
   * divisible by <tt>n_intervals</tt>. The
   * type of the iterators has to
   * fulfill the requirements of a
   * forward iterator,
   * i.e. <tt>operator++</tt> must be
   * available, and of course it must
   * be assignable.
   *
   * A list of subintervals is
   * returned as a vector of pairs of
   * iterators, where each pair
   * denotes the range
   * <tt>[begin[i],end[i])</tt>.
   *
   * @ingroup threads
   */
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator,ForwardIterator> >
  split_range (const ForwardIterator &begin,
               const ForwardIterator &end,
               const unsigned int n_intervals);

  /**
   * Split the interval <tt>[begin,end)</tt>
   * into subintervals of (almost)
   * equal size. This function works
   * mostly as the one before, with
   * the difference that instead of
   * iterators, now values are taken
   * that define the whole interval.
   *
   * @ingroup threads
   */
  std::vector<std::pair<unsigned int,unsigned int> >
  split_interval (const unsigned int begin,
                  const unsigned int end,
                  const unsigned int n_intervals);

  /**
   * @cond internal
   */

  /**
   * A namespace in which helper
   * functions and the like for the
   * threading subsystem are
   * implemented. The members of this
   * namespace are not meant for
   * public use.
   *
   * @author Wolfgang Bangerth, 2003
   */
  namespace internal
  {
    /**
     * @internal
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
     * type <tt>std::exception</tt> is
     * caught, it passes over control
     * to this function, which will
     * then provide some output.
     */
    void handle_std_exception (const std::exception &exc);

    /**
     * @internal
     * Same as above, but the type of
     * the exception is not derived
     * from <tt>std::exception</tt>, so
     * there is little way to provide
     * something more useful.
     */
    void handle_unknown_exception ();

    /**
     * @internal
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
     * @internal
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

  /**
   * @endcond
   */

}   // end declarations of namespace Threads

/* ----------- implementation of functions in namespace Threads ---------- */
#ifndef DOXYGEN
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
    const unsigned int n_elements              = std::distance (begin, end);
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
}

#endif // DOXYGEN

namespace Threads
{
  namespace internal
  {
    /**
     * @internal
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

      inline RT get () const
      {
        return value;
      }
      inline void set (RT v)
      {
        value = v;
      }
    };


    /**
     * @internal
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
      RT *value;
    public:
      inline return_value () : value(0) {}

      inline RT &get () const
      {
        return *value;
      }
      inline void set (RT &v)
      {
        value = &v;
      }
    };


    /**
     * @internal
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
    template <> struct return_value<void>
    {
      static inline void get () {}
    };
  }



  namespace internal
  {
    template <typename RT>
    inline void call (const std_cxx11::function<RT ()> &function,
                      internal::return_value<RT> &ret_val)
    {
      ret_val.set (function());
    }


    inline void call (const std_cxx11::function<void ()> &function,
                      internal::return_value<void> &)
    {
      function();
    }
  }



  namespace internal
  {
    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments, and whether the
     * second argument is a const or
     * non-const class, depending on
     * which the member function will
     * also me const or
     * non-const. There are
     * specializations of this class
     * for each number of arguments.
     */
    template <typename RT, typename ArgList,
              int length = std_cxx11::tuple_size<ArgList>::value>
    struct fun_ptr_helper;


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 0 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 0>
    {
      typedef RT (type) ();
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 1 argument.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 1>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 2 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 2>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 3 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 3>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 4 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 4>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 5 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 5>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type,
                         typename std_cxx11::tuple_element<4,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 6 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 6>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type,
                         typename std_cxx11::tuple_element<4,ArgList>::type,
                         typename std_cxx11::tuple_element<5,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 7 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 7>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type,
                         typename std_cxx11::tuple_element<4,ArgList>::type,
                         typename std_cxx11::tuple_element<5,ArgList>::type,
                         typename std_cxx11::tuple_element<6,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 8 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 8>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type,
                         typename std_cxx11::tuple_element<4,ArgList>::type,
                         typename std_cxx11::tuple_element<5,ArgList>::type,
                         typename std_cxx11::tuple_element<6,ArgList>::type,
                         typename std_cxx11::tuple_element<7,ArgList>::type);
    };


    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 9 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 9>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type,
                         typename std_cxx11::tuple_element<4,ArgList>::type,
                         typename std_cxx11::tuple_element<5,ArgList>::type,
                         typename std_cxx11::tuple_element<6,ArgList>::type,
                         typename std_cxx11::tuple_element<7,ArgList>::type,
                         typename std_cxx11::tuple_element<8,ArgList>::type);
    };



    /**
     * @internal
     * Construct a pointer to non-member
     * function based on the template
     * arguments. This is the
     * specialization for 10 arguments.
     */
    template <typename RT, typename ArgList>
    struct fun_ptr_helper<RT, ArgList, 10>
    {
      typedef RT (type) (typename std_cxx11::tuple_element<0,ArgList>::type,
                         typename std_cxx11::tuple_element<1,ArgList>::type,
                         typename std_cxx11::tuple_element<2,ArgList>::type,
                         typename std_cxx11::tuple_element<3,ArgList>::type,
                         typename std_cxx11::tuple_element<4,ArgList>::type,
                         typename std_cxx11::tuple_element<5,ArgList>::type,
                         typename std_cxx11::tuple_element<6,ArgList>::type,
                         typename std_cxx11::tuple_element<7,ArgList>::type,
                         typename std_cxx11::tuple_element<8,ArgList>::type,
                         typename std_cxx11::tuple_element<9,ArgList>::type);
    };



    /**
     * @internal
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
#ifdef DEAL_II_WITH_THREADS

    /**
     * A class that represents threads. For
     * each thread, we create exactly one of
     * these objects -- exactly one because
     * it carries the returned value of the
     * function called on the thread.
     *
     * While we have only one of these
     * objects per thread, several
     * Threads::Thread objects may refer to
     * this descriptor. If all Thread
     * objects go out of scope the
     * ThreadDescriptor will detach from
     * the thread before being destructed.
     */
    template <typename RT>
    struct ThreadDescriptor
    {
      /**
       * An object that represents the
       * thread started.
       */
      std_cxx11::thread thread;

      /**
       * An object that will hold the value
       * returned by the function called on
       * the thread.
       *
       * The return value is stored
       * in a shared_ptr because we might
       * abandon the the ThreadDescriptor.
       * This makes sure the object stays
       * alive until the thread exits.
       */
      std_cxx11::shared_ptr<return_value<RT> > ret_val;

      /**
       * A bool variable that is initially
       * false, is set to true when a new
       * thread is started, and is set back
       * to false once join() has been
       * called.
       *
       * We use this variable to make sure
       * we can call join() twice on the
       * same thread. For some reason, the
       * C++ standard library throws a
       * std::system_error exception if one
       * tries to call std::thread::join
       * twice (and in fact, before the
       * second call, std::thread::joinable
       * returns false) but this is a
       * somewhat desirable thing to do
       * because one doesn't have to keep
       * track whether join() has been
       * called before. Using this
       * variable, whenever we have called
       * join() before, the variable is set
       * to true and we can skip over
       * calling std::thread::join() a
       * second time. Access to this variable
       * is guarded by the following mutex.
       *
       * @note Historically, we did not need the
       * mutex for this variable: threads
       * can only be joined from the thread
       * that created it
       * originally. Consequently,
       * everything that happens in a
       * function that does not create
       * threads (such as the join()
       * function below) looks atomic to
       * the outside world. Since we clear
       * and test thread_is_active in the
       * same function as we call
       * std::thread::join, these actions
       * are atomic and need no mutex. Of
       * course, two threads may call
       * join() on the same thread object
       * at the same time, but this action
       * is undefined anyway since they can
       * not both join the same thread. That said,
       * more recent C++ standards do not appear to
       * have the requirement any more that the only
       * thread that can call join() is the one that
       * created the thread. Neither does pthread_join
       * appear to have this requirement any more.
       * Consequently, we can in fact join from different
       * threads and we test this in base/thread_validity_07.
       */
      bool thread_is_active;

      /**
       * Mutex guarding access to the previous variable.
       */
      Mutex thread_is_active_mutex;

      /**
       * Default constructor.
       */
      ThreadDescriptor ()
        :
        thread_is_active (false)
      {}

      ~ThreadDescriptor ()
      {
        if (!thread_is_active)
          return;
        thread.detach();
        thread_is_active = false;
      }

      /**
       * Start the thread and let
       * it put its return value
       * into the ret_val object.
       */
      void start (const std_cxx11::function<RT ()> &function)
      {
        thread_is_active = true;
        ret_val.reset(new return_value<RT>());
        thread = std_cxx11::thread (thread_entry_point, function, ret_val);
      }


      /**
       * Wait for the thread to end.
       */
      void join ()
      {
        // see if the thread hasn't been
        // joined yet. if it has, then
        // join() is a no-op. use schmidt's double-checking
        // strategy to use the mutex only when necessary
        if (thread_is_active == false)
          return;

        Mutex::ScopedLock lock (thread_is_active_mutex);
        if (thread_is_active == true)
          {
            Assert (thread.joinable(), ExcInternalError());
            thread.join ();
            thread_is_active = false;
          }
      }

    private:

      /**
       * The function that runs on the
       * thread.
       */
      static
      void thread_entry_point (const std_cxx11::function<RT ()> function,
                               std_cxx11::shared_ptr<return_value<RT> > ret_val)
      {
        // call the function
        // in question. since an
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
            call (function, *ret_val);
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
      }
    };

#else
    /**
     * A class that represents threads. For
     * each thread, we create exactly one of
     * these objects -- exactly one because
     * it carries the returned value of the
     * function called on the thread.
     *
     * While we have only one of these
     * objects per thread, several
     * Threads::Thread objects may refer to
     * this descriptor.
     */
    template <typename RT>
    struct ThreadDescriptor
    {
      /**
       * An object that will hold the value
       * returned by the function called on
       * the thread.
       */
      std_cxx11::shared_ptr<return_value<RT> > ret_val;

      /**
       * Start the thread and
       * let it put its return value into
       * the ret_val object.
       */
      void start (const std_cxx11::function<RT ()> &function)
      {
        ret_val.reset(new return_value<RT>());
        call (function, *ret_val);
      }

      /**
       * Wait for the thread to end.
       */
      void join ()
      {}
    };

#endif
  }


  /**
   * An object that represents a
   * spawned thread. This object can
   * be freely copied around in user
   * space, and all instances will
   * represent the same thread and
   * can require to wait for its
   * termination and access its
   * return value.
   *
   * Threads can be abandoned,
   * i.e. if you just call
   * Threads::new_thread but don't
   * care about the returned object,
   * or if you assign the return
   * Threads::Thread object to an
   * object that subsequently goes
   * out of scope, then the thread
   * previously created will still
   * continue to do work. You will
   * simply not be able to access its
   * return value any more, and it
   * may also happen that your
   * program terminates before the
   * thread has finished its work.
   *
   * The default value of the
   * template argument is <tt>void</tt>,
   * so if the function you are
   * calling on a new thread has no
   * return value, you can omit the
   * template argument.
   *
   * @author Wolfgang Bangerth, 2003, 2009
   * @ingroup threads
   * @ingroup threads
   */
  template <typename RT = void>
  class Thread
  {
  public:
    /**
     * Construct a thread object
     * with a function object.
     */
    Thread (const std_cxx11::function<RT ()> &function)
      :
      thread_descriptor (new internal::ThreadDescriptor<RT>())
    {
      // in a second step, start
      // the thread.
      thread_descriptor->start (function);
    }

    /**
     * Default constructor. You
     * can't do much with a thread
     * object constructed this way,
     * except for assigning it a
     * thread object that holds
     * data created by the
     * <tt>new_thread</tt> functions.
     */
    Thread () {}

    /**
     * Copy constructor.
     */
    Thread (const Thread<RT> &t)
      :
      thread_descriptor (t.thread_descriptor)
    {}

    /**
     * Join the thread represented by this
     * object, i.e. wait for it to
     * finish. If you have used the default
     * constructor of this class and have
     * not assigned a thread object to it,
     * then this function is a no-op.
     */
    void join () const
    {
      if (thread_descriptor)
        thread_descriptor->join ();
    }

    /**
     * Get the return value of the
     * function of the
     * thread. Since this is only
     * available once the thread
     * finishes, this implicitly
     * also calls <tt>join()</tt>.
     */
    RT return_value ()
    {
      join ();
      return thread_descriptor->ret_val->get();
    }

    /**
     * Return true if this object
     * has had a thread associated
     * with it, either by using the
     * non-default constructor or
     * by assignment.
     */
    bool valid () const
    {
      return static_cast<bool>(thread_descriptor);
    }


    /**
     * Check for equality of thread
     * objects. Since objects of
     * this class store an implicit
     * pointer to an object that
     * exists exactly once for each
     * thread, the check is simply
     * to compare these pointers.
     */
    bool operator == (const Thread &t)
    {
      return thread_descriptor == t.thread_descriptor;
    }

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
    std_cxx11::shared_ptr<internal::ThreadDescriptor<RT> > thread_descriptor;
  };


  namespace internal
  {
    /**
     * A general template that returns
     * std_cxx11::ref(t) if t is of reference
     * type, and t otherwise.
     *
     * The case that t is of reference type
     * is handled in a partial specialization
     * declared below.
     */
    template <typename T>
    struct maybe_make_ref
    {
      static T act (T &t)
      {
        return t;
      }
    };



    /**
     * A general template that returns
     * std_cxx11::ref(t) if t is of reference
     * type, and t otherwise.
     *
     * The case that t is of reference type
     * is handled in this partial
     * specialization.
     */
    template <typename T>
    struct maybe_make_ref<T &>
    {
      static std_cxx11::reference_wrapper<T> act (T &t)
      {
        return std_cxx11::ref(t);
      }
    };
  }



  namespace internal
  {
    /**
     * @internal
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
    template <typename RT, typename ArgList, int length>
    class fun_encapsulator;


// ----------- encapsulators for function objects

    /**
     * @internal
     * Encapsulator class for
     * functions with no arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 0>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() ()
      {
        return Thread<RT> (function);
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 1 argument.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 1>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 2 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 2>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 3 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 3>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 4 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 4>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3,
                  typename std_cxx11::tuple_element<3,ArgList>::type arg4)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<3,ArgList>::type>::act(arg4)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 5 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 5>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3,
                  typename std_cxx11::tuple_element<3,ArgList>::type arg4,
                  typename std_cxx11::tuple_element<4,ArgList>::type arg5)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<3,ArgList>::type>::act(arg4),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<4,ArgList>::type>::act(arg5)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 6 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 6>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3,
                  typename std_cxx11::tuple_element<3,ArgList>::type arg4,
                  typename std_cxx11::tuple_element<4,ArgList>::type arg5,
                  typename std_cxx11::tuple_element<5,ArgList>::type arg6)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<3,ArgList>::type>::act(arg4),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<4,ArgList>::type>::act(arg5),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<5,ArgList>::type>::act(arg6)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 7 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 7>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3,
                  typename std_cxx11::tuple_element<3,ArgList>::type arg4,
                  typename std_cxx11::tuple_element<4,ArgList>::type arg5,
                  typename std_cxx11::tuple_element<5,ArgList>::type arg6,
                  typename std_cxx11::tuple_element<6,ArgList>::type arg7)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<3,ArgList>::type>::act(arg4),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<4,ArgList>::type>::act(arg5),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<5,ArgList>::type>::act(arg6),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<6,ArgList>::type>::act(arg7)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 8 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 8>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3,
                  typename std_cxx11::tuple_element<3,ArgList>::type arg4,
                  typename std_cxx11::tuple_element<4,ArgList>::type arg5,
                  typename std_cxx11::tuple_element<5,ArgList>::type arg6,
                  typename std_cxx11::tuple_element<6,ArgList>::type arg7,
                  typename std_cxx11::tuple_element<7,ArgList>::type arg8)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<3,ArgList>::type>::act(arg4),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<4,ArgList>::type>::act(arg5),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<5,ArgList>::type>::act(arg6),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<6,ArgList>::type>::act(arg7),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<7,ArgList>::type>::act(arg8)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };

    /**
     * @internal
     * Encapsulator class for
     * functions with 9 arguments.
     */
    template <typename RT, typename ArgList>
    class fun_encapsulator<RT, ArgList, 9>
    {
    public:
      fun_encapsulator (typename internal::fun_ptr<RT,ArgList>::type *function)
        : function (*function)
      {}

      fun_encapsulator (const std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> &function)
        : function (function)
      {}

      inline
      Thread<RT>
      operator() (typename std_cxx11::tuple_element<0,ArgList>::type arg1,
                  typename std_cxx11::tuple_element<1,ArgList>::type arg2,
                  typename std_cxx11::tuple_element<2,ArgList>::type arg3,
                  typename std_cxx11::tuple_element<3,ArgList>::type arg4,
                  typename std_cxx11::tuple_element<4,ArgList>::type arg5,
                  typename std_cxx11::tuple_element<5,ArgList>::type arg6,
                  typename std_cxx11::tuple_element<6,ArgList>::type arg7,
                  typename std_cxx11::tuple_element<7,ArgList>::type arg8,
                  typename std_cxx11::tuple_element<8,ArgList>::type arg9)
      {
        return
          Thread<RT>
          (std_cxx11::bind (function,
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<0,ArgList>::type>::act(arg1),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<1,ArgList>::type>::act(arg2),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<2,ArgList>::type>::act(arg3),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<3,ArgList>::type>::act(arg4),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<4,ArgList>::type>::act(arg5),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<5,ArgList>::type>::act(arg6),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<6,ArgList>::type>::act(arg7),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<7,ArgList>::type>::act(arg8),
                            internal::maybe_make_ref<typename std_cxx11::tuple_element<8,ArgList>::type>::act(arg9)));
      }

    private:
      std_cxx11::function<typename internal::fun_ptr<RT,ArgList>::type> function;
    };
  }




// ----------- encapsulators for functions not taking any parameters

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with no arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<>,0>
  spawn (RT (*fun_ptr)()) DEAL_II_DEPRECATED;


  template <typename RT>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<>,0>
  spawn (RT (*fun_ptr)())
  {
    return fun_ptr;
  }


  /**
   * Overload of the non-const spawn
   * function for member functions with
   * no arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<>,0>
  spawn (C &c, RT (C::*fun_ptr)()) DEAL_II_DEPRECATED;


  template <typename RT, typename C>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<>,0>
  spawn (C &c, RT (C::*fun_ptr)())
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c)));
  }

  /**
   * Overload of the spawn function for
   * const member functions with no
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<>,0>
  spawn (const C &c, RT (C::*fun_ptr)() const) DEAL_II_DEPRECATED;


  template <typename RT, typename C>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<>,0>
  spawn (const C &c, RT (C::*fun_ptr)() const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c)));
  }




// ----------- encapsulators for unary functions

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 1 argument.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename Arg1>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1>,1>
  spawn (RT (*fun_ptr)(Arg1)) DEAL_II_DEPRECATED;


  template <typename RT, typename Arg1>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1>,1>
  spawn (RT (*fun_ptr)(Arg1))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 1 argument.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C, typename Arg1>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1>,1>
  spawn (C &c, RT (C::*fun_ptr)(Arg1)) DEAL_II_DEPRECATED;


  template <typename RT, typename C, typename Arg1>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1>,1>
  spawn (C &c, RT (C::*fun_ptr)(Arg1))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 1
   * argument.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C, typename Arg1>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1>,1>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C, typename Arg1>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1>,1>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1));
  }


// ----------- encapsulators for binary functions

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 2 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2>,2>
  spawn (RT (*fun_ptr)(Arg1,Arg2)) DEAL_II_DEPRECATED;


  template <typename RT, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2>,2>
  spawn (RT (*fun_ptr)(Arg1,Arg2))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 2 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2>,2>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2)) DEAL_II_DEPRECATED;


  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2>,2>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1,Arg2> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 2
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2>,2>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2>,2>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1,Arg2> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2));
  }


// ----------- encapsulators for ternary functions

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 3 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3>,3>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3>,3>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 3 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3>,3>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3>,3>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1,Arg2,Arg3> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 3
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3>,3>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3>,3>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1,Arg2,Arg3> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3));
  }



// ----------- encapsulators for functions with 4 arguments

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 4 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4>,4>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4>,4>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 4 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4>,4>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4>,4>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 4
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4>,4>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4>,4>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4));
  }


// ----------- encapsulators for functions with 5 arguments

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 5 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5>,5>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5>,5>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 5 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5>,5>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5>,5>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 5
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5>,5>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5>,5>
  spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5));
  }


// ----------- encapsulators for functions with 6 arguments

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 6 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,6>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,6>
  spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 6 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,6>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,6>
  spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 6
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,6>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,6>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6));
  }


// ----------- encapsulators for functions with 7 arguments

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 7 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6, Arg7>,7>
           spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6, Arg7>,7>
           spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 7 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6, Arg7>,7>
           spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6, Arg7>,7>
           spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6, std_cxx11::_7));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 7
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6, Arg7>,7>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6, Arg7>,7>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6, std_cxx11::_7));
  }


// ----------- encapsulators for functions with 8 arguments

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 8 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8>,8>
           spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8>,8>
           spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 8 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8>,8>
           spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                         Arg6,Arg7,Arg8)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8>,8>
           spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                         Arg6,Arg7,Arg8))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, Arg8> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6, std_cxx11::_7, std_cxx11::_8));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 8
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8>,8>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                               Arg6,Arg7,Arg8) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8>,8>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                               Arg6,Arg7,Arg8) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, Arg8> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6, std_cxx11::_7, std_cxx11::_8));
  }


// ----------- encapsulators for functions with 9 arguments

  /**
   * Overload of the spawn function for
   * non-member or static member
   * functions with 9 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8, Arg9>,9>
           spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8,Arg9)) DEAL_II_DEPRECATED;


  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8, Arg9>,9>
           spawn (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                Arg6,Arg7,Arg8,Arg9))
  {
    return fun_ptr;
  }



  /**
   * Overload of the non-const spawn
   * function for member functions with
   * 9 arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8, Arg9>,9>
           spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                         Arg6,Arg7,Arg8,Arg9)) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8, Arg9>,9>
           spawn (C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                         Arg6,Arg7,Arg8,Arg9))
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, Arg8, Arg9> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6, std_cxx11::_7, std_cxx11::_8, std_cxx11::_9));
  }

  /**
   * Overload of the spawn function for
   * const member functions with 9
   * arguments.
   *
   * @deprecated Use new_thread() instead.
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8, Arg9>,9>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                               Arg6,Arg7,Arg8,Arg9) const) DEAL_II_DEPRECATED;


  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  internal::fun_encapsulator<RT,
           std_cxx11::tuple<Arg1, Arg2, Arg3,
           Arg4, Arg5, Arg6,
           Arg7, Arg8, Arg9>,9>
           spawn (const C &c, RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                                               Arg6,Arg7,Arg8,Arg9) const)
  {
    return
      std_cxx11::function<typename internal::fun_ptr<RT,std_cxx11::tuple<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, Arg8, Arg9> >::type>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c), std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4, std_cxx11::_5, std_cxx11::_6, std_cxx11::_7, std_cxx11::_8, std_cxx11::_9));
  }



// ----------- thread starters for functions not taking any parameters

  /**
   * Overload of the new_thread function for
   * objects that can be converted to
   * std_cxx11::function<RT ()>, i.e. anything
   * that can be called like a function
   * object without arguments and returning
   * an object of type RT (or void).
   *
   * @ingroup threads
   */
  template <typename RT>
  inline
  Thread<RT>
  new_thread (const std_cxx11::function<RT ()> &function)
  {
    return Thread<RT>(function);
  }

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with no arguments.
   *
   * @ingroup threads
   */
  template <typename RT>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)())
  {
    return Thread<RT>(fun_ptr);
  }


  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * no arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(),
              typename identity<C>::type &c)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c)));
  }


#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with no
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)() const,
              const typename identity<C>::type &c)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c)));
  }
#endif



// ----------- thread starters for unary functions

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 1 argument.
   *
   * @ingroup threads
   */
  template <typename RT, typename Arg1>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1),
              typename identity<Arg1>::type arg1)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 1 argument.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 1
   * argument.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1)));
  }
#endif

// ----------- thread starters for binary functions

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 2 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename Arg1, typename Arg2>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 2 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 2
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2)));
  }
#endif

// ----------- thread starters for ternary functions

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 3 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 3 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 3
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3)));
  }
#endif


// ----------- thread starters for functions with 4 arguments

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 4 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 4 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 4
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4)));
  }
#endif

// ----------- thread starters for functions with 5 arguments

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 5 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 5 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 5
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5)));
  }
#endif

// ----------- thread starters for functions with 6 arguments

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 6 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 6 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 6
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6)));
  }
#endif

// ----------- thread starters for functions with 7 arguments

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 7 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 7 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 7
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7)));
  }
#endif

// ----------- thread starters for functions with 8 arguments

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 8 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                            Arg6,Arg7,Arg8),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7,
              typename identity<Arg8>::type arg8)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 8 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                               Arg6,Arg7,Arg8),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7,
              typename identity<Arg8>::type arg8)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 8
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                               Arg6,Arg7,Arg8) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7,
              typename identity<Arg8>::type arg8)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8)));
  }
#endif

// ----------- thread starters for functions with 9 arguments

  /**
   * Overload of the new_thread function for
   * non-member or static member
   * functions with 9 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  Thread<RT>
  new_thread (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                            Arg6,Arg7,Arg8,Arg9),
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7,
              typename identity<Arg8>::type arg8,
              typename identity<Arg9>::type arg9)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8),
                       internal::maybe_make_ref<Arg9>::act(arg9)));
  }



  /**
   * Overload of the non-const new_thread
   * function for member functions with
   * 9 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                               Arg6,Arg7,Arg8,Arg9),
              typename identity<C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7,
              typename identity<Arg8>::type arg8,
              typename identity<Arg9>::type arg9)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8),
                       internal::maybe_make_ref<Arg9>::act(arg9)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_thread function for
   * const member functions with 9
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  Thread<RT>
  new_thread (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                               Arg6,Arg7,Arg8,Arg9) const,
              typename identity<const C>::type &c,
              typename identity<Arg1>::type arg1,
              typename identity<Arg2>::type arg2,
              typename identity<Arg3>::type arg3,
              typename identity<Arg4>::type arg4,
              typename identity<Arg5>::type arg5,
              typename identity<Arg6>::type arg6,
              typename identity<Arg7>::type arg7,
              typename identity<Arg8>::type arg8,
              typename identity<Arg9>::type arg9)
  {
    return
      Thread<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8),
                       internal::maybe_make_ref<Arg9>::act(arg9)));
  }
#endif

// ------------------------ ThreadGroup -------------------------------------

  /**
   * A container for thread
   * objects. Allows to add new
   * thread objects and wait for them
   * all together. The thread objects
   * need to have the same return
   * value for the called function.
   *
   * @author Wolfgang Bangerth, 2003
   * @ingroup threads
   */
  template <typename RT = void>
  class ThreadGroup
  {
  public:
    /**
     * Add another thread object to
     * the collection.
     */
    ThreadGroup &operator += (const Thread<RT> &t)
    {
      threads.push_back (t);
      return *this;
    }

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
    void join_all () const
    {
      for (typename std::list<Thread<RT> >::const_iterator
           t=threads.begin(); t!=threads.end(); ++t)
        t->join ();
    }

  private:
    /**
     * List of thread objects.
     */
    std::list<Thread<RT> > threads;
  };


  template <typename> class Task;


  namespace internal
  {
#ifdef DEAL_II_WITH_THREADS

    template <typename> struct TaskDescriptor;

    /**
     * The task class for TBB that is
     * used by the TaskDescriptor
     * class.
     */
    template <typename RT>
    struct TaskEntryPoint : public tbb::task
    {
      TaskEntryPoint (TaskDescriptor<RT> &task_descriptor)
        :
        task_descriptor (task_descriptor)
      {}

      virtual tbb::task *execute ()
      {
        // call the function object
        // and put the return value
        // into the proper place
        try
          {
            call (task_descriptor.function, task_descriptor.ret_val);
          }
        catch (const std::exception &exc)
          {
            internal::handle_std_exception (exc);
          }
        catch (...)
          {
            internal::handle_unknown_exception ();
          }
        return 0;
      }

      /**
       * A reference to the descriptor
       * object of this task.
       */
      TaskDescriptor<RT> &task_descriptor;
    };

    /**
     * @internal
     * Base class describing a
     * task. This is the basic
     * class abstracting the
     * Threading Building Blocks
     * implementation of tasks.
     * It provides a mechanism
     * to start a new task, as well
     * as for joining it.
     *
     * Internally, the way things are
     * implemented is that all Task<>
     * objects keep a shared pointer
     * to the task descriptor. When
     * the last Task<> object goes
     * out of scope, the destructor
     * of the descriptor is
     * called. Since tasks can not be
     * abandoned, the destructor
     * makes sure that the task is
     * finished before it can
     * continue to destroy the
     * object.
     *
     * Note that unlike threads,
     * tasks are not always started
     * right away, and so the
     * starting thread can't rely on
     * the fact that the started task
     * can copy things off the
     * spawning thread's stack
     * frame. As a consequence, the
     * task description needs to
     * include a way to store the
     * function and its arguments
     * that shall be run on the task.
     *
     * @author Wolfgang Bangerth, 2009
     */
    template <typename RT>
    struct TaskDescriptor
    {
    private:
      /**
       * The function and its arguments
       * that are to be run on the task.
       */
      std_cxx11::function<RT ()> function;

      /**
       * Variable holding the data the TBB
       * needs to work with a task. Set by
       * the queue_up_task() function. Note
       * that the object behind this
       * pointer will be deleted upon
       * termination of the task, so we do
       * not have to do so ourselves. In
       * particular, if all objects with
       * pointers to this task_description
       * object go out of scope then no
       * action is needed on our behalf.
       */
      tbb::task *task;

      /**
       * A place where the task will
       * deposit its return value.
       */
      return_value<RT> ret_val;

      /**
       * A flag indicating whether the task
       * has terminated.
       */
      bool task_is_done;

    public:

      /**
       * Constructor. Take the function to
       * be run on this task as argument.
       */
      TaskDescriptor (const std_cxx11::function<RT ()> &function);

      /**
       * Default
       * constructor. Throws an
       * exception since we want to
       * queue a task immediately
       * upon construction of these
       * objects to make sure that
       * each TaskDescriptor object
       * corresponds to exactly one
       * task.
       */
      TaskDescriptor ();

      /**
       * Copy constructor. Throws
       * an exception since we want
       * to make sure that each
       * TaskDescriptor object
       * corresponds to exactly one
       * task.
       */
      TaskDescriptor (const TaskDescriptor &);

      /**
       * Destructor.
       */
      ~TaskDescriptor ();

      /**
       * Queue up the task to the
       * scheduler. We need to do
       * this in a separate
       * function since the new
       * tasks needs to access
       * objects from the current
       * object and that can only
       * reliably happen if the
       * current object is
       * completely constructed
       * already.
       */
      void queue_task ();

      /**
       * Join a task, i.e. wait
       * for it to finish. This
       * function can safely be
       * called from different
       * threads at the same time,
       * and can also be called
       * more than once.
       */
      void join ();


      template <typename> friend struct TaskEntryPoint;
      friend class dealii::Threads::Task<RT>;
    };



    template <typename RT>
    inline
    TaskDescriptor<RT>::TaskDescriptor (const std_cxx11::function<RT ()> &function)
      :
      function (function),
      task_is_done (false)
    {}


    template <typename RT>
    inline
    void
    TaskDescriptor<RT>::queue_task ()
    {
      // use the pattern described in
      // the TBB book on pages
      // 230/231 ("Start a large task
      // in parallel with the main
      // program")
      task = new (tbb::task::allocate_root()) tbb::empty_task;
      task->set_ref_count (2);

      tbb::task *worker = new (task->allocate_child()) TaskEntryPoint<RT>(*this);
      task->spawn (*worker);
    }



    template <typename RT>
    TaskDescriptor<RT>::TaskDescriptor ()
    {
      Assert (false, ExcInternalError());
    }



    template <typename RT>
    TaskDescriptor<RT>::TaskDescriptor (const TaskDescriptor &)
    {
      Assert (false, ExcInternalError());
    }



    template <typename RT>
    inline
    TaskDescriptor<RT>::~TaskDescriptor ()
    {
      // wait for the task to
      // complete for sure
      join ();

      // now destroy the empty task
      // structure. the book recommends to
      // spawn it as well and let the
      // scheduler destroy the object when
      // done, but this has the disadvantage
      // that the scheduler may not get to
      // actually finishing the task before
      // it goes out of scope (at the end of
      // the program, or if a thread is done
      // on which it was run) and then we
      // would get a hard-to-decipher warning
      // about unfinished tasks when the
      // scheduler "goes out of the
      // arena". rather, let's explicitly
      // destroy the empty task
      // object. before that, make sure that
      // the task has been shut down,
      // expressed by a zero reference count
      Assert (task != 0, ExcInternalError());
      Assert (task->ref_count()==0, ExcInternalError());
      task->destroy (*task);
    }


    template <typename RT>
    inline
    void
    TaskDescriptor<RT>::join ()
    {
      // if the task is already done, just
      // return. this makes sure we call
      // tbb::Task::wait_for_all() exactly
      // once, as required by TBB. we could
      // also get the reference count of task
      // for doing this, but that is usually
      // slower. note that this does not work
      // when the thread calling this
      // function is not the same as the one
      // that initialized the task.
      //
      // TODO: can we assert that no other
      // thread tries to end the task?
      if (task_is_done == true)
        return;

      // let TBB wait for the task to
      // complete.
      task_is_done = true;
      task->wait_for_all();
    }



#else        // no threading enabled

    /**
     * A way to describe tasks. Since
     * we are in non-MT mode at this
     * place, things are a lot
     * simpler than in MT mode.
     */
    template <typename RT>
    struct TaskDescriptor
    {
      /**
       * A place where the task will
       * deposit its return value.
       */
      return_value<RT> ret_val;

      /**
       * Constructor. Call the
       * given function and emplace
       * the return value into the
       * slot reserved for this
       * purpose.
       */
      TaskDescriptor (const std_cxx11::function<RT ()> &function)
      {
        call (function, ret_val);
      }

      /**
       * Wait for the task to
       * return. Since we are in
       * non-MT mode here, there is
       * nothing to do.
       */
      static void join () {}

      /**
       * Run the task. Since we are
       * here in non-MT mode, there
       * is nothing to do that the
       * constructor hasn't already
       * done.
       */
      static void queue_task () {}
    };

#endif

  }



  /**
   * Describes one task object based on the
   * Threading Building Blocks' Task. Note
   * that the call to join() must be executed
   * on the same thread as the call to the
   * constructor. Otherwise, there might be a
   * deadlock. In other words, a Task
   * object should never passed on to another
   * task for calling the join() method.
   *
   * @author Wolfgang Bangerth, 2009
   * @ingroup threads
   */
  template <typename RT = void>
  class Task
  {
  public:
    /**
     * Construct a task object given a function object to execute on
     * the task, and then schedule this function for execution.
     *
     * @post Using this constructor automatically makes the task
     * object joinable().
     */
    Task (const std_cxx11::function<RT ()> &function_object)
    {
      // create a task descriptor and tell it
      // to queue itself up with the scheduling
      // system
      task_descriptor =
        std_cxx11::shared_ptr<internal::TaskDescriptor<RT> >
        (new internal::TaskDescriptor<RT>(function_object));
      task_descriptor->queue_task ();
    }


    /**
     * Copy constructor.
     *
     * @post Using this constructor automatically makes the task
     * object joinable().
     */
    Task (const Task<RT> &t)
      :
      task_descriptor (t.task_descriptor)
    {}


    /**
     * Default constructor. You can't do much with a task object
     * constructed this way, except for assigning it a task object
     * that holds data created by the Threads::new_task() functions.
     *
     * @post Using this constructor leaves the object in an unjoinable
     *   state, i.e., joinable() will return false.
     */
    Task () {}

    /**
     * Join the task represented by this object, i.e. wait for it to
     * finish.
     *
     * @pre You can't call this function if you have used the default
     *   constructor of this class and have not assigned a task object
     *   to it. In other words, the function joinable() must return
     *   true.
     */
    void join () const
    {
      AssertThrow (joinable(), ExcNoTask());
      task_descriptor->join ();
    }

    /**
     * Return whether the current object can be joined. You can join a
     * task object once a task (typically created with
     * Threads::new_task()) has actually been assigned to it. On the
     * other hand, the function returns false if it has been default
     * constructed.
     */
    bool joinable () const
    {
      return (task_descriptor !=
              std_cxx11::shared_ptr<internal::TaskDescriptor<RT> >());
    }


    /**
     * Get the return value of the function of the task. Since this is
     * only available once the task finishes, this implicitly also
     * calls <tt>join()</tt>.
     *
     * @pre You can't call this function if you have used the default
     *   constructor of this class and have not assigned a task object
     *   to it. In other words, the function joinable() must return
     *   true.
     */
    RT return_value ()
    {
      join ();
      return task_descriptor->ret_val.get();
    }


    /**
     * Check for equality of task
     * objects. Since objects of
     * this class store an implicit
     * pointer to an object that
     * exists exactly once for each
     * task, the check is simply
     * to compare these pointers.
     */
    bool operator == (const Task &t)
    {
      AssertThrow (joinable(), ExcNoTask());
      return task_descriptor == t.task_descriptor;
    }

    /** @addtogroup Exceptions
     * @{ */

    /**
     * Exception
     */
    DeclException0 (ExcNoTask);
    //@}
  private:
    /**
     * Shared pointer to the object
     * representing the task. Boost's
     * shared pointer implementation will
     * make sure that that object lives as
     * long as there is at least one
     * subscriber to it.
     */
    std_cxx11::shared_ptr<internal::TaskDescriptor<RT> > task_descriptor;
  };



  /**
   * Overload of the new_task function for
   * objects that can be converted to
   * std_cxx11::function<RT ()>, i.e. anything
   * that can be called like a function
   * object without arguments and returning
   * an object of type RT (or void).
   *
   * @ingroup threads
   */
  template <typename RT>
  inline
  Task<RT>
  new_task (const std_cxx11::function<RT ()> &function)
  {
    return Task<RT>(function);
  }

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with no arguments.
   *
   * @ingroup threads
   */
  template <typename RT>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)())
  {
    return new_task<RT>(std_cxx11::function<RT ()>(fun_ptr));
  }


  /**
   * Overload of the non-const new_task
   * function for member functions with
   * no arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(),
            typename identity<C>::type &c)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with no
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)() const,
            const typename identity<C>::type &c)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c)));
  }
#endif



// ----------- thread starters for unary functions

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 1 argument.
   *
   * @ingroup threads
   */
  template <typename RT, typename Arg1>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1),
            typename identity<Arg1>::type arg1)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 1 argument.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 1
   * argument.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1)));
  }
#endif

// ----------- thread starters for binary functions

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 2 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename Arg1, typename Arg2>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 2 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 2
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename Arg1, typename Arg2>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2)));
  }
#endif

// ----------- thread starters for ternary functions

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 3 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 3 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 3
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3)));
  }
#endif


// ----------- thread starters for functions with 4 arguments

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 4 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 4 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 4
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4)));
  }
#endif

// ----------- thread starters for functions with 5 arguments

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 5 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 5 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 5
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5)));
  }
#endif

// ----------- thread starters for functions with 6 arguments

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 6 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 6 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 6
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6)));
  }
#endif

// ----------- thread starters for functions with 7 arguments

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 7 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 7 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 7
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7)));
  }
#endif

// ----------- thread starters for functions with 8 arguments

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 8 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                          Arg6,Arg7,Arg8),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7,
            typename identity<Arg8>::type arg8)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 8 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                             Arg6,Arg7,Arg8),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7,
            typename identity<Arg8>::type arg8)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 8
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                             Arg6,Arg7,Arg8) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7,
            typename identity<Arg8>::type arg8)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8)));
  }
#endif

// ----------- thread starters for functions with 9 arguments

  /**
   * Overload of the new_task function for
   * non-member or static member
   * functions with 9 arguments.
   *
   * @ingroup threads
   */
  template <typename RT,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  Task<RT>
  new_task (RT (*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                          Arg6,Arg7,Arg8,Arg9),
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7,
            typename identity<Arg8>::type arg8,
            typename identity<Arg9>::type arg9)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr,
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8),
                       internal::maybe_make_ref<Arg9>::act(arg9)));
  }



  /**
   * Overload of the non-const new_task
   * function for member functions with
   * 9 arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                             Arg6,Arg7,Arg8,Arg9),
            typename identity<C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7,
            typename identity<Arg8>::type arg8,
            typename identity<Arg9>::type arg9)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::ref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8),
                       internal::maybe_make_ref<Arg9>::act(arg9)));
  }

#ifndef DEAL_II_CONST_MEMBER_DEDUCTION_BUG
  /**
   * Overload of the new_task function for
   * const member functions with 9
   * arguments.
   *
   * @ingroup threads
   */
  template <typename RT, typename C,
            typename Arg1, typename Arg2, typename Arg3,
            typename Arg4, typename Arg5, typename Arg6,
            typename Arg7, typename Arg8, typename Arg9>
  inline
  Task<RT>
  new_task (RT (C::*fun_ptr)(Arg1,Arg2,Arg3,Arg4,Arg5,
                             Arg6,Arg7,Arg8,Arg9) const,
            typename identity<const C>::type &c,
            typename identity<Arg1>::type arg1,
            typename identity<Arg2>::type arg2,
            typename identity<Arg3>::type arg3,
            typename identity<Arg4>::type arg4,
            typename identity<Arg5>::type arg5,
            typename identity<Arg6>::type arg6,
            typename identity<Arg7>::type arg7,
            typename identity<Arg8>::type arg8,
            typename identity<Arg9>::type arg9)
  {
    return
      new_task<RT>
      (std_cxx11::bind(fun_ptr, std_cxx11::cref(c),
                       internal::maybe_make_ref<Arg1>::act(arg1),
                       internal::maybe_make_ref<Arg2>::act(arg2),
                       internal::maybe_make_ref<Arg3>::act(arg3),
                       internal::maybe_make_ref<Arg4>::act(arg4),
                       internal::maybe_make_ref<Arg5>::act(arg5),
                       internal::maybe_make_ref<Arg6>::act(arg6),
                       internal::maybe_make_ref<Arg7>::act(arg7),
                       internal::maybe_make_ref<Arg8>::act(arg8),
                       internal::maybe_make_ref<Arg9>::act(arg9)));
  }
#endif


// ------------------------ TaskGroup -------------------------------------

  /**
   * A container for task
   * objects. Allows to add new
   * task objects and wait for them
   * all together. The task objects
   * need to have the same return
   * value for the called function.
   *
   * Note that the call to join_all() must be
   * executed on the same thread as the calls
   * that add subtasks. Otherwise, there
   * might be a deadlock. In other words, a
   * Task object should never passed on to
   * another task for calling the join()
   * method.
   *
   * @author Wolfgang Bangerth, 2003
   * @ingroup tasks
   */
  template <typename RT = void>
  class TaskGroup
  {
  public:
    /**
     * Add another task object to
     * the collection.
     */
    TaskGroup &operator += (const Task<RT> &t)
    {
      tasks.push_back (t);
      return *this;
    }

    /**
     * Wait for all tasks in the
     * collection to finish. It is
     * not a problem if some of
     * them have already been
     * waited for, i.e. you may
     * call this function more than
     * once, and you can also add
     * new task objects between
     * subsequent calls to this
     * function if you want.
     */
    void join_all () const
    {
      for (typename std::list<Task<RT> >::const_iterator
           t=tasks.begin(); t!=tasks.end(); ++t)
        t->join ();
    }

  private:
    /**
     * List of task objects.
     */
    std::list<Task<RT> > tasks;
  };

}   // end of implementation of namespace Threads

/**
 * @}
 */


//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef __deal2__thread_management_h
#endif
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__thread_local_storage_h
#define __deal2__thread_local_storage_h


#include <deal.II/base/config.h>

#if DEAL_II_USE_MT == 1
#  include <tbb/enumerable_thread_specific.h>
#endif



DEAL_II_NAMESPACE_OPEN

/*!@addtogroup threads */
/*@{*/


namespace Threads
{
  template <typename T>
  /**
   * @brief A class that provides a separate storage location on each thread that accesses the object.
   * 
   * This class offers ways so that every thread that accesses it has its
   * own copy of an object of type T. In essence, accessing this object
   * can never result in race conditions in multithreaded programs since
   * no other thread than the current one can ever access it.
   * 
   * The class builds on the Threading Building Blocks's tbb::enumerable_thread_specific
   * class but wraps it in such a way that this class can also be used when deal.II
   * is configured not to use threads at all -- in that case, this class simply
   * stores a single copy of an object of type T.
   * 
   * <h3>Construction and destruction</h3>
   * 
   * Objects of this class can either be default constructed or by providing an
   * "exemplar", i.e. an object of type T so that everytime we need to create
   * a T on a thread that doesn't already have such an object, it is copied from
   * the exemplar.
   * 
   * Upon destruction of objects of this class, all T objects that correspond
   * to threads that have accessed this object are destroyed. Note that this may be
   * before the time when a thread is terminated.
   * 
   * <h3>Access</h3>
   * 
   * The T object stored by this object can be accessed using the get() function. It
   * provides a reference to a unique object when accessed from different threads.
   * Objects of type T are created lazily, i.e. they are only created whenever a
   * thread actually calls get().
   **/
  class ThreadLocalStorage 
  {
  public:
    /**
     * Default constructor. Initialize each thread local object
     * using its default constructor.
     **/
    ThreadLocalStorage ();
    
    /**
     * A kind of copy constructor. Initialize each thread local object
     * by copying the given object.
     **/
    ThreadLocalStorage (const T &t);

    /**
     * Copy constructor. Initialize each thread local object
     * with the corresponding object of the given object.
     **/
    ThreadLocalStorage (const ThreadLocalStorage<T> &t);

    /**
     * Return a reference to the data stored by this object for the current
     * thread this function is called on.
     * 
     * Note that there is no member function get() that is const and
     * returns a const reference as one would expect. The reason is that
     * if such a member function were called on a thread for which no
     * thread-local object has been created yet, then one has to create
     * such an object first which would certainly be a non-constant
     * operation. If you need to call the get() function for a member
     * variable of a class from a const member function, then you
     * need to declare the member variable <code>mutable</code> to
     * allow such access.
     */
    T & get ();
    
  private:
#if DEAL_II_USE_MT == 1
    /**
     * The data element we store. If we support threads, then this
     * object will be of a type that provides a separate object
     * for each thread. Otherwise, it is simply a single object 
     * of type T.
     **/
    tbb::enumerable_thread_specific<T> data;
#else
    T data;
#endif    
  };
  
// ----------------- inline and template functions ----------------------------
  
  template <typename T>
  inline
  ThreadLocalStorage<T>::ThreadLocalStorage()
  {}


  template <typename T>
  inline
  ThreadLocalStorage<T>::ThreadLocalStorage(const T &t)
  :
  data (t)
  {}


  template <typename T>
  inline
  ThreadLocalStorage<T>::ThreadLocalStorage(const ThreadLocalStorage<T> &t)
  :
  data (t)
  {}


  template <typename T>
  inline
  T &
  ThreadLocalStorage<T>::get ()
  {
#if DEAL_II_USE_MT == 1
    return data.local();
#else
    return data;
#endif        
  }

  
}   // end of implementation of namespace Threads

/**
 * @}
 */


//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef __deal2__thread_local_storage_h
#endif
//---------------------------------------------------------------------------

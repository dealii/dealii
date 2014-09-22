// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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

#ifndef __deal2__thread_local_storage_h
#define __deal2__thread_local_storage_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_THREADS
#  include <tbb/enumerable_thread_specific.h>
#endif



DEAL_II_NAMESPACE_OPEN

/*!@addtogroup threads */
/*@{*/


namespace Threads
{
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
   * "exemplar", i.e. an object of type T so that every time we need to create
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
   */
  template <typename T>
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
    explicit ThreadLocalStorage (const T &t);

    /**
     * Copy constructor. Initialize each thread local object
     * with the corresponding object of the given object.
     */
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
    T &get ();

    /**
     * Same as above, except that @p exists is set to true if an element
     * was already present for the current thread; false otherwise.
     */
    T &get (bool &exists);

    /**
     * Conversion operator that simply converts the thread-local object
     * to the data type that it stores. This function is equivalent to
     * calling the get() member function; it's purpose is to make the
     * TLS object look more like the object it is storing.
     */
    operator T &();

    /**
     * Copy the given argument into the storage space used to represent
     * the current thread. Calling this function as <code>tls_data = object</code>
     * is equivalent to calling <code>tls_data.get() = object</code>. The
     * intent of this operator is to make the ThreadLocalStorage object
     * look more like the object it represents on the current thread.
     *
     * @param t The object to be copied into the storage space used
     * for the current thread.
     *
     * @return The current object, after the changes have been made
     **/
    ThreadLocalStorage<T> &operator = (const T &t);

    /**
     * Remove the thread-local objects stored for all threads that have
     * created one with this object (i.e., that have called get()
     * at least once on this thread. This includes the current thread. If you
     * call get() subsequently on this or any other thread, new objects will
     * again be created.
     *
     * If deal.II has been configured to not use multithreading, then this function
     * does not do anything at all. Note that this of course has different semantics
     * as in the multithreading context the objects are deleted and created again
     * (possible by copying from a sample object, if the appropriate constructor
     * of this class was called), whereas in the multithreaded context the object
     * is simply not touched at all. At the same time, the purpose of this function
     * is to release memory other threads may have allocated for their own thread
     * local objects after which every use of this object will require some kind
     * of initialization. This is necessary both in the multithreaded or
     * non-multithreaded case.
     */
    void clear ();

    /**
     * Returns a reference to the internal Threading Building Blocks
     * implementation. This function is really only useful if deal.II
     * has been configured with multithreading and has no useful
     * purpose otherwise.
     */
#ifdef DEAL_II_WITH_THREADS
    tbb::enumerable_thread_specific<T> &
#else
    T &
#endif
    get_implementation();

  private:
#ifdef DEAL_II_WITH_THREADS
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
#ifdef DEAL_II_WITH_THREADS
    return data.local();
#else
    return data;
#endif
  }


  template <typename T>
  inline
  T &
  ThreadLocalStorage<T>::get (bool &exists)
  {
#ifdef DEAL_II_WITH_THREADS
    return data.local(exists);
#else
    exists = true;
    return data;
#endif
  }


  template <typename T>
  inline
  ThreadLocalStorage<T>::operator T &()
  {
    return get();
  }


  template <typename T>
  inline
  ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator = (const T &t)
  {
    get() = t;
    return *this;
  }


  template <typename T>
  inline
#ifdef DEAL_II_WITH_THREADS
  tbb::enumerable_thread_specific<T> &
#else
  T &
#endif
  ThreadLocalStorage<T>::get_implementation()
  {
    return data;
  }



  template <typename T>
  inline
  void
  ThreadLocalStorage<T>::clear ()
  {
#ifdef DEAL_II_WITH_THREADS
    data.clear ();
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

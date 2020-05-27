// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_thread_local_storage_h
#  define dealii_thread_local_storage_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>

#  include <list>
#  include <map>
#  include <shared_mutex>
#  include <thread>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup threads */
/*@{*/

#  ifndef DOXYGEN
class LogStream;
#  endif

namespace Threads
{
#  ifndef DOXYGEN
  namespace internal
  {
    /*
     * Workaround: The standard unfortunately has an unfortunate design
     * "flaw" in the std::is_copy_constructible type trait
     * when it comes to STL containers and containing non-copyable objects
     * T. The type trait is true even though any attempted invocation leads
     * to a compilation error. Work around this issue by unpacking some
     * commonly used containers:
     */
    template <typename T>
    struct unpack_container
    {
      using type = T;
    };

    template <typename T, typename A>
    struct unpack_container<std::vector<T, A>>
    {
      using type = T;
    };

    template <typename T, typename A>
    struct unpack_container<std::list<T, A>>
    {
      using type = T;
    };
  } // namespace internal
#  endif

  /**
   * @brief A class that provides a separate storage location on each thread
   * that accesses the object.
   *
   * This class offers ways so that every thread that accesses it has its own
   * copy of an object of type T. In essence, accessing this object can never
   * result in race conditions in multithreaded programs since no other thread
   * than the current one can ever access it.
   *
   * <h3>Construction and destruction</h3>
   *
   * Objects of this class can either be default constructed or by providing
   * an "exemplar", i.e. an object of type T so that every time we need to
   * create a T on a thread that doesn't already have such an object, it is
   * copied from the exemplar.
   *
   * Upon destruction of objects of this class, all T objects that correspond
   * to threads that have accessed this object are destroyed. Note that this
   * may be before the time when a thread is terminated.
   *
   * <h3>Access</h3>
   *
   * The T object stored by this object can be accessed using the get()
   * function. It provides a reference to a unique object when accessed from
   * different threads. Objects of type T are created lazily, i.e. they are
   * only created whenever a thread actually calls get().
   */
  template <typename T>
  class ThreadLocalStorage
  {
    static_assert(
      std::is_copy_constructible<
        typename internal::unpack_container<T>::type>::value ||
        std::is_default_constructible<T>::value,
      "The stored type must be either copyable, or default constructible");

  public:
    /**
     * Default constructor. Initialize each thread local object using its
     * default constructor.
     */
    ThreadLocalStorage() = default;

    /**
     * A kind of copy constructor. Initializes an internal exemplar by the
     * given object. The exemplar is in turn used to initialize each thread
     * local object instead of invoking the default constructor.
     */
    explicit ThreadLocalStorage(const T &t);

    /**
     * A kind of move constructor. Moves the given object into an internal
     * exemplar. The exemplar is in turn used to initialize each thread
     * local object instead of invoking the default constructor.
     */
    explicit ThreadLocalStorage(T &&t);

    /**
     * Copy constructor. Initialize each thread local object with the
     * corresponding object of the given object.
     */
    ThreadLocalStorage(const ThreadLocalStorage<T> &t) = default;

    /**
     * Move constructor. Copies the internal state over from the given
     * object.
     */
    ThreadLocalStorage(ThreadLocalStorage<T> &&t) noexcept = default;

    /**
     * Return a reference to the data stored by this object for the current
     * thread this function is called on.
     *
     * Note that there is no member function get() that is const and returns a
     * const reference as one would expect. The reason is that if such a
     * member function were called on a thread for which no thread-local
     * object has been created yet, then one has to create such an object
     * first which would certainly be a non-constant operation. If you need to
     * call the get() function for a member variable of a class from a const
     * member function, then you need to declare the member variable
     * <code>mutable</code> to allow such access.
     */
    T &
    get();

    /**
     * Same as above, except that @p exists is set to true if an element was
     * already present for the current thread; false otherwise.
     */
    T &
    get(bool &exists);

    /**
     * Conversion operator that simply converts the thread-local object to the
     * data type that it stores. This function is equivalent to calling the
     * get() member function; it's purpose is to make the TLS object look more
     * like the object it is storing.
     */
    operator T &();

    /**
     * Copy the given argument into the storage space used to represent the
     * current thread. Calling this function as <code>tls_data = object</code>
     * is equivalent to calling <code>tls_data.get() = object</code>. The
     * intent of this operator is to make the ThreadLocalStorage object look
     * more like the object it represents on the current thread.
     *
     * @param t The object to be copied into the storage space used for the
     * current thread.
     *
     * @return The current object, after the changes have been made
     */
    ThreadLocalStorage<T> &
    operator=(const T &t);

    /**
     * Move the given argument into the storage space used to represent the
     * current thread. Calling this function as <code>tls_data =
     * object</code> is equivalent to calling <code>tls_data.get() =
     * object</code>. The intent of this operator is to make the
     * ThreadLocalStorage object look more like the object it represents on
     * the current thread. Move assignment operator.
     *
     * @param t The object to be copied into the storage space used for the
     * current thread.
     *
     * @return The current object, after the changes have been made
     */
    ThreadLocalStorage<T> &
    operator=(T &&t);

    /**
     * Remove the thread-local objects stored for all threads that have
     * created one with this object (i.e., that have called get() at least
     * once on this thread. This includes the current thread. If you call
     * get() subsequently on this or any other thread, new objects will again
     * be created.
     *
     * If deal.II has been configured to not use multithreading, then this
     * function does not do anything at all. Note that this of course has
     * different semantics as in the multithreading context the objects are
     * deleted and created again (possible by copying from a sample object, if
     * the appropriate constructor of this class was called), whereas in the
     * multithreaded context the object is simply not touched at all. At the
     * same time, the purpose of this function is to release memory other
     * threads may have allocated for their own thread local objects after
     * which every use of this object will require some kind of
     * initialization. This is necessary both in the multithreaded or non-
     * multithreaded case.
     */
    void
    clear();

  private:
    /**
     * The data element we store.
     */
    std::map<std::thread::id, T> data;

    /**
     * A mutex to guard insertion into the data object.
     *
     * We use a std::shared_timed_mutex (or std::shared_mutex if available)
     * here to be able to use std::unique_lock and std::shared_lock for a
     * readers-writer lock
     * (https://en.wikipedia.org/wiki/Readers%E2%80%93writer_lock).
     */
#  ifdef DEAL_II_HAVE_CXX17
    std::shared_mutex insertion_mutex;
#  else
    std::shared_timed_mutex insertion_mutex;
#  endif

    /**
     * An exemplar for creating a new (thread specific) copy.
     */
    std::shared_ptr<const T> exemplar;

    friend class dealii::LogStream;
  };
} // namespace Threads
/**
 * @}
 */

#  ifndef DOXYGEN
namespace Threads
{
  // ----------------- inline and template functions --------------------------

  template <typename T>
  inline ThreadLocalStorage<T>::ThreadLocalStorage(const T &t)
    : exemplar(std::make_shared<const T>(t))
  {}


  template <typename T>
  inline ThreadLocalStorage<T>::ThreadLocalStorage(T &&t)
    : exemplar(std::make_shared<T>(std::forward<T>(t)))
  {}


#    ifndef DOXYGEN
  namespace internal
  {
    /*
     * We have to make sure not to call "data.emplace(id, *exemplar)" if
     * the corresponding element is not copy constructible. We use some
     * SFINAE magic to work around the fact that C++14 does not have
     * "if constexpr".
     */
    template <typename T>
    typename std::enable_if_t<
      std::is_copy_constructible<typename unpack_container<T>::type>::value,
      T &>
    construct_element(std::map<std::thread::id, T> &  data,
                      const std::thread::id &         id,
                      const std::shared_ptr<const T> &exemplar)
    {
      if (exemplar)
        {
          const auto it = data.emplace(id, *exemplar).first;
          return it->second;
        }
      return data[id];
    }

    template <typename T>
    typename std::enable_if_t<
      !std::is_copy_constructible<typename unpack_container<T>::type>::value,
      T &>
    construct_element(std::map<std::thread::id, T> &data,
                      const std::thread::id &       id,
                      const std::shared_ptr<const T> &)
    {
      return data[id];
    }
  } // namespace internal
#    endif


  template <typename T>
  inline T &
  ThreadLocalStorage<T>::get(bool &exists)
  {
    const std::thread::id my_id = std::this_thread::get_id();

    // Note that std::map<..>::emplace guarantees that no iterators or
    // references to stored objects are invalidated. We thus only have to
    // ensure that we do not perform a lookup while writing, and that we
    // do not write concurrently. This is precisely the "reader-writer
    // lock" paradigm supported by C++14 by means of the std::shared_lock
    // and the std::unique_lock.

    {
      // Take a shared ("reader") lock for lookup and record the fact
      // whether we could find an entry in the boolean exists.
      std::shared_lock<decltype(insertion_mutex)> lock(insertion_mutex);

      const auto it = data.find(my_id);
      if (it != data.end())
        {
          exists = true;
          return it->second;
        }
      else
        {
          exists = false;
        }
    }

    {
      // Take a unique ("writer") lock for manipulating the std::map. This
      // lock ensures that no other threat does a lookup at the same time.
      std::unique_lock<decltype(insertion_mutex)> lock(insertion_mutex);

      return internal::construct_element(data, my_id, exemplar);
    }
  }


  template <typename T>
  inline T &
  ThreadLocalStorage<T>::get()
  {
    bool exists;
    return get(exists);
  }


  template <typename T>
  inline ThreadLocalStorage<T>::operator T &()
  {
    return get();
  }


  template <typename T>
  inline ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator=(const T &t)
  {
    get() = t;
    return *this;
  }


  template <typename T>
  inline ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator=(T &&t)
  {
    get() = std::forward<T>(t);
    return *this;
  }


  template <typename T>
  inline void
  ThreadLocalStorage<T>::clear()
  {
    std::unique_lock<decltype(insertion_mutex)> lock(insertion_mutex);
    data.clear();
  }
} // namespace Threads

#  endif // DOXYGEN

//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_thread_local_storage_h
#endif
//---------------------------------------------------------------------------

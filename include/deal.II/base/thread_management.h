// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2022 by the deal.II authors
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

#ifndef dealii_thread_management_h
#  define dealii_thread_management_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/multithread_info.h>
#  include <deal.II/base/std_cxx17/tuple.h>
#  include <deal.II/base/template_constraints.h>

#  include <atomic>
#  include <functional>
#  include <future>
#  include <list>
#  include <memory>
#  include <mutex>
#  include <thread>
#  include <utility>
#  include <vector>

#  ifdef DEAL_II_WITH_TBB
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#    include <tbb/task_group.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#  endif

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup threads */
/*@{*/


/**
 * A namespace for the implementation of thread management in deal.II. Most of
 * the content of this namespace is discussed in detail in one of the reports
 * linked to from the documentation page of deal.II.
 *
 * @ingroup threads
 */
namespace Threads
{
  /**
   * A class implementing a <a
   * href="https://en.wikipedia.org/wiki/Lock_(computer_science)">mutex</a>.
   * Mutexes are used to lock data structures to ensure that only a
   * single thread of execution can access them at the same time.
   *
   * This class is a thin wrapper around `std::mutex`. The only difference
   * is that this class is copyable when `std::mutex` is not.  Indeed, when
   * copied, the receiving object does not copy any state from the object
   * being copied, i.e. an entirely new mutex is created. These semantics
   * are consistent with the common use case if a mutex is used as a member
   * variable to lock the other member variables of a class: in that case,
   * the mutex of the copied-to object should only guard the members of the
   * copied-to object, not the members of both the copied-to and
   * copied-from object. Since at the time when the class is copied, the
   * destination's member variable is not used yet, its corresponding mutex
   * should also remain in its original state.
   */
  class Mutex : public std::mutex
  {
  public:
    /**
     * Default constructor.
     */
    Mutex() = default;

    /**
     * Copy constructor. As discussed in this class's documentation, no state
     * is copied from the object given as argument.
     */
    Mutex(const Mutex &)
      : std::mutex()
    {}

    /**
     * Copy operators. As discussed in this class's documentation, no state
     * is copied from the object given as argument.
     */
    Mutex &
    operator=(const Mutex &)
    {
      return *this;
    }
  };
} // namespace Threads


namespace Threads
{
  /**
   * Split the range <code>[begin,end)</code> into <code>n_intervals</code>
   * subintervals of equal size. The last interval will be a little bit
   * larger, if the number of elements in the whole range is not exactly
   * divisible by <code>n_intervals</code>. The type of the iterators has to
   * fulfill the requirements of a forward iterator, i.e.
   * <code>operator++</code> must be available, and of course it must be
   * assignable.
   *
   * A list of subintervals is returned as a vector of pairs of iterators,
   * where each pair denotes the range <code>[begin[i],end[i])</code>.
   *
   * @ingroup threads
   */
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator, ForwardIterator>>
  split_range(const ForwardIterator &begin,
              const ForwardIterator &end,
              const unsigned int     n_intervals);

  /**
   * Split the interval <code>[begin,end)</code> into subintervals of (almost)
   * equal size. This function works mostly as the one before, with the
   * difference that instead of iterators, now values are taken that define
   * the whole interval.
   *
   * @ingroup threads
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  split_interval(const unsigned int begin,
                 const unsigned int end,
                 const unsigned int n_intervals);

  /**
   * @cond internal
   */

  /**
   * A namespace in which helper functions and the like for the threading
   * subsystem are implemented. The members of this namespace are not meant
   * for public use.
   */
  namespace internal
  {
    /**
     * @internal
     *
     * If in a sub-thread an exception is thrown, it is not propagated to the
     * main thread. Therefore, the exception handler that is provided by the
     * applications main function or some of its other parts will not be able
     * to catch these exceptions. Therefore, we have to provide an exception
     * handler in the top function of each sub-thread that at least catches
     * the exception and prints some information, rather than letting the
     * operating system to just kill the program without a message. In each of
     * the functions we use as entry points to new threads, we therefore
     * install a try-catch block, and if an exception of type
     * <code>std::exception</code> is caught, it passes over control to this
     * function, which will then provide some output.
     */
    [[noreturn]] void
    handle_std_exception(const std::exception &exc);

    /**
     * @internal
     *
     * Same as above, but the type of the exception is not derived from
     * <code>std::exception</code>, so there is little way to provide
     * something more useful.
     */
    [[noreturn]] void
    handle_unknown_exception();
  } // namespace internal

  /**
   * @endcond
   */

} // namespace Threads

/* ----------- implementation of functions in namespace Threads ---------- */
#  ifndef DOXYGEN
namespace Threads
{
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator, ForwardIterator>>
  split_range(const ForwardIterator &begin,
              const ForwardIterator &end,
              const unsigned int     n_intervals)
  {
    using IteratorPair = std::pair<ForwardIterator, ForwardIterator>;

    // in non-multithreaded mode, we often have the case that this
    // function is called with n_intervals==1, so have a shortcut here
    // to handle that case efficiently

    if (n_intervals == 1)
      return (std::vector<IteratorPair>(1, IteratorPair(begin, end)));

    // if more than one interval requested, do the full work
    const unsigned int n_elements              = std::distance(begin, end);
    const unsigned int n_elements_per_interval = n_elements / n_intervals;
    const unsigned int residual                = n_elements % n_intervals;

    std::vector<IteratorPair> return_values(n_intervals);

    return_values[0].first = begin;
    for (unsigned int i = 0; i < n_intervals; ++i)
      {
        if (i != n_intervals - 1)
          {
            return_values[i].second = return_values[i].first;
            // note: the cast is performed to avoid a warning of gcc
            // that in the library `dist>=0' is checked (dist has a
            // template type, which here is unsigned if no cast is
            // performed)
            std::advance(return_values[i].second,
                         static_cast<signed int>(n_elements_per_interval));
            // distribute residual in division equally among the first
            // few subintervals
            if (i < residual)
              ++return_values[i].second;

            return_values[i + 1].first = return_values[i].second;
          }
        else
          return_values[i].second = end;
      }
    return return_values;
  }
} // namespace Threads

#  endif // DOXYGEN

namespace Threads
{
  namespace internal
  {
    /**
     * @internal
     *
     * Given an arbitrary type RT, store an element of it and grant
     * access to it through functions get() and set(). There are
     * specializations for reference types (which need to be stored as
     * pointers to the object being referenced), and for type void.
     *
     * This function is not dissimilar to the `std::promise`/`std::future`
     * combination of classes. The difference is that a `std::promise`
     * can only be read once via `std::future::get()` (presumably this
     * design is due to the fact that `std::future::get()` can throw
     * an exception previously stored in the `std::promise`). On
     * the other hand, this class makes the result available for
     * as many times as desired. It also doesn't store any exceptions
     * (though they will be forwarded by the classes using the current
     * class).
     */
    template <typename RT>
    struct return_value
    {
    private:
      RT value;

    public:
      using reference_type = RT &;

      inline return_value()
        : value()
      {}

      inline reference_type
      get()
      {
        return value;
      }

      inline void
      set(RT &&v)
      {
        value = std::move(v);
      }

      inline void
      set_from(std::future<RT> &v)
      {
        value = std::move(v.get());
      }
    };


    /**
     * @internal
     *
     * Given an arbitrary type RT, store an element of it and grant access to
     * it through functions get() and set(). This is the specialization for
     * reference types: since references cannot be set after construction time,
     * we store a pointer instead, which holds the address of the object being
     * referenced.
     *
     * This function is not dissimilar to the `std::promise`/`std::future`
     * combination of classes. The difference is that a `std::promise`
     * can only be read once via `std::future::get()` (presumably this
     * design is due to the fact that `std::future::get()` can throw
     * an exception previously stored in the `std::promise`). On
     * the other hand, this class makes the result available for
     * as many times as desired. It also doesn't store any exceptions
     * (though they will be forwarded by the classes using the current
     * class).
     */
    template <typename RT>
    struct return_value<RT &>
    {
    private:
      RT *value;

    public:
      using reference_type = RT &;

      inline return_value()
        : value(nullptr)
      {}

      inline reference_type
      get() const
      {
        return *value;
      }

      inline void
      set(RT &v)
      {
        value = &v;
      }

      inline void
      set_from(std::future<RT &> &v)
      {
        value = &v.get();
      }
    };


    /**
     * @internal
     *
     * Given an arbitrary type RT, store an element of it and grant access to
     * it through functions get() and set(). This is the specialization for
     * type void: there is obviously nothing to store, so no function set(),
     * and a function get() that returns void.
     *
     * This function is not dissimilar to the `std::promise`/`std::future`
     * combination of classes. The difference is that a `std::promise`
     * can only be read once via `std::future::get()` (presumably this
     * design is due to the fact that `std::future::get()` can throw
     * an exception previously stored in the `std::promise`). On
     * the other hand, this class makes the result available for
     * as many times as desired. It also doesn't store any exceptions
     * (though they will be forwarded by the classes using the current
     * class).
     */
    template <>
    struct return_value<void>
    {
      using reference_type = void;

      static inline void
      get()
      {}


      inline void
      set_from(std::future<void> &)
      {}
    };
  } // namespace internal



  namespace internal
  {
    template <typename RT>
    inline void
    call(const std::function<RT()> & function,
         internal::return_value<RT> &ret_val)
    {
      ret_val.set(function());
    }


    inline void
    call(const std::function<void()> &function, internal::return_value<void> &)
    {
      function();
    }
  } // namespace internal



  namespace internal
  {
    /**
     * A class that represents threads. For each thread, we create exactly one
     * of these objects -- exactly one because it carries the returned value
     * of the function called on the thread.
     *
     * While we have only one of these objects per thread, several
     * Threads::Thread objects may refer to this descriptor. If all Thread
     * objects go out of scope the ThreadDescriptor will detach from the
     * thread before being destroyed.
     */
    template <typename RT>
    struct ThreadDescriptor
    {
      /**
       * An object that represents the thread started.
       */
      std::thread thread;

      /**
       * An object that will hold the value returned by the function called on
       * the thread.
       *
       * The return value is stored in a shared_ptr because we might abandon
       * the ThreadDescriptor.  This makes sure the object stays alive
       * until the thread exits.
       */
      std::shared_ptr<return_value<RT>> ret_val;

      /**
       * An atomic  bool variable that is initially false, is set to true
       * when a new thread is started, and is set back to false once join()
       * has been called.
       *
       * We use this variable to make sure we can call join() twice on the
       * same thread. For some reason, the C++ standard library throws a
       * std::system_error exception if one tries to call std::thread::join
       * twice (and in fact, before the second call, std::thread::joinable
       * returns false) but this is a somewhat desirable thing to do because
       * one doesn't have to keep track whether join() has been called before.
       * Using this variable, whenever we have called join() before, the
       * variable is set to true and we can skip over calling
       * std::thread::join() a second time. Access to this variable is guarded
       * by the following mutex.
       *
       * @note Historically, we did not need the mutex for this variable:
       * threads can only be joined from the thread that created it
       * originally. Consequently, everything that happens in a function that
       * does not create threads (such as the join() function below) looks
       * atomic to the outside world. Since we clear and test thread_is_active
       * in the same function as we call std::thread::join, these actions are
       * atomic and need no mutex. Of course, two threads may call join() on
       * the same thread object at the same time, but this action is undefined
       * anyway since they can not both join the same thread. That said, more
       * recent C++ standards do not appear to have the requirement any more
       * that the only thread that can call join() is the one that created the
       * thread. Neither does `pthread_join` appear to have this requirement any
       * more.  Consequently, we can in fact join from different threads and
       * we test this in base/thread_validity_07.
       *
       * @note The reason why we need to use a std::atomic<bool> is
       * discussed in detail in the documentation of
       * Task::task_has_finished.
       */
      std::atomic<bool> thread_is_active;

      /**
       * Mutex guarding access to the previous variable.
       */
      Mutex thread_is_active_mutex;

      /**
       * Default constructor.
       */
      ThreadDescriptor()
        : thread_is_active(false)
      {}

      ~ThreadDescriptor()
      {
        if (!thread_is_active)
          return;
        thread.detach();
        thread_is_active = false;
      }

      /**
       * Start the thread and let it put its return value into the ret_val
       * object.
       */
      void
      start(const std::function<RT()> &function)
      {
        thread_is_active = true;
        ret_val          = std::make_shared<return_value<RT>>();
        thread           = std::thread(thread_entry_point, function, ret_val);
      }


      /**
       * Wait for the thread to end.
       */
      void
      join()
      {
        // see if the thread hasn't been joined yet. if it has, then
        // join() is a no-op. use schmidt's double-checking strategy
        // to use the mutex only when necessary
        if (thread_is_active == false)
          return;

        std::lock_guard<std::mutex> lock(thread_is_active_mutex);
        if (thread_is_active == true)
          {
            Assert(thread.joinable(), ExcInternalError());
            thread.join();
            thread_is_active = false;
          }
      }

    private:
      /**
       * The function that runs on the thread.
       */
      static void
      thread_entry_point(const std::function<RT()> &       function,
                         std::shared_ptr<return_value<RT>> ret_val)
      {
        // call the function in question. since an exception that is
        // thrown from one of the called functions will not propagate
        // to the main thread, it will kill the program if not treated
        // here before we return to the operating system's thread
        // library
        try
          {
            call(function, *ret_val);
          }
        catch (const std::exception &exc)
          {
            internal::handle_std_exception(exc);
          }
        catch (...)
          {
            internal::handle_unknown_exception();
          }
      }
    };
  } // namespace internal


  /**
   * An object that represents a spawned thread. This object can be freely
   * copied around in user space, and all instances will represent the same
   * thread and can require to wait for its termination and access its return
   * value.
   *
   * Threads can be abandoned, i.e. if you just call Threads::new_thread but
   * don't care about the returned object, or if you assign the return
   * Threads::Thread object to an object that subsequently goes out of scope,
   * then the thread previously created will still continue to do work. You
   * will simply not be able to access its return value any more, and it may
   * also happen that your program terminates before the thread has finished
   * its work.
   *
   * The default value of the template argument is <code>void</code>, so if
   * the function you are calling on a new thread has no return value, you can
   * omit the template argument.
   *
   * @ingroup threads
   *
   * @deprecated Use std::thread or std::jthread instead.
   *
   * @note Since this class is used in ThreadGroup, its constructors, rather
   * than the class itself, are deprecated to allow compilation with
   * -Werror=deprecated-declarations.
   */
  template <typename RT = void>
  class Thread
  {
  public:
    /**
     * Construct a thread object with a function object.
     */
    DEAL_II_DEPRECATED
    Thread(const std::function<RT()> &function)
      : thread_descriptor(new internal::ThreadDescriptor<RT>())
    {
      // in a second step, start the thread.
      thread_descriptor->start(function);
    }

    /**
     * Default constructor. You can't do much with a thread object constructed
     * this way, except for assigning it a thread object that holds data
     * created by the new_thread() functions.
     */
    DEAL_II_DEPRECATED
    Thread() = default;

    /**
     * Copy constructor.
     */
    DEAL_II_DEPRECATED
    Thread(const Thread<RT> &t)
      : thread_descriptor(t.thread_descriptor)
    {}

    /**
     * Join the thread represented by this object, i.e. wait for it to finish.
     * If you have used the default constructor of this class and have not
     * assigned a thread object to it, then this function is a no-op.
     */
    void
    join() const
    {
      if (thread_descriptor)
        thread_descriptor->join();
    }

    /**
     * Get the return value of the function of the thread. Since it
     * is only available once the thread finishes, this function
     * internally also calls join(). You can call this function
     * multiple times as long as the object refers to the same task,
     * and expect to get the same return value every time. (With the
     * exception of the case where the returned object has been moved;
     * see below.)
     *
     * @note The function returns a <i>non-@p const reference</i> to
     * the returned object, instead of the returned object. This
     * allows writing code such as
     * @code
     *   Threads::Thread<int> t = Threads::new_thread (
     *     ...function returning an int...);
     *   t.return_value() = 42;      // overwrite returned value
     *   int i = t.return_value();   // i is now 42
     * @endcode
     * You will rarely have a need to write such code. On the other hand,
     * the function needs to return a writable (non-@p const) reference to
     * support code such as this:
     * @code
     *   std::unique_ptr<int> create_int (const std::string &s)
     *   {
     *     ...
     *   }
     *
     *   void f()
     *   {
     *     Threads::Thread<std::unique_ptr<int>>
     *       t = Threads::new_thread (&create_int, "42");
     *
     *     std::unique_ptr<int> i = std::move(t.return_value());
     *     ...
     *   }
     * @endcode
     * Here, it is necessary to `std::move` the returned object (namely,
     * the <code>std::unique_ptr</code> object) because
     * <code>std::unique_ptr</code> objects can not be copied. In other words,
     * to get the pointer out of the object returned from the thread, it needs
     * to be moved, and in order to be moved, the current function needs to
     * return a writable (non-@p const) reference.
     */
    typename internal::return_value<RT>::reference_type
    return_value()
    {
      join();
      return thread_descriptor->ret_val->get();
    }

    /**
     * Return true if this object has had a thread associated with it, either
     * by using the non-default constructor or by assignment.
     */
    bool
    valid() const
    {
      return static_cast<bool>(thread_descriptor);
    }


    /**
     * Check for equality of thread objects. Since objects of this class store
     * an implicit pointer to an object that exists exactly once for each
     * thread, the check is simply to compare these pointers.
     */
    bool
    operator==(const Thread &t) const
    {
      return thread_descriptor == t.thread_descriptor;
    }

  private:
    /**
     * Shared pointer to the object representing the thread, and abstracting
     * operating system functions to work on it. This also makes sure that the
     * object lives as long as there is at least one subscriber to it.
     */
    std::shared_ptr<internal::ThreadDescriptor<RT>> thread_descriptor;
  };


  namespace internal
  {
    /**
     * A general template that returns std::ref(t) if t is of reference
     * type, and t otherwise.
     *
     * The case that t is of reference type is handled in a partial
     * specialization declared below.
     */
    template <typename T>
    struct maybe_make_ref
    {
      static T
      act(T &t)
      {
        return t;
      }
    };



    /**
     * A general template that returns std::ref(t) if t is of reference
     * type, and t otherwise.
     *
     * The case that t is of reference type is handled in this partial
     * specialization.
     */
    template <typename T>
    struct maybe_make_ref<T &>
    {
      static std::reference_wrapper<T>
      act(T &t)
      {
        return std::ref(t);
      }
    };
  } // namespace internal



  // ----------- thread starters for functions not taking any parameters

  /**
   * Overload of the new_thread function for objects that can be converted to
   * std::function<RT ()>, i.e. anything that can be called like a
   * function object without arguments and returning an object of type RT (or
   * void).
   *
   * @deprecated Use std::thread or std::jthread instead.
   *
   * @ingroup threads
   */
  template <typename RT>
  DEAL_II_DEPRECATED inline Thread<RT>
  new_thread(const std::function<RT()> &function)
  {
    // Here and below we need to disable deprecation warnings for calling the
    // constructor in this function - as this function itself is deprecated
    // these warnings are not helpful. This problem only appears in some
    // configurations (e.g., Debian 11 with GCC-10).
    DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
    return Thread<RT>(function);
    DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
  }



  /**
   * Overload of the new_thread() function for objects that can be called like a
   * function object without arguments. In particular, this function allows
   * calling Threads::new_thread() with either objects that result from using
   * std::bind, or using lambda functions. For example, this function is called
   * when writing code such as
   * @code
   * Threads::Thread<int>
   *   thread = Threads::new_thread ( [] () {
   *                                          do_this();
   *                                          then_do_that();
   *                                          return 42;
   *                                        });
   * @endcode
   * Here, we run the sequence of functions
   * <code>do_this()</code> and <code>then_do_that()</code> on
   * a separate thread, by making the lambda function declared here the
   * function to execute on the thread. The lambda function then returns
   * 42 (which is a bit pointless here, but it could of course be some
   * computed number), and this is going to be the returned value you
   * can later retrieve via <code>thread.return_value()</code> once the
   * thread (i.e., the body of the lambda function) has completed.
   *
   * @note Every lambda function (or whatever else it is you pass to
   *   the new_thread() function here, for example the result of a
   *   std::bind() expression) has a return type and consequently
   *   returns an object of this type. This type can be inferred
   *   using the C++11 <code>decltype</code> statement used in the
   *   declaration of this function, and it is then used as the template
   *   argument of the Threads::Thread object returned by the current function.
   *   In the example above, because the lambda function returns 42
   *   (which in C++ has data type <code>int</code>), the inferred
   *   type is <code>int</code> and the task object will have type
   *   <code>Task@<int@></code>. In other words, it is not <i>necessary</i>
   *   to explicitly specify in user code what that return type
   *   of the lambda or std::bind expression will be, though it is
   *   possible to explicitly do so by (entirely equivalently) writing
   *   @code
   *   Threads::Thread<int>
   *     thread = Threads::new_thread ( [] () -> int {
   *                                                   do_this();
   *                                                   then_do_that();
   *                                                   return 42;
   *                                                 });
   *   @endcode
   *
   * @note In practice, the lambda functions you will pass to
   *   new_thread() will of course typically be more complicated.
   *   In particular, they will likely <i>capture</i> variables
   *   from the surrounding context and use them within the lambda.
   *   See
   * https://en.wikipedia.org/wiki/Anonymous_function#C.2B.2B_.28since_C.2B.2B11.29
   *   for more on how lambda functions work.
   *
   * @note If you pass a lambda function as an argument to the
   *   current function that captures a variable <i>by reference</i>,
   *   or if you use a std::bind that binds a function argument to
   *   a reference variable using std::ref() or std::cref(), then
   *   obviously you can only do this if the variables you reference
   *   or capture have a lifetime that extends at least until the time
   *   where the thread finishes.
   *
   * @deprecated Use std::thread or std::jthread instead.
   *
   * @ingroup CPP11
   */
  template <typename FunctionObjectType>
  DEAL_II_DEPRECATED inline auto
  new_thread(FunctionObjectType function_object)
    -> Thread<decltype(function_object())>
  {
    // See the comment in the first new_thread() implementation
    DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
    using return_type = decltype(function_object());
    return Thread<return_type>(std::function<return_type()>(function_object));
    DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
  }



  /**
   * Overload of the new_thread function for non-member or static member
   * functions.
   *
   * @deprecated Use std::thread or std::jthread instead.
   *
   * @ingroup threads
   */
  template <typename RT, typename... Args>
  DEAL_II_DEPRECATED inline Thread<RT>
  new_thread(RT (*fun_ptr)(Args...), typename identity<Args>::type... args)
  {
    // See the comment in the first new_thread() implementation
    DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
    auto dummy = std::make_tuple(internal::maybe_make_ref<Args>::act(args)...);
    return new_thread(
      [dummy, fun_ptr]() -> RT { return std_cxx17::apply(fun_ptr, dummy); });
    DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
  }



  /**
   * Overload of the non-const new_thread function for member functions.
   *
   * @deprecated Use std::thread or std::jthread instead.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename... Args>
  DEAL_II_DEPRECATED inline Thread<RT>
  new_thread(RT (C::*fun_ptr)(Args...),
             typename identity<C>::type &c,
             typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_thread(std::function<RT()>(std::bind(
      fun_ptr, std::ref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }

  /**
   * Overload of the new_thread function for const member functions.
   *
   * @deprecated Use std::thread or std::jthread instead.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename... Args>
  DEAL_II_DEPRECATED inline Thread<RT>
  new_thread(RT (C::*fun_ptr)(Args...) const,
             typename identity<const C>::type &c,
             typename identity<Args>::type... args)
  {
    // See the comment in the first new_thread() implementation
    DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_thread(std::function<RT()>(std::bind(
      fun_ptr, std::cref(c), internal::maybe_make_ref<Args>::act(args)...)));
    DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
  }

  // ------------------------ ThreadGroup -------------------------------------

  /**
   * A container for thread objects. Allows to add new thread objects and wait
   * for them all together. The thread objects need to have the same return
   * value for the called function.
   *
   * @ingroup threads
   *
   * @deprecated Use TaskGroup instead.
   */
  template <typename RT = void>
  class DEAL_II_DEPRECATED ThreadGroup
  {
  public:
    /**
     * Add another thread object to the collection.
     */
    ThreadGroup &
    operator+=(const Thread<RT> &t)
    {
      threads.push_back(t);
      return *this;
    }

    /**
     * Wait for all threads in the collection to finish. It is not a problem
     * if some of them have already been waited for, i.e. you may call this
     * function more than once, and you can also add new thread objects
     * between subsequent calls to this function if you want.
     */
    void
    join_all() const
    {
      for (auto &t : threads)
        t.join();
    }

  private:
    /**
     * List of thread objects.
     */
    std::list<Thread<RT>> threads;
  };


  namespace internal
  {
    /**
     * Set the value of a std::promise object by evaluating the action.
     */
    template <typename RT, typename Function>
    void
    evaluate_and_set_promise(Function &function, std::promise<RT> &promise)
    {
      promise.set_value(function());
    }


    /**
     * Set the value of a std::promise object by evaluating the
     * action. This function is a specialization of the previous one
     * for the case where the return type is `void`. Consequently, we
     * can't set a value. But we do evaluate the function object and
     * call `std::promise::set_value()` without argument.
     */
    template <typename Function>
    void
    evaluate_and_set_promise(Function &function, std::promise<void> &promise)
    {
      function();
      promise.set_value();
    }
  } // namespace internal



  /**
   * This class describes a task object, i.e., what one obtains by calling
   * Threads::new_task(). The idea is that Threads::new_task() allows one to run
   * a function whenever the C++ run-time system finds it convenient --
   * typically, when there is an idle processor available. This can be used to
   * run things in the background when there is no immediate need for the
   * result, or if there are other things that could well be done in parallel.
   * Whenever the result of that background task is needed, one can call either
   * join() to just wait for the task to finish, or return_value() to obtain the
   * value that was returned by the function that was run on that background
   * task.
   *
   * This class is conceptually similar to the
   * [`std::future`](https://en.cppreference.com/w/cpp/thread/future) class that
   * is returned by
   * [`std::async`](https://en.cppreference.com/w/cpp/thread/async) (which is
   * itself similar to what Threads::new_task() does). The principal conceptual
   * difference is that one can only call `std::future::get()` once, whereas one
   * can call Threads::Task::return_value() as many times as desired. It is,
   * thus, comparable to the
   * [`std::shared_future`](https://en.cppreference.com/w/cpp/thread/shared_future)
   * class. However, `std::shared_future` can not be used for types that can not
   * be copied -- a particular restriction for `std::unique_ptr`, for example.
   *
   * @ingroup threads
   */
  template <typename RT = void>
  class Task
  {
  public:
    /**
     * Construct a task object, given a function object to execute on
     * the task, and then schedule this function for
     * execution. However, when MultithreadInfo::n_threads() returns
     * 1, i.e., if the deal.II runtime system has been configured to
     * only use one thread, then just execute the given function
     * object.
     *
     * @post Using this constructor automatically makes the task object
     * joinable().
     */
    Task(const std::function<RT()> &function_object)
    {
      if (MultithreadInfo::n_threads() > 1)
        {
#  ifdef DEAL_II_WITH_TBB
          // Create a promise object and from it extract a future that
          // we can use to refer to the outcome of the task. For reasons
          // explained below, we can't just create a std::promise object,
          // but have to make do with a pointer to such an object.
          std::unique_ptr<std::promise<RT>> promise =
            std::make_unique<std::promise<RT>>();
          task_data =
            std::make_shared<TaskData>(std::move(promise->get_future()));

          // Then start the task, using a task_group object (for just this one
          // task) that is associated with the TaskData object. Note that we
          // have to *copy* the function object being executed so that it is
          // guaranteed to live on the called thread as well -- the copying is
          // facilitated by capturing the 'function_object' variable by value.
          //
          // We also have to *move* the promise object into the new task's
          // memory space because promises can not be copied and we can't refer
          // to it by reference because it's a local variable of the current
          // (surrounding) function that may go out of scope before the promise
          // is ultimately set. This leads to a conundrum: if we had just
          // declared 'promise' as an object of type std::promise, then we could
          // capture it in the lambda function via
          //     [..., promise=std::move(promise)]() {...}
          // and set the promise in the body of the lambda. But setting a
          // promise is a non-const operation on the promise, and so we would
          // actually have to declare the lambda function as 'mutable' because
          // by default, lambda captures are 'const'. That is, we would have
          // to write
          //     [..., promise=std::move(promise)]() mutable {...}
          // But this leads to other problems: It turns out that the
          // tbb::task_group::run() function cannot take mutable lambdas as
          // argument :-(
          //
          // We work around this issue by not declaring the 'promise' variable
          // as an object of type std::promise, but as a pointer to such an
          // object. This pointer we can move, and the *pointer* itself can
          // be 'const' (meaning we can leave the lambda as non-mutable)
          // even though we modify the object *pointed to*. One would think
          // that a std::unique_ptr would be the right choice for this, but
          // that's not true: the resulting lambda function can then be
          // non-mutable, but the lambda function object is not copyable
          // and at least some TBB variants require that as well. So
          // instead we move the std::unique_ptr used above into a
          // std::shared_ptr to be stored within the lambda function object.
          task_data->task_group.run(
            [function_object,
             promise =
               std::shared_ptr<std::promise<RT>>(std::move(promise))]() {
              try
                {
                  internal::evaluate_and_set_promise(function_object, *promise);
                }
              catch (...)
                {
                  try
                    {
                      // store anything thrown in the promise
                      promise->set_exception(std::current_exception());
                    }
                  catch (...)
                    {
                      // set_exception() may throw too. But ignore this on
                      // the task.
                    }
                }
            });

#  else
          // If no threading library is supported, just fall back onto C++11
          // facilities. The problem with this is that the standard does
          // not actually say what std::async should do. The first
          // argument to that function can be std::launch::async or
          // std::launch::deferred, or both. The *intent* of the standard's
          // authors was probably that if one sets it to
          //   std::launch::async | std::launch::deferred,
          // that the task is run in a thread pool. But at least as of
          // 2021, GCC doesn't do that: It just runs it on a new thread.
          // If one chooses std::launch::deferred, it runs the task on
          // the same thread but only when one calls join() on the task's
          // std::future object. In the former case, this leads to
          // oversubscription, in the latter case to undersubscription of
          // resources. We choose oversubscription here.
          //
          // The issue illustrates why relying on external libraries
          // with task schedulers is the way to go.
          task_data = std::make_shared<TaskData>(
            std::async(std::launch::async | std::launch::deferred,
                       function_object));
#  endif
        }
      else
        {
          // Only one thread allowed. So let the task run to completion
          // and just emplace a 'ready' future.
          //
          // The design of std::promise/std::future is unclear, but it
          // seems that the intent is to obtain the std::future before
          // we set the std::promise. So create the TaskData object at
          // the top and then run the task and set the returned
          // value. Since everything here happens sequentially, it
          // really doesn't matter in which order all of this is
          // happening.
          std::promise<RT> promise;
          task_data = std::make_shared<TaskData>(promise.get_future());
          try
            {
              internal::evaluate_and_set_promise(function_object, promise);
            }
          catch (...)
            {
              try
                {
                  // store anything thrown in the promise
                  promise.set_exception(std::current_exception());
                }
              catch (...)
                {
                  // set_exception() may throw too. But ignore this on
                  // the task.
                }
            }
        }
    }

    /**
     * Default constructor. You can't do much with a task object constructed
     * this way, except for assigning it a task object that holds data created
     * by the Threads::new_task() functions.
     *
     * @post Using this constructor leaves the object in an unjoinable state,
     * i.e., joinable() will return false.
     */
    Task() = default;

    /**
     * Join the task represented by this object, i.e. wait for it to finish.
     *
     * A task can be joined multiple times (while the first join() operation
     * may block until the task has completed running, all successive attempts
     * to join will return immediately).
     *
     * If the operation that was executed on the task with which this
     * object was initialized throws an exception instead of returning
     * regularly, then calling the current join() function will first
     * wait for that task to finish, and then in turn throw the
     * exception that the task operation had thrown originally. This
     * allows for the propagation of exceptions from tasks executed on
     * a separate thread to the calling thread.
     *
     * (This behavior differs from that of
     * [`std::future`](https://en.cppreference.com/w/cpp/thread/future),
     * where the `std::future::wait()` function only waits for
     * completion of the operation, whereas the exception is
     * propagated only once one calls `std::future::get()`. However,
     * this is awkward when putting `void` functions onto separate
     * tasks because these do not actually return anything;
     * consequently, it is more natural to call `std::task::wait()`
     * for such tasks than the `std::task::get()` function since the
     * latter does not, actually, return anything that could be
     * gotten.)
     *
     * @pre You can't call this function if you have used the default
     * constructor of this class and have not assigned a task object to it. In
     * other words, the function joinable() must return true.
     */
    void
    join() const
    {
      // Make sure we actually have a task that we can wait for.
      AssertThrow(joinable(), ExcNoTask());

      task_data->wait();
    }

    /**
     * Return whether the current object can be joined. You can join a task
     * object once a task (typically created with Threads::new_task()) has
     * actually been assigned to it. On the other hand, the function returns
     * false if the object has been default constructed.
     *
     * A task can be joined multiple times (while the first join() operation
     * may block until the task has completed running, all successive attempts
     * to join will return immediately). Consequently, if this function
     * returns true, it will continue to return true until the task object it
     * reports on is assigned to from another object.
     */
    bool
    joinable() const
    {
      return (task_data != nullptr);
    }


    /**
     * Get the return value of the function of the task. Since it is
     * only available once the thread finishes, this function
     * internally also calls join(). You can call this function
     * multiple times as long as the object refers to the same task,
     * and expect to get the same return value every time. (With the
     * exception of the case where the returned object has been moved;
     * see below.)
     *
     * @note The function returns a <i>non-@p const reference</i> to
     * the returned object, instead of the returned object. This
     * allows writing code such as
     * @code
     *   Threads::Task<int> t = Threads::new_task (...function returning an
     * int...); t.return_value() = 42;      // overwrite returned value int i =
     * t.return_value();   // i is now 42
     * @endcode
     * You will rarely have a need to write such code. On the other hand,
     * the function needs to return a writable (non-@p const) reference to
     * support code such as this:
     * @code
     *   std::unique_ptr<int> create_int (const std::string &s) { ... }
     *
     *   void f()
     *   {
     *     Threads::Task<std::unique_ptr<int>>
     *       t = Threads::new_task (&create_int, "42");
     *
     *     std::unique_ptr<int> i = std::move(t.return_value());
     *     ...
     *   }
     * @endcode
     * Here, it is necessary to `std::move` the returned object (namely,
     * the <code>std::unique_ptr</code> object) because
     * <code>std::unique_ptr</code> objects can not be copied. In other words,
     * to get the pointer out of the object returned from the task, it needs
     * to be moved, and in order to be moved, the current function needs to
     * return a writable (non-@p const) reference.
     *
     * This function internally calls the join() member function. As a
     * consequence, and as explained there, if the packaged task
     * throws an exception that is then re-thrown by the join()
     * function and consequently also the current function if you have
     * not previously called join().
     *
     * @pre You can't call this function if you have used the default
     * constructor of this class and have not assigned a task object to it. In
     * other words, the function joinable() must return true.
     */
    typename internal::return_value<RT>::reference_type
    return_value()
    {
      // Make sure we actually have a task that we can wait for.
      AssertThrow(joinable(), ExcNoTask());

      // Then return the promised object. If necessary, wait for the promise to
      // be set.
      return task_data->get();
    }


    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception
     */
    DeclExceptionMsg(ExcNoTask,
                     "The current object is not associated with a task that "
                     "can be joined. It may have been detached, or you "
                     "may have already joined it in the past.");
    //@}
  private:
    /**
     * A data structure that holds a std::future into which the task deposits
     * its return value. Since one can only call std::future::get() once,
     * we do so in the get() member function and then move the returned object
     * into the `returned_object` member variable from where we can read it
     * multiple times and from where it can also be moved away if it is not
     * copyable.
     */
    class TaskData
    {
    public:
      /**
       * Constructor. Initializes an std::future object and assumes
       * that the task so set has not finished yet.
       */
      TaskData(std::future<RT> &&future) noexcept
        : future(std::move(future))
        , task_has_finished(false)
      {}

      /**
       * There can only be one TaskData object referencing
       * a task. Make sure that these objects are not copied.
       */
      TaskData(const TaskData &) = delete;

      /**
       * There can only be one TaskData object referencing
       * a task. Make sure that these objects are not moved.
       */
      TaskData(TaskData &&) = delete;

      /**
       * There can only be one TaskData object referencing
       * a task. Make sure that these objects are not copied.
       */
      TaskData &
      operator=(const TaskData &) = delete;

      /**
       * There can only be one TaskData object referencing
       * a task. Make sure that these objects are not moved.
       */
      TaskData &
      operator=(TaskData &&) = delete;

      /**
       * Destructor. Wait for the results to be ready. This ensures that the
       * last Task object holding a shared pointer to the current TaskData
       * object blocks until the task has actually finished -- in essence,
       * this makes sure that one cannot just abandon a task completely
       * by letting all Task objects that point to it go out of scope.
       */
      ~TaskData() noexcept
      {
        // Explicitly wait for the results to be ready. This class stores
        // a std::future object, and we could just let the compiler generate
        // the destructor which would then call the destructor of std::future
        // which *may* block until the future is ready. As explained in
        //   https://en.cppreference.com/w/cpp/thread/future/~future
        // this is only a *may*, not a *must*. (The standard does not
        // appear to say anything about it at all.) As a consequence,
        // let's be explicit about waiting.
        //
        // One of the corner cases we have to worry about is that if a task
        // ends by throwing an exception, then wait() will re-throw that
        // exception on the thread that calls it, the first time around
        // someone calls wait() (or the return_value() function of the
        // surrounding class). So if we get to this constructor and an exception
        // is thrown by wait(), then that means that the last Task object
        // referring to a task is going out of scope with nobody having
        // ever checked the return value of the task itself. In that case,
        // one could argue that they would also not have cared about whether
        // an exception is thrown, and that we should simply ignore the
        // exception. This is what we do here. It is also the simplest solution,
        // because we don't know what one should do with the exception to begin
        // with: destructors aren't allowed to throw exceptions, so we can't
        // just rethrow it here if one had been triggered.
        try
          {
            wait();
          }
        catch (...)
          {}
      }

      /**
       * Wait for the std::future object to be ready, i.e., for the
       * time when the std::promise receives its value. If this has
       * already happened, this function can follow a fast path.
       */
      void
      wait()
      {
        // If we have previously already moved the result, then we don't
        // need a lock and can just return.
        if (task_has_finished)
          return;

        // Else, we need to go under a lock and try again. A different thread
        // may have waited and finished the task since then, so we have to try
        // a second time. (This is Schmidt's double-checking pattern.)
        std::lock_guard<std::mutex> lock(mutex);
        if (task_has_finished)
          return;
        else
          {
#  ifdef DEAL_II_WITH_TBB
            // If we build on the TBB, then we can't just wait for the
            // std::future object to get ready. Apparently the TBB happily
            // enqueues a task into an arena and then just sits on it without
            // ever executing it unless someone expresses an interest in the
            // task. The way to avoid this is to add the task to a
            // tbb::task_group, and then here wait for the single task
            // associated with that task group.
            task_group.wait();
#  endif

            // Wait for the task to finish and then move its
            // result. (We could have made the set_from() function
            // that we call here wait for the future to be ready --
            // which happens implicitly when it calls future.get() --
            // but that would have required putting an explicit
            // future.wait() into the implementation of
            // internal::return_value<void>::set_from(), which is a
            // bit awkward: that class doesn't actually need to set
            // anything, and so it looks odd to have the explicit call
            // to future.wait() in the set_from() function. Avoid the
            // issue by just explicitly calling future.wait() here.)
            future.wait();
            returned_object.set_from(future);

            // Now we can safely set the flag and return.
            task_has_finished = true;
          }
      }



      typename internal::return_value<RT>::reference_type
      get()
      {
        wait();
        return returned_object.get();
      }

    private:
      /**
       * A mutex used to synchronize access to the data structures of this
       * class.
       */
      std::mutex mutex;

      /**
       * The promise associated with the task that is represented by the current
       * class.
       */
      std::future<RT> future;

      /**
       * A boolean indicating whether the task in question has finished.
       *
       * @note We are using a `std::atomic_bool` here because we have
       * to make sure that concurrent reads and stores between threads are
       * properly synchronized, and that sequential reads on a given thread
       * are not reordered or optimized away. A std::atomic [1] achieves
       * this because (if not otherwise annotated) reads and stores to the
       * boolean are subject to the std::memory_order_seq_cst memory
       * ordering [2]. This ensures that Schmidt's double checking does
       * indeed work. For additional information (and a potentially more
       * efficient implementation) see [3].
       *
       * [1] https://en.cppreference.com/w/cpp/atomic/atomic
       * [2] https://en.cppreference.com/w/cpp/atomic/memory_order
       * [3]
       * https://preshing.com/20130930/double-checked-locking-is-fixed-in-cpp11/
       */
      std::atomic<bool> task_has_finished;

      /**
       * The place where the returned value is moved to once the std::future
       * has delivered.
       */
      internal::return_value<RT> returned_object;

#  ifdef DEAL_II_WITH_TBB
      /**
       * A task group object we can wait for.
       */
      tbb::task_group task_group;

      friend class Task<RT>;
#  endif
    };

    /**
     * A pointer to a descriptor of the object that described the task
     * and its return value.
     */
    std::shared_ptr<TaskData> task_data;
  };



  /**
   * Overload of the new_task function for objects that can be converted to
   * std::function<RT ()>, i.e. anything that can be called like a
   * function object without arguments and returning an object of type RT (or
   * void).
   *
   * @note When MultithreadInfo::n_threads() returns 1, i.e., if the
   *   deal.II runtime system has been configured to only use one
   *   thread, then this function just executes the given function
   *   object immediately and stores the return value in the Task
   *   object returned by this function.
   *
   * @note Threads::new_task() is, in essence, equivalent to calling
   *   `std::async(std::launch::async, ...)` in that it runs the given task
   *   in the background. See https://en.cppreference.com/w/cpp/thread/async
   *   for more information.
   *
   * @ingroup threads
   */
  template <typename RT>
  inline Task<RT>
  new_task(const std::function<RT()> &function)
  {
    return Task<RT>(function);
  }



  /**
   * Overload of the new_task function for objects that can be called like a
   * function object without arguments. In particular, this function allows
   * calling Threads::new_task() with either objects that result from using
   * std::bind, or using lambda functions. For example, this function is called
   * when writing code such as
   * @code
   * Threads::Task<int>
   *   task = Threads::new_task ( [] () {
   *                                      do_this();
   *                                      then_do_that();
   *                                      return 42;
   *                                    });
   * @endcode
   * Here, we schedule the call to the sequence of functions
   * <code>do_this()</code> and <code>then_do_that()</code> on
   * a separate task, by making the lambda function declared here the
   * function to execute on the task. The lambda function then returns
   * 42 (which is a bit pointless here, but it could of course be some
   * computed number), and this is going to be the returned value you
   * can later retrieve via <code>task.return_value()</code> once the
   * task (i.e., the body of the lambda function) has completed.
   *
   * @note When MultithreadInfo::n_threads() returns 1, i.e., if the
   *   deal.II runtime system has been configured to only use one
   *   thread, then this function just executes the given function
   *   object immediately and stores the return value in the Task
   *   object returned by this function.
   *
   * @note Every lambda function (or whatever else it is you pass to
   *   the new_task() function here, for example the result of a
   *   std::bind() expression) has a return type and consequently
   *   returns an object of this type. This type can be inferred
   *   using the C++11 <code>decltype</code> statement used in the
   *   declaration of this function, and it is then used as the template
   *   argument of the Threads::Task object returned by the current function.
   *   In the example above, because the lambda function returns 42
   *   (which in C++ has data type <code>int</code>), the inferred
   *   type is <code>int</code> and the task object will have type
   *   <code>Task@<int@></code>. In other words, it is not <i>necessary</i>
   *   to explicitly specify in user code what that return type
   *   of the lambda or std::bind expression will be, though it is
   *   possible to explicitly do so by (entirely equivalently) writing
   *   @code
   *   Threads::Task<int>
   *     task = Threads::new_task ( [] () -> int {
   *                                               do_this();
   *                                               then_do_that();
   *                                               return 42;
   *                                             });
   *   @endcode
   *
   * @note In practice, the lambda functions you will pass to
   *   new_task() will of course typically be more complicated.
   *   In particular, they will likely <i>capture</i> variables
   *   from the surrounding context and use them within the lambda.
   *   See
   * https://en.wikipedia.org/wiki/Anonymous_function#C.2B.2B_.28since_C.2B.2B11.29
   *   for more on how lambda functions work.
   *
   * @note If you pass a lambda function as an argument to the
   *   current function that captures a variable <i>by reference</i>,
   *   or if you use a std::bind that binds a function argument to
   *   a reference variable using std::ref() or std::cref(), then
   *   obviously you can only do this if the variables you reference
   *   or capture have a lifetime that extends at least until the time
   *   where the task finishes.
   *
   * @note Threads::new_task() is, in essence, equivalent to calling
   *   `std::async(std::launch::async, ...)` in that it runs the given task
   *   in the background. See https://en.cppreference.com/w/cpp/thread/async
   *   for more information.
   *
   * @ingroup CPP11
   */
  template <typename FunctionObjectType>
  inline auto
  new_task(FunctionObjectType function_object)
    -> Task<decltype(function_object())>
  {
    using return_type = decltype(function_object());
    dealii::MultithreadInfo::initialize_multithreading();
    return new_task(std::function<return_type()>(function_object));
  }



  /**
   * Overload of the new_task function for non-member or static member
   * functions. See the other functions of same name for more information.
   *
   * @ingroup threads
   */
  template <typename RT, typename... Args>
  inline Task<RT>
  new_task(RT (*fun_ptr)(Args...), typename identity<Args>::type... args)
  {
    auto dummy = std::make_tuple(internal::maybe_make_ref<Args>::act(args)...);
    return new_task(
      [dummy, fun_ptr]() -> RT { return std_cxx17::apply(fun_ptr, dummy); });
  }



  /**
   * Overload of the non-const new_task function. See the other functions of
   * same name for more information.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename... Args>
  inline Task<RT>
  new_task(RT (C::*fun_ptr)(Args...),
           typename identity<C>::type &c,
           typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_task(std::function<RT()>(std::bind(
      fun_ptr, std::ref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }

  /**
   * Overload of the new_task function. See the other functions of same name for
   * more information.
   *
   * @ingroup threads
   */
  template <typename RT, typename C, typename... Args>
  inline Task<RT>
  new_task(RT (C::*fun_ptr)(Args...) const,
           typename identity<const C>::type &c,
           typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_task(std::function<RT()>(std::bind(
      fun_ptr, std::cref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }


  // ------------------------ TaskGroup -------------------------------------

  /**
   * A container for task objects. Allows to add new task objects and wait for
   * them all together. The task objects need to have the same return value
   * for the called function.
   *
   * Note that the call to join_all() must be executed on the same thread as
   * the calls that add subtasks. Otherwise, there might be a deadlock. In
   * other words, a Task object should never passed on to another task for
   * calling the join() method.
   *
   * @ingroup tasks
   */
  template <typename RT = void>
  class TaskGroup
  {
  public:
    /**
     * Add another task object to the collection.
     */
    TaskGroup &
    operator+=(const Task<RT> &t)
    {
      tasks.push_back(t);
      return *this;
    }


    /**
     * Return how many tasks have been put into this group. This
     * function does not distinguish how many of these tasks have
     * already run and have finished, are still waiting to be
     * scheduled to a CPU resource, or are currently running. Tasks
     * that have been joined already are also still counted.
     */
    std::size_t
    size() const
    {
      return tasks.size();
    }


    /**
     * Wait for all tasks in the collection to finish. It is not a problem if
     * some of them have already been waited for, i.e. you may call this
     * function more than once, and you can also add new task objects between
     * subsequent calls to this function if you want.
     */
    void
    join_all() const
    {
      for (auto &t : tasks)
        t.join();
    }

  private:
    /**
     * List of task objects.
     */
    std::list<Task<RT>> tasks;
  };

} // namespace Threads

/**
 * @}
 */


//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_thread_management_h
#endif
//---------------------------------------------------------------------------

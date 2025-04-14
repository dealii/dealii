// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_task_result_h
#define dealii_task_result_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>

#include <atomic>
#include <mutex>
#include <optional>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup threads
 * @{
 */

namespace Threads
{
  /**
   * A class that represents the outcome of a Threads::Task. The class is used
   * as a member variable or local variable when something is computed on a
   * separate task in the background. For example, a class's constructor may
   * want to set up some expensive data structures that are not typically
   * used right away in the next line of the function that created the object,
   * but will be useful later; in such a case, it may defer computation of that
   * data to a separate task that will be scheduled whenever a CPU core is
   * available, and the result of that deferred computation will become
   * available through a variable of the current type.
   *
   * In some regard, class TaskResult is similar to class Lazy in that it
   * defers construction of an object to a later point while already giving it
   * a home. The difference is simply that TaskResult has the creation job
   * already scheduled whereas Lazy defers creation to the first use. Lazy is
   * therefore a better choice whenever the object referenced may or may not
   * ever be used, say for some obscure functionality of a class that is rarely
   * required. On the other hand, TaskResult is used for member variables that
   * are definitely needed, but perhaps not right away.
   *
   * This class has an interface that is quite similar to the Threads::Task
   * class itself, and object of which is provided to the constructor. The key
   * difference is that the Task class describes the *task* that computes a
   * result, whereas the current TaskResult class describes the *result* that is
   * being computed by the task. The main practical differences between these
   * perspectives are what happens when you want to copy an object: Copying a
   * Task object results in two objects that are both referencing the same task,
   * with the same returned object when the task has finished. On the other
   * hand, copying a TaskResult object after the computing task has finished
   * results in two copies of the returned object.
   *
   * This class can also be compared with
   * [`std::future`](https://en.cppreference.com/w/cpp/thread/future). That
   * class also represents the result of a pending operation, but it lacks the
   * specifics of what kind of operation that is (in particular, it does not
   * know whether a task, a thread, or a `std::packaged_task` will eventually
   * produce the value). As a consequence, it can *wait* for
   * the result to become available, but it lacks the knowledge to detect
   * certain common programming mistakes such as those described in the
   * documentation of the destructor of this class, or of this class's
   * `operator=()`. Another significant difference is that
   * one can only call `std::future::get()` once, whereas one
   * can call Threads::Task::return_value() as many times as desired. It is,
   * thus, comparable to the
   * [`std::shared_future`](https://en.cppreference.com/w/cpp/thread/shared_future)
   * class. However, `std::shared_future` can not be used for types that can not
   * be copied -- a particular restriction for `std::unique_ptr`, for example.
   */
  template <typename T>
  class TaskResult
  {
  public:
    /**
     * Default constructor. An object of this kind is not (presently)
     * associated with a task, and so cannot be asked for that task's
     * result.
     */
    TaskResult()
      : result_is_available(false)
    {}

    /**
     * A constructor that takes a Task object and initializes the
     * current object to be the result of that task.
     */
    TaskResult(const Task<T> &task)
      : result_is_available(false)
      , task(task)
    {
      // It is conceivable that the task has already finished if we
      // come to this point, but that is not important to us here:
      // we will simply find out once someone calls get()
    }

    /**
     * Copy constructor. Because the result of a currently still
     * running task is a unique object, TaskResult objects cannot be
     * copied.
     */
    TaskResult(const TaskResult<T> &) = delete;

    /**
     * Move constructor. Following this call, the newly created object
     * represents the task's result, whereas the old object no longer
     * represents anything and is left as if default-constructed.
     */
    TaskResult(TaskResult<T> &&other) noexcept DEAL_II_CXX20_REQUIRES(
      std::is_move_constructible_v<T> &&std::is_move_assignable_v<T>);

    /**
     * Destructor. If the current object was associated with a task,
     * then the destructor will throw an error if that task is still
     * running. This is because typically, the task will still be
     * working on data that is already gone away, or will shortly go away, and
     * that means that that task will likely encounter unpredictable
     * outcomes. As an example, consider the following code:
     * @code
     *   class Complicated {
     *     public:
     *       Complicated () {
     *         ...initialize members...;
     *         hash_value = Threads::new_task(
     *           [&]() { compute has value for the current object; });
     *       }
     *
     *       ~Complicated() { ... }
     *
     *     private:
     *       ...complicated data members...
     *       Threads::TaskResult<int> hash_value;
     *       ...more data members...
     *   };
     * @endcode
     * Here, the constructor `~Complicated()` destroys the object's member
     * variables, and then calls the destructor of `hash_value`. If at the
     * time we get to the latter destructor the background task to compute
     * a hash value is still running, then this is strictly speaking not a
     * problem for the `hash_value` variable (we could simply abandon the
     * background task and its result), but it is indicative of a programming
     * bug where that background task is now working on invalid memory.
     * To catch these things, the destructor of this class will throw an error
     * if it encounters a situation where the task it is associated with is
     * still running.
     *
     * The correct solution is to wait in the surrounding class's destructor
     * for the background task to finish before tearing down data structures:
     * @code
     *       ~Complicated() {
     *         hash_value.join();
     *         ...tear down data structures...
     *       }  // compilers calls hash_value's destructor here
     * @endcode
     *
     * (One could have designed this class in such a way that the destructor
     * waits for the task to complete. But this is error prone, because C++
     * prescribes that destructors of member variables are called in reverse
     * order in which they are declared, and so if this destructor is called
     * by the compiler as part of a surrounding class's destructor, the
     * destructors of all member variables have already run by the time
     * we get to wait for the background task -- but the background task
     * is likely still working on these data members, and that presents
     * a probable bug. Instead, we force the destructor of the surrounding
     * class to explicitly wait for the task to complete before reaching
     * the end of the destructor at which the compiler inserts the calls
     * to the destructors of the class's member variables.)
     */
    ~TaskResult();

    /**
     * Copy operator. Like the copy constructor, this function is deleted
     * because TaskResult corresponds to the output of a specific task,
     * not a copy of that tasks's outcome.
     */
    TaskResult &
    operator=(const TaskResult &) = delete;

    /**
     * Move assignment operator. Following this call, the newly created object
     * represents the task's result, whereas the old object no longer
     * represents anything and is left as if default-constructed.
     */
    TaskResult &
    operator=(TaskResult &&other) noexcept DEAL_II_CXX20_REQUIRES(
      std::is_move_constructible_v<T> &&std::is_move_assignable_v<T>);

    /**
     * Copy assignment operator from a Task object. By assigning the Task
     * object to the current object, the current object is set to represent
     * the result of that task.
     *
     * When calling this operator, the current object will no longer
     * represent the result of a previous Task. If the previously associated
     * task is still running, this function throws an exception. This is
     * because assigning a new task to this object when the previous task
     * is still running is likely a bug: In most cases, the new and old
     * tasks are both operating on the state of another object; the fact
     * that the old task is still running often indicates that it is working
     * on data that has been changed underneath. To illustrate this, case,
     * consider the following code snippet:
     * @code
     *   class Complicated {
     *     public:
     *       Complicated () {
     *         ...initialize members...;
     *         hash_value = Threads::new_task(
     *           [&]() { compute has value for the current object; });
     *       }
     *
     *       void frobnicate () {
     *         ...do something with members...;
     *         hash_value = Threads::new_task(
     *           [&]() { compute has value for the current object; });
     *       }
     *
     *     private:
     *       ...complicated data members...
     *       Threads::TaskResult<int> hash_value;
     *   };
     * @endcode
     * Here, in `frobnicate()`, member variables are updated according to what
     * the `frobnicate()` operation represents, and then the `hash_value`
     * variable is updated as well, but in the background. The issue is that
     * if the background task started in the constructor is still running,
     * then that background task is working on data that is being changed
     * beneath it, likely resulting in unpredictable outcomes. The correct
     * approach would therefore be to write `frobnicate()` as follows:
     * @code
     *       void frobnicate () {
     *         hash_value.join();
     *         ...do something with members...;
     *         hash_value = Threads::new_task(
     *           [&]() { compute has value for the current object; });
     *       }
     * @endcode
     * The call to join() ensures that the program waits for the previous task
     * to finish before starting to modify the underlying data. (Calling
     * `clear()` might seem like a useful alternative in that it simply abandons
     * the first background task, but while that ensures that the first task's
     * result are not put into the current object, it still means that that task
     * is working on data that is changing as it is working, with obviously
     * unpredictable results.)
     *
     * @note As mentioned above, it is a considered a bug to assign a task to
     *   a TaskResult object that already has a task running. That means that
     *   you will get in trouble if multiple threads (or multiple tasks running
     *   concurrently) call this operator at the same time: One of these
     *   threads will set a task, and the other threads will try the same but
     *   because the background task is likely still running will encounter
     *   an error. As a consequence, you cannot easily use this operator
     *   from multiple threads. Use try_emplace_task() in that case.
     */
    void
    operator=(const Task<T> &t);

    /**
     * This function is similar to `operator=()` in that it associates a
     * task with the current object if one has not been associated so far,
     * but does not do so if a task is already assigned. For this to work,
     * the object provided as argument to this function must be a "callable"
     * (i.e., a object `creator` that can be called on a separate task via
     * its `operator()`), rather than a Task object itself.
     *
     * As a consequence, code such as the following will work:
     * @code
     *   class LazyInt
     *   {
     *     public:
     *       LazyInt () {} // no task assigned to task_result
     *
     *       int get () {
     *         task_result.try_emplace_task( []() { return 42; } );
     *         return task_result.value();
     *       }
     *
     *     private:
     *       TaskResult<int> task_result;
     *   }
     * @endcode
     * In this context, the `LazyInt::get()` function is thread-safe, i.e., it
     * can be called more than once from multiple threads. One of these threads
     * -- namely, the first one to get into `try_emplace_task()` -- will create
     * a task that calls the lambda function that returns `42` whereas all of
     * the others will simply proceed to the call to `value()` that waits for
     * the task to finish.
     *
     * On the other hand, implementing the `get()` function as
     * @code
     *       int get () {
     *         task_result = Threads::new_task( []() { return 42; } );
     *         return task_result.value();
     *       }
     * @endcode
     * would have led to the errors mentioned above because `operator=` called
     * from separate threads would have emplaced a task while another task
     * is (possibly) still running.
     *
     * @note There is nothing that prevents you from concurrently calling
     *   this function with *different* callables as arguments, i.e., with
     *   functions that create non-identical objects. That is obviously
     *   not the intent here since you can't control which callable will
     *   eventually be turned into a task.
     */
    template <typename Callable>
    void
    try_emplace_task(const Callable &creator) const
      DEAL_II_CXX20_REQUIRES((std::is_invocable_r_v<T, Callable>));

    /**
     * Instead of letting a task compute the object stored by this
     * instance of Lazy, just copy the object passed as an argument
     * and use it instead.
     *
     * As for several other member functions, this function can only be called
     * if there is not currently a running task whose result is supposed to
     * be used.
     */
    void
    emplace_object(const T &t)
      DEAL_II_CXX20_REQUIRES((std::is_copy_constructible_v<T> ||
                              std::is_copy_assignable_v<T>));

    /**
     * Instead of letting a task compute the object stored by this
     * instance of Lazy, just move the object passed as an argument
     * and use it instead.
     *
     * As for several other member functions, this function can only be called
     * if there is not currently a running task whose result is supposed to
     * be used.
     */
    void
    emplace_object(T &&t)
      DEAL_II_CXX20_REQUIRES((std::is_copy_constructible_v<T> ||
                              std::is_copy_assignable_v<T>));

    /**
     * Reset the current object to a state as if it had been
     * default-constructed. For the same reasons as outlined
     * in the documentation of the destructor and of the
     * assignment operator, this function considers it an
     * error (and throws an exception) if there is still a
     * currently running task associated with the object.
     * Because you cannot know when a task associated with an
     * object finishes, the practical realization of there not being
     * a currently still running task is that you can only call this
     * function on an object that has an associated task *after* you
     * have either called join() or asked for the return value by
     * having called value() (which internally calls join()).
     */
    void
    clear();

    /**
     * If this object is associated with a task, wait for it to finish. If
     * it isn't associated with a task, just return.
     */
    void
    join() const;

    /**
     * Return a reference to the object computed by the task. This
     * is always a `const` reference to reflect the semantics that this
     * object represents what the Threads::Task computed (which is what it
     * is, and cannot be changed later on).
     */
    const T &
    value() const;

    /**
     * Return whether the object has an associated task result object or not.
     * Only objects that are default-initialized, or for which clear() has
     * been called, or that have been moved from, will return `true` in
     * this function.
     */
    bool
    empty() const;

  private:
    /**
     * An atomic flag that allows us to test whether the task has finished
     * and the result is available.
     */
    mutable std::atomic<bool> result_is_available;

    /**
     * An object that references the task that computes the result TaskResult
     * stores. Once the task has finished, and the result has been moved out
     * of the `task` object and into `task_result`, the `task` object is
     * reset -- that is, we destroy all traces of references to the Thread::Task
     * that computed the result.
     */
    mutable std::optional<Task<T>> task;

    /**
     * The object that results from running the background task.
     */
    mutable std::optional<T> task_result;

    /**
     * A lock object that guards access to all of the `mutable` objects above.
     */
    mutable std::mutex mutex;
  };


  // ------------------------------- inline functions --------------------------

#ifndef DOXYGEN

  template <typename T>
  inline TaskResult<T>::TaskResult(TaskResult<T> &&other) noexcept
    DEAL_II_CXX20_REQUIRES(
      std::is_move_constructible_v<T> &&std::is_move_assignable_v<T>)
  {
    // First lock the other object, then move the members of the other
    // object and reset it. Note that we do not have to wait for
    // the other object's task to finish (nor should we).
    std::lock_guard<std::mutex> lock(other.mutex);

    result_is_available       = other.result_is_available.load();
    other.result_is_available = false;

    task = std::move(other.task);
    other.task.reset();

    task_result = std::move(other.task_result);
    other.task_result.reset();
  }



  template <typename T>
  inline TaskResult<T>::~TaskResult()
  {
    // Ensure that there is no currently running task. As
    // documented, we consider this an error. Since clear()
    // also checks for this error, we can just defer to that function:
    clear();
  }



  template <typename T>
  inline void
  TaskResult<T>::operator=(const Task<T> &t)
  {
    // First ensure that there is no currently running task. As
    // documented, we consider this an error.  Since clear()
    // also checks for this error, we can just defer to that function:
    clear();

    // Having established that there is no previous task still running,
    // set the current task as the one we're waiting for:
    {
      std::lock_guard<std::mutex> lock(mutex);
      task = t;
    }
  }


  template <typename T>
  inline TaskResult<T> &
  TaskResult<T>::operator=(TaskResult<T> &&other) noexcept
    DEAL_II_CXX20_REQUIRES(
      std::is_move_constructible_v<T> &&std::is_move_assignable_v<T>)
  {
    // First clear the current object before we put new content into it:
    clear();

    // Then lock the other object and move the members of the other
    // object, and finally reset it. Note that we do not have to wait for
    // the other object's task to finish (nor should we): We may simply
    // inherit the other object's task.
    std::lock_guard<std::mutex> lock(other.mutex);

    result_is_available       = other.result_is_available.load();
    other.result_is_available = false;

    task = std::move(other.task);
    other.task.reset();

    task_result = std::move(other.task_result);
    other.task_result.reset();

    return *this;
  }



  template <typename T>
  template <typename Callable>
  void
  TaskResult<T>::try_emplace_task(const Callable &creator) const
    DEAL_II_CXX20_REQUIRES((std::is_invocable_r_v<T, Callable>))
  {
    // If the result is already available, simply return.
    if (result_is_available)
      return;

    // If the result was not available above, we need to go under a lock
    // to check that perhaps it has appeared in the meantime. We again use
    // the double-checking pattern:
    {
      std::lock_guard<std::mutex> lock(mutex);
      if (result_is_available)
        return;
      else
        // If there is no result, but there is a task, some other thread has
        // emplaced it in the meantime and we can simply return
        if (task.has_value())
          return;
        else
          // If there is no task object, emplace one:
          task = Threads::new_task(creator);
    }
  }



  template <typename T>
  inline void
  TaskResult<T>::emplace_object(const T &t)
    DEAL_II_CXX20_REQUIRES((std::is_copy_constructible_v<T> ||
                            std::is_copy_assignable_v<T>))
  {
    clear();
    task_result         = t;
    result_is_available = true;
  }


  template <typename T>
  inline void
  TaskResult<T>::emplace_object(T &&t)
    DEAL_II_CXX20_REQUIRES((std::is_copy_constructible_v<T> ||
                            std::is_copy_assignable_v<T>))
  {
    clear();
    task_result         = std::move(t);
    result_is_available = true;
  }


  template <typename T>
  inline void
  TaskResult<T>::clear()
  {
    std::lock_guard<std::mutex> lock(mutex);

    if (result_is_available)
      {
        // First make clear that the result is no longer available, then
        // reset the object:
        result_is_available = false;
        task_result.reset();
      }
    else
      Assert(task.has_value() == false,
             ExcMessage("You cannot destroy a TaskResult object "
                        "while it is still waiting for its associated task "
                        "to finish. See the documentation of this class' "
                        "destructor for more information."));
  }



  template <typename T>
  inline void
  TaskResult<T>::join() const
  {
    Assert(empty() == false,
           ExcMessage("You can't join a TaskResult object that has not "
                      "been associated with a task."));

    // If we have waited before, then return immediately:
    if (result_is_available)
      return;
    else
      // If we have not waited, wait now. We need to use the double-checking
      // pattern to ensure that if two threads get to this place at the same
      // time, one returns right away while the other does the work. Note
      // that this happens under the lock, so only one thread gets to be in
      // this code block at the same time:
      {
        std::lock_guard<std::mutex> lock(mutex);
        if (result_is_available)
          return;
        else
          {
            // The object is not empty and it has not received its result yet.
            // So it must have a task object:
            Assert(task.has_value(), ExcInternalError());

            task.value().join();
            task_result = std::move(task.value().return_value());
            task.reset();

            result_is_available = true;
          }
      }
  }



  template <typename T>
  inline bool
  TaskResult<T>::empty() const
  {
    // If we have waited for a task to complete, then the object is not empty:
    if (result_is_available)
      return false;
    // Otherwise, if result_is_available has not been set, but we have a task
    // associated (i.e., the task is still running, or at least we haven't
    // waited for it to complete), then the object is also not empty:
    else if (task.has_value())
      return false;
    else
      // If when we asked above we had not joined a task, and if there was
      // no task currently associated with the object, then one of two cases
      // could have happened: either, there never was a task, and the object
      // is consequently empty. Or there was a task and somewhere between the
      // checks above and now, join() has flipped the state to
      // result_is_available==true and task.has_value()==false. We can
      // check that, but only under a lock.
      {
        std::lock_guard<std::mutex> lock(mutex);
        if (result_is_available)
          return false;
        else
          // We know from getting into the above 'else that no task was
          // associated with this object at the time. This cannot have
          // changed since then in a way that is thread-safe (i.e., by
          // way of other 'const' functions), so if the result is still
          // not available, then the object must necessarily be empty:
          return true;
      }
  }



  template <typename T>
  inline const T &
  TaskResult<T>::value() const
  {
    Assert(empty() == false,
           ExcMessage(
             "You can't ask for the result of a TaskResult object that "
             "has not been associated with a task."));

    if (!result_is_available)
      join();
    return task_result.value();
  }

#endif

} // namespace Threads


/**
 * @}
 */



DEAL_II_NAMESPACE_CLOSE
#endif

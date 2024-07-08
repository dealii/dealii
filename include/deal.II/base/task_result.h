// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2024 by the deal.II authors
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
    TaskResult(TaskResult<T> &&other) noexcept;

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
     */
    void
    operator=(const Task<T> &t);

    /**
     * Reset the current object to a state as if it had been
     * default-constructed. For the same reasons as outlined
     * in the documentation of the destructor and of the
     * assignment operator, this function considers it an
     * error (and throws an exception) if there is still a
     * currently running task associated with the object.
     */
    void
    clear();

    /**
     * If this object is associated with a task, wait for it to finish. If
     * it isn't associated with a task, just return.
     */
    void
    join();

    /**
     * Return a reference to the object computed by the task. This
     * is always a `const` reference to reflect the semantics that this
     * object represents what the Threads::Task computed (which is what it
     * is, and cannot be changed later on).
     */
    const T &
    value() const;

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

    /**
     * Wait for the task to finish, move its result into the `task_result`
     * object, and then release all information still associated with the
     * task that originally computed the result.
     */
    void
    wait_and_move_result() const;
  };


  // ------------------------------- inline functions --------------------------


  template <typename T>
  inline TaskResult<T>::TaskResult(TaskResult<T> &&other) noexcept
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
  inline void
  TaskResult<T>::clear()
  {
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
          Assert(task.has_value() == false,
                 ExcMessage("You cannot destroy a TaskResult object "
                            "while it is still waiting for its associated task "
                            "to finish. See the documentation of this class' "
                            "destructor for more information."));
      }
    std::lock_guard<std::mutex> lock(mutex);
    // First make clear that the result is no longer available:
    result_is_available = false;
    // Then abandon a previous task, should there have been one. Also abandon
    // any previously available returned object
    task.reset();
    task_result.reset();
  }



  template <typename T>
  inline void
  TaskResult<T>::join()
  {
    // If we have waited before, then return immediately:
    if (result_is_available)
      return;
    else // If we have not waited, wait now. We need to use the double-checking
         // pattern to ensure that if two threads get to this place at the same
         // time, one returns right away while the other does the work. Note
         // that this happens under the lock, so only one thread gets to be in
         // this code block at the same time:
      {
        std::lock_guard<std::mutex> lock(mutex);
        if (result_is_available)
          return;
        else
          // If there is a task, wait for it to finish. We could then move
          // the result, but it's fine to postpone that until someone actually
          // asks for the result.
          if (task.has_value())
            task.value().join();
      }
  }



  template <typename T>
  inline const T &
  TaskResult<T>::value() const
  {
    if (!result_is_available)
      wait_and_move_result();
    return task_result.value();
  }



  template <typename T>
  inline void
  TaskResult<T>::wait_and_move_result() const
  {
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
            Assert(task.has_value(),
                   ExcMessage("You cannot wait for the result of a TaskResult "
                              "object that has no task associated with it."));
            task.value().join();
            task_result = std::move(task.value().return_value());
            task.reset();

            result_is_available = true;
          }
      }
  }

} // namespace Threads

/**
 * @}
 */



DEAL_II_NAMESPACE_CLOSE
#endif

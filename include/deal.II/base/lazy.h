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

#ifndef dealii_lazy_h
#define dealii_lazy_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/task_result.h>

#include <atomic>
#include <mutex>
#include <optional>
#include <type_traits>


DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup threads
 * @{
 */

/**
 * This class is a wrapper that provides a convenient mechanism for lazy
 * initialization of the contained object on first use. The class ensures
 * that on-demand initialization of some expensive data structure happens
 * (a) exactly once in a thread-safe manner, and that (b) subsequent checks
 * in hot paths are cheap.
 *
 * Lazy<T> is closely modeled after the `std::optional` interface providing
 * a `reset()` and `value()` method, but also and extending it with two
 * methods: `ensure_initialized(creator)` which, as the name suggests,
 * ensures that the contained object is properly initialized. If the
 * `Lazy<T>` happens happens to contain no value yet, it initializes the
 * wrapped object by calling the `creator()` function object and storing
 * the return value. In addition a `value_or_initialize(creator)` function
 * is provided that, similarly, ensures that the object is properly
 * initialized and then returns a reference to the contained value.
 *
 * Example usage could look like the following, where the `FE` class stores a
 * matrix that is expensive to compute and so we do not want to do it unless it
 * is actually needed. As a consequence, rather than storing a matrix, we store
 * a `Lazy<FullMatrix<double>>` that by default is empty; whenever the matrix is
 * first requested, we create it and store it for later reuse:
 * ```
 * template<...>
 * class FE
 * {
 * public:
 *   const FullMatrix<double> & get_prolongation_matrix() const
 *   {
 *     prolongation_matrix.ensure_initialized([&](){
 *       // Some expensive operation initializing the prolongation matrix
 *       // that we only want to perform once and when necessary.
 *       });
 *     return prolongation_matrix.value();
 *   }
 *
 * private:
 *   Lazy<FullMatrix<double>> prolongation_matrix;
 * };
 * ```
 *
 * @note Conceptually, this class is not so different from
 *   [std::future](https://en.cppreference.com/w/cpp/thread/future), which
 *   can also be used to represent a possibly-not-yet-available value on which
 *   one can wait when used with the "deferred" policy of
 *   [std::async](https://en.cppreference.com/w/cpp/thread/async).
 *   In particular, the following code could be used in place
 *   of the one above:
 * ```
 * template<...>
 * class FE
 * {
 * public:
 *   FE () {
 *     prolongation_matrix = std::async(std::launch::deferred,
 *       [&](){
 *       // Some expensive operation initializing the prolongation matrix
 *       // that we only want to perform once and when necessary.
 *       });
 *   }
 *
 *   FullMatrix<double> get_prolongation_matrix() const
 *   {
 *     return prolongation_matrix.get();
 *   }
 *
 * private:
 *   std::future<FullMatrix<double>> prolongation_matrix;
 * };
 * ```
 *   The difference to what Lazy does is that for Lazy, the action must be
 *   specified in the place where we want to access the deferred computation's
 *   result. In contrast, in the scheme with `std::future` and `std::async`,
 *   the action has to be provided at the point where the `std::future`
 *   object is initialized. Both are valid approaches and, depending on
 *   context, can usefully be employed. The difference is simply in what
 *   kind of information the provided lambda function can capture: Is it
 *   the environment available at the time the constructor is run, or the
 *   environment available at the time the access function is run. The latter
 *   has the advantage that the information captured is always up to date,
 *   whereas in the scheme with `std::async`, one has to be careful not to
 *   capture information in the lambda function that could be changed by later
 *   calls to member functions but before the lambda function is finally
 *   evaluated in the getter function. (There is another difference:
 *   `std::future::get()` can only be called once, as the function returns
 *   the computed object by value and may move the object out of its internal
 *   storage. As a consequence, the call to `FE::get_prolongation_matrix()`
 *   is only valid the first time around. Lazy does not have this restriction.)
 *
 * @dealiiConceptRequires{std::is_move_constructible_v<T> &&
                          std::is_move_assignable_v<T >}
 */
template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
class Lazy
{
public:
  /**
   * Default Constructor.
   */
  Lazy();


  /**
   * Copy constructor. If the `other` object contains an initialized
   * value, then that value will be copied into the current object. If the
   * `other` object is uninitialized, then the current object will be as well.
   */
  Lazy(const Lazy &other);


  /**
   * Move constructor. If the `other` object contains an initialized
   * value, then that value will be moved into the current object, and the
   * `other` object will end up being empty (as if default initialized). If the
   * `other` object is uninitialized, then the current object will be as well.
   */
  Lazy(Lazy &&other) noexcept;


  /**
   * Copy assignment. If the `other` object contains an initialized
   * value, then that value will be copied into the current object. If the
   * `other` object is uninitialized, then the current object will be as well.
   *
   * Any content of the current object is lost in the assignment.
   */
  Lazy &
  operator=(const Lazy &other);


  /**
   * Move assignment. If the `other` object contains an initialized
   * value, then that value will be moved into the current object, and the
   * `other` object will end up being empty (as if default initialized). If the
   * `other` object is uninitialized, then the current object will be as well.
   *
   * Any content of the current object is lost in the move assignment.
   */
  Lazy &
  operator=(Lazy &&other) noexcept;


  /**
   * Reset the Lazy<T> object to an uninitialized state.
   */
  void
  reset() noexcept;


  /**
   * Initialize the wrapped object.
   *
   * If the contained object is already initialized this function simply
   * returns and does nothing.
   *
   * If, instead, the object has not yet been initialized then the @p
   * creator function object (oftentimes a lambda function) is called to
   * initialize the contained object.
   *
   * This operation is thread safe: The ensure_initialized() method
   * guarantees that the creator function object is only called once on one
   * of the calling threads and that after completion the initialization
   * result (which is stored in the std::optional) is visible on all
   * threads.
   *
   * @dealiiConceptRequires{std::is_invocable_r_v<T, Callable>}
   */
  template <typename Callable>
  void
  ensure_initialized(const Callable &creator) const
    DEAL_II_CXX20_REQUIRES((std::is_invocable_r_v<T, Callable>));


  /**
   * Returns true if the contained object has been initialized, otherwise
   * false.
   */
  bool
  has_value() const;


  /**
   * Return a const reference to the contained object.
   *
   * @pre The object has been initialized with a call to
   * ensure_initialized() or value_or_initialized().
   */
  const T &
  value() const;


  /**
   * If the underlying object is initialized the function simply returns a
   * const reference to the contained value. Otherwise, the @p creator()
   * function object is called to initialize the object first.
   *
   * This function mimics the syntax of the std::optional<T> interface and
   * is functionally equivalent to calling ensure_initialized() followed by
   * value(). It returns a `const` reference to make clear that the object
   * created by the `creator` function is what it is, and is not subject
   * to later modification unless one calls reset() and creates a new
   * object.
   *
   * @post The underlying object is initialized, meaning, has_value()
   * returns true.
   *
   * @dealiiConceptRequires{std::is_invocable_r_v<T, Callable>}
   */
  template <typename Callable>
  const T &
  value_or_initialize(const Callable &creator) const
    DEAL_II_CXX20_REQUIRES((std::is_invocable_r_v<T, Callable>));

  /**
   * Compute the memory consumption of this structure.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * A handle to the task used to create the lazily initialized object.
   */
  mutable Threads::TaskResult<T> task_result;


  /**
   * An atomic bool used for checking whether the object is initialized in
   * a thread-safe manner.
   */
  mutable std::atomic<bool> object_is_initialized;
};

/**
 * @}
 */


// ------------------------------- inline functions --------------------------

#ifndef DOXYGEN

template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline Lazy<T>::Lazy()
  : object_is_initialized(false)
{}



template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline Lazy<T>::Lazy(const Lazy &other)
{
  // If the other object has a value stored, then get it and set our
  // own object to a copy of it:
  if (other.has_value())
    {
      object_is_initialized.store(true);
      task_result.emplace_object(other.value());
    }
  else
    object_is_initialized.store(false);
}



template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline Lazy<T>::Lazy(Lazy &&other) noexcept
  : task_result(std::move(other.task_result))
{
  object_is_initialized.store(other.object_is_initialized.load());

  // Mark the other object as uninitialized.
  other.object_is_initialized.store(false);
}



template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline Lazy<T> &Lazy<T>::operator=(const Lazy &other)
{
  // If the other object has a value stored, then get it and set our
  // own object to a copy of it:
  if (other.has_value())
    {
      object_is_initialized.store(true);
      task_result.emplace_object(other.value());
    }
  else
    {
      object_is_initialized.store(false);
      task_result.clear();
    }

  return *this;
}



template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline Lazy<T> &Lazy<T>::operator=(Lazy &&other) noexcept
{
  task_result = std::move(other.task_result);
  object_is_initialized.store(other.object_is_initialized.load());

  // Mark the other object as uninitialized.
  other.object_is_initialized.store(false);

  return *this;
}



template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline void Lazy<T>::reset() noexcept
{
  object_is_initialized.store(false);
  task_result.clear();
}


template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
template <typename Callable>
inline DEAL_II_ALWAYS_INLINE
  void Lazy<T>::ensure_initialized(const Callable &creator) const
  DEAL_II_CXX20_REQUIRES((std::is_invocable_r_v<T, Callable>))
{
  //
  // Use Schmidt's double checking [1] for checking and initializing the
  // object.
  //
  // [1] https://en.wikipedia.org/wiki/Double-checked_locking
  //

  //
  // Check the object_is_initialized atomic with "acquire" semantics [1].
  //
  // This ensures that (a) all subsequent reads (of the object) are
  // ordered after this check, and that (b) all writes to the object
  // before the atomic bool was set to true with "release" semantics are
  // visible on this thread.
  //
  // [1]
  // https://en.cppreference.com/w/cpp/atomic/memory_order#Release-Acquire_ordering
  //
  if (!object_is_initialized.load(std::memory_order_acquire))
#  ifdef DEAL_II_HAVE_CXX20
    [[unlikely]]
#  endif
    {
      // Check again. If this thread won the race to the lock then we
      // would like to initialize the object. Otherwise another thread has
      // already initialized the object and flipped the object_is_initialized
      // bit. (Here, the initialization_mutex ensures consistent ordering
      // with a memory fence, so we will observe the updated bool without
      // acquire semantics.)
      //
      // Naively, we should think that we can check again for
      // object_is_initialized to be true, and if it is false just execute
      // the 'creator' function object. The problem with this approach is that
      // if we have N worker threads and spawn N+1 tasks that all want to
      // end querying the Lazy object, then N will be run right away of which
      // N-1 will block. The one remaining task will have gotten into the
      // locked section and run the 'creator' object. But if the 'creator'
      // itself spawns tasks, the scheduler may decide to first execute
      // the (N+1)st task from above which will then promptly get stuck
      // here as well -- and we're in a deadlock situation.
      //
      // The solution to the problem is to use a scheme whereby the
      // work we do in setting up 'creator' cannot block, and where we block
      // below is in a context where it's the *scheduler* that blocks to ensure
      // that it can continue to schedule tasks.
      //
      // This is what TaskResult::try_emplace_task() does:
      if (!object_is_initialized.load(std::memory_order_relaxed))
        task_result.try_emplace_task(creator);

      // At this point, either this or another thread have emplaced a
      // task and we wait for it to complete. If we do need to wait, then
      // the waiting happens in the task scheduler, which can run other
      // tasks in the meantime, ensuring progress:
      task_result.join();

      // At this point, we know that the task has completed (whether set by
      // this or another thread), and we can flip the object_is_initialized
      // boolean with "release" semantics [1].
      //
      // This ensures that the above move is visible on all threads
      // before checking the atomic bool with acquire semantics.
      // If another thread that had also been waiting in the join() call
      // above has gotten here first and set the flag to 'true', that
      // ok.
      object_is_initialized.store(true, std::memory_order_release);
    }

  Assert(has_value(),
         ExcMessage("The current object does not contain a valid object "
                    "even though we have just initialized it."));
}


template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline DEAL_II_ALWAYS_INLINE bool Lazy<T>::has_value() const
{
  //
  // In principle it would be sufficient to solely check the atomic<bool>
  // object_is_initialized because the load() is performed with "acquire"
  // semantics. But just in case let's check the object.has_value() boolean
  // as well:
  //
  return (object_is_initialized && (task_result.empty() == false));
}


template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
inline DEAL_II_ALWAYS_INLINE const T &Lazy<T>::value() const
{
  Assert(
    has_value(),
    ExcMessage(
      "value() has been called but the contained object has not been "
      "initialized. Did you forget to call 'ensure_initialized()' first?"));

  return task_result.value();
}


template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
template <typename Callable>
inline DEAL_II_ALWAYS_INLINE const T &Lazy<T>::value_or_initialize(
  const Callable &creator) const
  DEAL_II_CXX20_REQUIRES((std::is_invocable_r_v<T, Callable>))
{
  ensure_initialized(creator);
  return task_result.value();
}


template <typename T>
DEAL_II_CXX20_REQUIRES((std::is_move_constructible_v<T> &&
                        std::is_move_assignable_v<T>))
std::size_t Lazy<T>::memory_consumption() const
{
  return MemoryConsumption::memory_consumption(task_result) + //
         sizeof(*this) - sizeof(task_result);
}

#endif

DEAL_II_NAMESPACE_CLOSE
#endif

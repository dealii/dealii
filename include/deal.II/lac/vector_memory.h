// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_vector_memory_h
#define dealii_vector_memory_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/vector.h>

#include <iostream>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup VMemory */
/*@{*/

/**
 * Memory management base class for vectors. This is an abstract base class
 * used, among other places, by all iterative methods to allocate space for
 * auxiliary vectors.
 *
 * The purpose of this class is as follows: in iterative solvers and other
 * places, one needs to allocate temporary storage for vectors, for example
 * for auxiliary vectors. One could allocate and release them anew every time,
 * but this may be expensive in some situations if it has to happen very
 * frequently. A common case for this is when an iterative method is used to
 * invert a matrix in each iteration of an outer solver, such as when inverting
 * a matrix block for a Schur complement solver. (step-20 does this, for
 * example, but instead just keeps a vector around permanently for temporary
 * storage.)
 *
 * In such situations, allocating and deallocating vectors anew in each call
 * to the inner solver is expensive and leads to memory fragmentation. The
 * present class allows to avoid this by offering an interface that other
 * classes can use to allocate and deallocate vectors. Different derived
 * classes then implement different strategies to provide temporary storage
 * vectors to using classes.
 *
 * For example, the PrimitiveVectorMemory class simply allocates and
 * deallocates vectors via the operating system facilities (i.e., using
 * @p new and @p delete) each time it is asked for a vector. It is an
 * appropriate implementation to use for iterative solvers that are called
 * only once, or very infrequently.
 *
 * On the other hand, the GrowingVectorMemory class never returns memory space
 * to the operating system memory management subsystem during its lifetime; it
 * only marks them as unused and allows them to be reused next time a vector
 * is requested.
 *
 *
 * <h3> Practical use </h3>
 *
 * Classes derived from this base class return pointers to new vectors
 * via the VectorMemory::alloc() function, and re-claim the vector
 * when it is returned via VectorMemory::free(). These two functions
 * therefore play a similar role as @p new and @p delete. This
 * includes the usual drawbacks: It is simple to forget to call
 * VectorMemory::free() at the end of a function that uses this
 * facility, or to forget it in an @p if branch of the function where
 * one has an early @p return from the function. In both cases, this
 * results in a memory leak: a correct piece of code has to call
 * VectorMemory::free() for all allocated vectors at <i>all</i>
 * possible exit points. This includes places where a function is left
 * because an exception is thrown further down in the call stack and
 * not explicitly handled here.
 *
 * In other words, vectors allocated via VectorMemory::alloc() have
 * the same issue as raw pointers allocated via @p new: It is easy to
 * write code that has memory leaks. In the case of raw pointers, the
 * common solution is to use the std::unique_ptr class instead (see
 * http://en.cppreference.com/w/cpp/memory/unique_ptr). In the case of
 * the current class, the VectorMemory::Pointer class is the solution:
 * it is a class that for all practical purposes looks like a pointer,
 * but upon destruction also returns the vector back to the
 * VectorMemory object from which it got it. Since destruction of the
 * VectorMemory::Pointer class happens whenever it goes out of scope
 * (whether because the function explicitly returns, or because
 * control flow leaves it due to an exception), a memory leak cannot
 * happen: the vector the VectroMemory::Pointer object points to is
 * <i>always</i> returned.
 *
 *
 * @author Guido Kanschat, 1998-2003; Wolfgang Bangerth, 2017.
 */
template <typename VectorType = dealii::Vector<double>>
class VectorMemory : public Subscriptor
{
public:
  /**
   * Virtual destructor. This destructor is declared @p virtual to allow
   * destroying objects of derived type through pointers to this base
   * class.
   */
  virtual ~VectorMemory() override = default;

  /**
   * Return a pointer to a new vector. The number of elements or their
   * subdivision into blocks (if applicable) is unspecified and users of this
   * function should reset vectors to their proper size. The same holds for
   * the contents of vectors: they are unspecified. In other words,
   * the place that calls this function will need to resize or reinitialize
   * it appropriately.
   *
   * @warning Just like using <code>new</code> and <code>delete</code>
   *   explicitly in code invites bugs where memory is leaked (either
   *   because the corresponding <code>delete</code> is forgotten
   *   altogether, or because of exception safety issues), using the
   *   alloc() and free() functions explicitly invites writing code
   *   that accidentally leaks memory. You should consider using
   *   the VectorMemory::Pointer class instead, which provides the
   *   same kind of service that <code>std::unique</code> provides
   *   for arbitrary memory allocated on the heap.
   */
  virtual VectorType *
  alloc() = 0;

  /**
   * Return a vector and indicate that it is not going to be used any further
   * by the place that called alloc() to get a pointer to it.
   *
   * @warning Just like using <code>new</code> and <code>delete</code>
   *   explicitly in code invites bugs where memory is leaked (either
   *   because the corresponding <code>delete</code> is forgotten
   *   altogether, or because of exception safety issues), using the
   *   alloc() and free() functions explicitly invites writing code
   *   that accidentally leaks memory. You should consider using
   *   the VectorMemory::Pointer class instead, which provides the
   *   same kind of service that <code>std::unique</code> provides
   *   for arbitrary memory allocated on the heap.
   */
  virtual void
  free(const VectorType *const) = 0;

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Vector was not allocated from this memory pool.
   */
  DeclExceptionMsg(
    ExcNotAllocatedHere,
    "You are trying to deallocate a vector from a memory pool, but this "
    "vector has not actually been allocated by the same pool before.");

  //@}

  /**
   * A class that looks like a pointer for all practical purposes and that
   * upon construction time allocates a vector from a VectorMemory object
   * (or an object of a class derived from VectorMemory) that is passed
   * to the constructor of this class. The destructor then automatically
   * returns the vector's ownership to the same VectorMemory object.
   *
   * Pointers of this type are therefore safe in the sense that they
   * automatically call VectorMemory::free() when they are destroyed, whether
   * that happens at the end of a code block or because local variables are
   * destroyed during exception unwinding. These kinds of object thus relieve
   * the user from using vector management functions explicitly.
   *
   * In many senses, this class acts like <code>std::unique_ptr</code> in that
   * it is the unique owner of a chunk of memory that it frees upon destruction.
   * The main differences to <code>std::unique_ptr</code> are (i) that it
   * allocates memory from a memory pool upon construction, and (ii) that the
   * memory is not destroyed using `operator delete` but returned to the
   * VectorMemory pool.
   *
   * @author Guido Kanschat, 2009; Wolfgang Bangerth, 2017.
   */
  class Pointer
    : public std::unique_ptr<VectorType, std::function<void(VectorType *)>>
  {
  public:
    /**
     * Default constructor. This constructor corresponds to a @p nullptr
     * object that does not own a vector. It can, however, later be
     * assigned another Pointer object via move assignment in which case
     * it will steal the vector owned by the other object
     * (as @p std::unique_ptr does).
     */
    Pointer() = default;

    /**
     * Move constructor: this creates a new Pointer by stealing the internal
     * data owned by @p p.
     */
    Pointer(Pointer &&p) noexcept = default;

    /**
     * Move operator: this releases the vector owned by the current Pointer
     * and then steals the internal data owned by @p p.
     */
    Pointer &
    operator=(Pointer &&p) noexcept = default;

    /**
     * Constructor. This constructor automatically allocates a vector from
     * the given vector memory object @p mem.
     */
    Pointer(VectorMemory<VectorType> &mem);

    /**
     * Destructor, automatically releasing the vector from the memory pool.
     */
    ~Pointer() = default;
  };
};



/**
 * Simple memory management. See the documentation of the base class for a
 * description of its purpose.
 *
 * This class allocates and deletes vectors as needed from the global heap,
 * i.e. performs no specially adapted actions for memory management.
 */
template <typename VectorType = dealii::Vector<double>>
class PrimitiveVectorMemory : public VectorMemory<VectorType>
{
public:
  /**
   * Return a pointer to a new vector. The number of elements or their
   * subdivision into blocks (if applicable) is unspecified and users of this
   * function should reset vectors to their proper size. The same holds for
   * the contents of vectors: they are unspecified. In other words,
   * the place that calls this function will need to resize or reinitialize
   * it appropriately.
   *
   * For the present class, calling this function will allocate a new vector
   * on the heap and returning a pointer to it. Later calling free() then
   * returns the memory to the global heap managed by the operating system.
   *
   * @warning Just like using <code>new</code> and <code>delete</code>
   *   explicitly in code invites bugs where memory is leaked (either
   *   because the corresponding <code>delete</code> is forgotten
   *   altogether, or because of exception safety issues), using the
   *   alloc() and free() functions explicitly invites writing code
   *   that accidentally leaks memory. You should consider using
   *   the VectorMemory::Pointer class instead, which provides the
   *   same kind of service that <code>std::unique</code> provides
   *   for arbitrary memory allocated on the heap.
   */
  virtual VectorType *
  alloc() override;

  /**
   * Return a vector and indicate that it is not going to be used any further
   * by the instance that called alloc() to get a pointer to it.
   *
   * For the present class, this means that the vector is returned to the
   * global heap.
   *
   * @warning Just like using <code>new</code> and <code>delete</code>
   *   explicitly in code invites bugs where memory is leaked (either
   *   because the corresponding <code>delete</code> is forgotten
   *   altogether, or because of exception safety issues), using the
   *   alloc() and free() functions explicitly invites writing code
   *   that accidentally leaks memory. You should consider using
   *   the VectorMemory::Pointer class instead, which provides the
   *   same kind of service that <code>std::unique</code> provides
   *   for arbitrary memory allocated on the heap.
   */
  virtual void
  free(const VectorType *const v) override;
};



/**
 * A pool based memory management class. See the documentation of the base
 * class for a description of its purpose.
 *
 * Each time a vector is requested from this class, it checks if it has one
 * available and returns its address, or allocates a new one on the heap. If a
 * vector is returned from its user, through the GrowingVectorMemory::free()
 * member function, it doesn't return the allocated memory to the operating
 * system memory subsystem, but keeps it around unused for later use if
 * GrowingVectorMemory::alloc() is called again. The
 * class therefore avoid the overhead of repeatedly allocating memory on the
 * heap if temporary vectors are required and released frequently; on the
 * other hand, it doesn't release once-allocated memory at the earliest
 * possible time and may therefore lead to an increased overall memory
 * consumption.
 *
 * All GrowingVectorMemory objects of the same vector type use the same memory
 * pool. (In other words: The pool of vectors from which this class draws is
 * <i>global</i>, rather than a regular member variable of the current class
 * that is destroyed at the time that the surrounding GrowingVectorMemory
 * object is destroyed.) Therefore, functions can create such a
 * GrowingVectorMemory object whenever needed without the performance penalty
 * of creating a new memory pool every time. A drawback of this policy is that
 * vectors once allocated are only released at the end of the program run.
 *
 * @author Guido Kanschat, 1999, 2007; Wolfgang Bangerth, 2017.
 */
template <typename VectorType = dealii::Vector<double>>
class GrowingVectorMemory : public VectorMemory<VectorType>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Constructor.  The argument allows to preallocate a certain number of
   * vectors. The default is not to do this.
   */
  GrowingVectorMemory(const size_type initial_size   = 0,
                      const bool      log_statistics = false);

  /**
   * Destructor. The destructor also checks that all vectors that have been
   * allocated through the current object have all been released again.
   * However, as discussed in the class documentation, this does not imply
   * that their memory is returned to the operating system.
   */
  virtual ~GrowingVectorMemory() override;

  /**
   * Return a pointer to a new vector. The number of elements or their
   * subdivision into blocks (if applicable) is unspecified and users of this
   * function should reset vectors to their proper size. The same holds for
   * the contents of vectors: they are unspecified. In other words,
   * the place that calls this function will need to resize or reinitialize
   * it appropriately.
   *
   * @warning Just like using <code>new</code> and <code>delete</code>
   *   explicitly in code invites bugs where memory is leaked (either
   *   because the corresponding <code>delete</code> is forgotten
   *   altogether, or because of exception safety issues), using the
   *   alloc() and free() functions explicitly invites writing code
   *   that accidentally leaks memory. You should consider using
   *   the VectorMemory::Pointer class instead, which provides the
   *   same kind of service that <code>std::unique</code> provides
   *   for arbitrary memory allocated on the heap.
   */
  virtual VectorType *
  alloc() override;

  /**
   * Return a vector and indicate that it is not going to be used any further
   * by the instance that called alloc() to get a pointer to it.
   *
   * For the present class, this means retaining the vector for later reuse by
   * the alloc() method.
   *
   * @warning Just like using <code>new</code> and <code>delete</code>
   *   explicitly in code invites bugs where memory is leaked (either
   *   because the corresponding <code>delete</code> is forgotten
   *   altogether, or because of exception safety issues), using the
   *   alloc() and free() functions explicitly invites writing code
   *   that accidentally leaks memory. You should consider using
   *   the VectorMemory::Pointer class instead, which provides the
   *   same kind of service that <code>std::unique</code> provides
   *   for arbitrary memory allocated on the heap.
   */
  virtual void
  free(const VectorType *const) override;

  /**
   * Release all vectors that are not currently in use.
   */
  static void
  release_unused_memory();

  /**
   * Memory consumed by this class and all currently allocated vectors.
   */
  virtual std::size_t
  memory_consumption() const;

private:
  /**
   * A type that describes this entries of an array that represents
   * the vectors stored by this object. The first component of the pair
   * is be a flag telling whether the vector is used, the second
   * a pointer to the vector itself.
   */
  using entry_type = std::pair<bool, std::unique_ptr<VectorType>>;

  /**
   * The class providing the actual storage for the memory pool.
   *
   * This is where the actual storage for GrowingVectorMemory is provided.
   * Only one of these pools is used for each vector type, thus allocating all
   * vectors from the same storage.
   *
   * @author Guido Kanschat, 2007, Wolfgang Bangerth 2017.
   */
  struct Pool
  {
    /**
     * Standard constructor creating an empty pool
     */
    Pool();

    /**
     * Destructor.
     */
    ~Pool();

    /**
     * Create data vector; does nothing after first initialization
     */
    void
    initialize(const size_type size);

    /**
     * Pointer to the storage object
     */
    std::vector<entry_type> *data;
  };

  /**
   * Return an array of allocated vectors.
   */
  static Pool &
  get_pool();

  /**
   * Overall number of allocations. Only used for bookkeeping and to generate
   * output at the end of an object's lifetime.
   */
  size_type total_alloc;

  /**
   * Number of vectors currently allocated in this object; used for detecting
   * memory leaks.
   */
  size_type current_alloc;

  /**
   * A flag controlling the logging of statistics by the destructor.
   */
  bool log_statistics;

  /**
   * Mutex to synchronize access to internal data of this object from multiple
   * threads.
   */
  static Threads::Mutex mutex;
};



namespace internal
{
  namespace GrowingVectorMemoryImplementation
  {
    void
    release_all_unused_memory();
  }
} // namespace internal

/*@}*/

#ifndef DOXYGEN
/* --------------------- inline functions ---------------------- */


template <typename VectorType>
inline VectorMemory<VectorType>::Pointer::Pointer(VectorMemory<VectorType> &mem)
  : std::unique_ptr<VectorType, std::function<void(VectorType *)>>(
      mem.alloc(),
      [&mem](VectorType *v) { mem.free(v); })
{}



template <typename VectorType>
VectorType *
PrimitiveVectorMemory<VectorType>::alloc()
{
  return new VectorType();
}



template <typename VectorType>
void
PrimitiveVectorMemory<VectorType>::free(const VectorType *const v)
{
  delete v;
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__petsc_vector_h
#define __deal2__petsc_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/subscriptor.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#  include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN


namespace PETScWrappers
{
  /*! @addtogroup PETScWrappers
   *@{
   */

  /**
   * Implementation of a sequential vector class based on PETSC. All the
   * functionality is actually in the base class, except for the calls to
   * generate a sequential vector. This is possible since PETSc only works on an
   * abstract vector type and internally distributes to functions that do the
   * actual work depending on the actual vector type (much like using virtual
   * functions). Only the functions creating a vector of specific type differ,
   * and are implemented in this particular class.
   *
   * @ingroup Vectors
   * @author Wolfgang Bangerth, 2004
   */
  class Vector : public VectorBase
  {
  public:
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * A variable that indicates whether this vector
     * supports distributed data storage. If true, then
     * this vector also needs an appropriate compress()
     * function that allows communicating recent set or
     * add operations to individual elements to be communicated
     * to other processors.
     *
     * For the current class, the variable equals
     * false, since it does not support parallel data storage.
     * If you do need parallel data storage, use
     * PETScWrappers::MPI::Vector.
     */
    static const bool supports_distributed_data = false;

    /**
     * Default constructor. Initialize the
     * vector as empty.
     */
    Vector ();

    /**
     * Constructor. Set dimension to
     * @p n and initialize all
     * elements with zero.
     *
     * The constructor is made explicit to
     * avoid accidents like this:
     * <tt>v=0;</tt>. Presumably, the user wants
     * to set every element of the vector to
     * zero, but instead, what happens is
     * this call: <tt>v=Vector@<number@>(0);</tt>,
     * i.e. the vector is replaced by one of
     * length zero.
     */
    explicit Vector (const size_type n);

    /**
     * Copy-constructor from deal.II
     * vectors. Sets the dimension to that
     * of the given vector, and copies all
     * elements.
     */
    template <typename Number>
    explicit Vector (const dealii::Vector<Number> &v);

    /**
     * Construct it from an existing PETSc
     * Vector of type Vec. Note: this does
     * not copy the contents and just keeps
     * a pointer. You need to make sure the
     * vector is not used twice at the same
     * time or destroyed while in use. This
     * class does not destroy the PETSc
     * object. Handle with care!
     */
    explicit Vector (const Vec &v);

    /**
     * Copy-constructor the values from a
     * PETSc wrapper vector class.
     */
    Vector (const Vector &v);

    /**
     * Copy-constructor: copy the values
     * from a PETSc wrapper parallel vector
     * class.
     *
     * Note that due to the communication
     * model of MPI, @em all processes have
     * to actually perform this operation,
     * even if they do not use the
     * result. It is not sufficient if only
     * one processor tries to copy the
     * elements from the other processors
     * over to its own process space.
     */
    explicit Vector (const MPI::Vector &v);

    /**
     * Copy the given vector. Resize the
     * present vector if necessary.
     */
    Vector &operator = (const Vector &v);

    /**
     * Copy all the elements of the
     * parallel vector @p v into this
     * local vector. Note that due to the
     * communication model of MPI, @em all
     * processes have to actually perform
     * this operation, even if they do not
     * use the result. It is not sufficient
     * if only one processor tries to copy
     * the elements from the other
     * processors over to its own process
     * space.
     */
    Vector &operator = (const MPI::Vector &v);

    /**
     * Set all components of the vector to
     * the given number @p s. Simply pass
     * this down to the base class, but we
     * still need to declare this function
     * to make the example given in the
     * discussion about making the
     * constructor explicit work.
     *
     * Since the semantics of assigning a
     * scalar to a vector are not
     * immediately clear, this operator
     * should really only be used if you
     * want to set the entire vector to
     * zero. This allows the intuitive
     * notation <tt>v=0</tt>. Assigning
     * other values is deprecated and may
     * be disallowed in the future.
     */
    Vector &operator = (const PetscScalar s);

    /**
     * Copy the values of a deal.II vector
     * (as opposed to those of the PETSc
     * vector wrapper class) into this
     * object.
     */
    template <typename number>
    Vector &operator = (const dealii::Vector<number> &v);

    /**
     * Change the dimension of the vector
     * to @p N. It is unspecified how
     * resizing the vector affects the
     * memory allocation of this object;
     * i.e., it is not guaranteed that
     * resizing it to a smaller size
     * actually also reduces memory
     * consumption, or if for efficiency
     * the same amount of memory is used
     * for less data.
     *
     * If @p fast is false, the vector is
     * filled by zeros. Otherwise, the
     * elements are left an unspecified
     * state.
     */
    void reinit (const size_type N,
                 const bool      fast = false);

    /**
     * Change the dimension to that of the
     * vector @p v. The same applies as
     * for the other reinit() function.
     *
     * The elements of @p v are not
     * copied, i.e.  this function is the
     * same as calling <tt>reinit (v.size(),
     * fast)</tt>.
     */
    void reinit (const Vector &v,
                 const bool    fast = false);

  protected:
    /**
     * Create a vector of length @p n. For
     * this class, we create a sequential
     * vector. @p n denotes the total
     * size of the vector to be
     * created.
     */
    void create_vector (const size_type n);
  };

  /*@}*/

// ------------------ template and inline functions -------------


  /**
   * Global function @p swap which overloads the default implementation
   * of the C++ standard library which uses a temporary object. The
   * function simply exchanges the data of the two vectors.
   *
   * @relates PETScWrappers::Vector
   * @author Wolfgang Bangerth, 2004
   */
  inline
  void swap (Vector &u, Vector &v)
  {
    u.swap (v);
  }


#ifndef DOXYGEN

  template <typename number>
  Vector::Vector (const dealii::Vector<number> &v)
  {
    Vector::create_vector (v.size());

    *this = v;
  }



  inline
  Vector::Vector (const Vec &v)
    :
    VectorBase(v)
  {}


  inline
  Vector &
  Vector::operator = (const PetscScalar s)
  {
    VectorBase::operator = (s);

    return *this;
  }


  inline
  Vector &
  Vector::operator = (const Vector &v)
  {
    // if the vectors have different sizes,
    // then first resize the present one
    if (size() != v.size())
      reinit (v.size(), true);

    const int ierr = VecCopy (v.vector, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  inline
  Vector &
  Vector::operator = (const MPI::Vector &v)
  {
    int ierr;
    if (attained_ownership)
      {
        // the petsc function we call wants to
        // generate the vector itself, so destroy
        // the old one first
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
        ierr = VecDestroy (vector);
#else
        ierr = VecDestroy (&vector);
#endif
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }

    attained_ownership = true;

    // then do the gather
    // operation. <rant>petsc has changed its
    // interface again, and replaced a single
    // function call by several calls that
    // are hard to understand. gets me all
    // annoyed at their development
    // model</rant>
#if DEAL_II_PETSC_VERSION_LT(2,2,0)
    ierr = VecConvertMPIToSeqAll (static_cast<const Vec &>(v),
                                  &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#else

    VecScatter ctx;

    ierr = VecScatterCreateToAll (static_cast<const Vec &>(v), &ctx, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if ((PETSC_VERSION_MAJOR == 2) && \
     ((PETSC_VERSION_MINOR < 3) || \
      ((PETSC_VERSION_MINOR == 3) &&            \
       (PETSC_VERSION_SUBMINOR < 3))))
    ierr = VecScatterBegin (static_cast<const Vec &>(v), vector,
                            INSERT_VALUES, SCATTER_FORWARD, ctx);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = VecScatterEnd (static_cast<const Vec &>(v), vector,
                          INSERT_VALUES, SCATTER_FORWARD, ctx);

#else

    ierr = VecScatterBegin (ctx,static_cast<const Vec &>(v), vector,
                            INSERT_VALUES, SCATTER_FORWARD);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = VecScatterEnd (ctx, static_cast<const Vec &>(v), vector,
                          INSERT_VALUES, SCATTER_FORWARD);

#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    ierr = VecScatterDestroy (ctx);
#else
    ierr = VecScatterDestroy (&ctx);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif

    return *this;
  }



  template <typename number>
  inline
  Vector &
  Vector::operator = (const dealii::Vector<number> &v)
  {
    reinit (v.size(), true);
    // the following isn't necessarily fast,
    // but this is due to the fact that PETSc
    // doesn't offer an inlined access
    // operator.
    //
    // if someone wants to contribute some
    // code: to make this code faster, one
    // could either first convert all values
    // to PetscScalar, and then set them all
    // at once using VecSetValues. This has
    // the drawback that it could take quite
    // some memory, if the vector is large,
    // and it would in addition allocate
    // memory on the heap, which is
    // expensive. an alternative would be to
    // split the vector into chunks of, say,
    // 128 elements, convert a chunk at a
    // time and set it in the output vector
    // using VecSetValues. since 128 elements
    // is small enough, this could easily be
    // allocated on the stack (as a local
    // variable) which would make the whole
    // thing much more efficient.
    //
    // a second way to make things faster is
    // for the special case that
    // number==PetscScalar. we could then
    // declare a specialization of this
    // template, and omit the conversion. the
    // problem with this is that the best we
    // can do is to use VecSetValues, but
    // this isn't very efficient either: it
    // wants to see an array of indices,
    // which in this case a) again takes up a
    // whole lot of memory on the heap, and
    // b) is totally dumb since its content
    // would simply be the sequence
    // 0,1,2,3,...,n. the best of all worlds
    // would probably be a function in Petsc
    // that would take a pointer to an array
    // of PetscScalar values and simply copy
    // n elements verbatim into the vector...
    for (size_type i=0; i<v.size(); ++i)
      (*this)(i) = v(i);

    compress (::dealii::VectorOperation::insert);

    return *this;
  }
#endif // DOXYGEN

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

/*----------------------------   petsc_vector.h     ---------------------------*/

#endif
/*----------------------------   petsc_vector.h     ---------------------------*/

// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
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

#ifndef dealii__dof_vector_h
#define dealii__dof_vector_h

#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN
/**
 * An object storing a vector describing a finite element function
 * together with its DoFHandler and the ConstraintMatrix associated
 * with the handler.
 *
 * This class encapsulates a data vector with the finite element space
 * information associated with it, such that every operation which
 * needs function values of a finite element function, namely every
 * operation which is not purely algebraic can associate the right
 * finite element space with the data vector. Thus, functions
 * implementing such operations need only a single argument of type
 * DoFVector instead of currently at least three arguments with
 * separate information on data and finite element space (the latter
 * being contributed by DoFHandler and ConstraintMatrix).  To begin,
 * it is implemented in a minimally intrusive way, giving access to
 * both structures.
 *
 * @todo Ensure that both objects are synchronized, that is, the
 * DoFHandler is not changed without the vector being resized. Answer
 * questions on unusual data structures that cannot be handled by
 * sync() now first.
 */
template <class DH, class VECTOR=Vector<double> >
class DoFVector :  public Subscriptor
{
public:
  /// @name Constructors and initialization
  ///@{
  /**
   * Construct an object with a DoFHandler and borrow a vector from
   * GrowingVectorMemory for its data. This vector is then owned by
   * this object and released when the object itself expires.
   *
   * The life span of both arguments of this constructor must exceed
   * the life span of the constructed object, otherwise, Subscriptor
   * will throw an exception complaining that a used object is being
   * deleted.
   */
  DoFVector (const DH &dh, const ConstraintMatrix &constraints);

  /**
   * Construct an object with a DoFHandler, but without constraints and
   * borrow a vector from GrowingVectorMemory for its data. This
   * vector is then owned by this object and released when the object
   * itself expires.
   *
   * @note This constructor uses an empty ConstraintMatrix. The
   * created object thus can output vectors in distributed format (see
   * @ref constraints on compressed and distributed vectors). But
   * writing into vectors will most likely fail, if there are
   * constraints associated with the dof handler, since this class
   * cannot know about them. It cannot even test for them, thus there
   * cannot be a warning. Exceptions are discontinuous Galerkin
   * methods which do not require constraints for hanging nodes or at
   * the boundary.
   *
   * The life span of the argument of this constructor must exceed
   * the life span of the constructed object, otherwise, Subscriptor
   * will throw an exception complaining that a used object is being
   * deleted.
   */
  explicit DoFVector (const DH &dh);

  /**
   * Construct an object with a DoFHandler and use the vector provided
   * as data vector without copying it. This vector is not owned by the
   * object and the constructed object will be writeable.
   *
   * The life span of all objects usesd as arguments must exceed the
   * life span of the constructed object, otherwise, Subscriptor will
   * throw an exception complaining that a used object is being
   * deleted.
   */
  DoFVector (const DH &dh, const ConstraintMatrix &constraints, VECTOR &v);

  /**
   * Construct an object with a DoFHandler and use the vector provided
   * as data vector without copying it. This vector is not owned by the
   * object and read-access only.
   *
   * @note An object created this way is effectively constant. In
   * order to avoid obscure access errors to the data vector, it
   * should also be declared `const`.
   *
   * The life span of all objects usesd as arguments must exceed the
   * life span of the constructed object, otherwise, Subscriptor will
   * throw an exception complaining that a used object is being
   * deleted.
   */
  DoFVector (const DH &dh, const ConstraintMatrix &constraints, const VECTOR &v);

  /**
   * Construct an object with a DoFHandler but without constraints and
  * use the vector provided as data vector without copying it. This
  * vector is not owned by the object, but the object is writeable.
  *
  * The life span of the arguments of this constructor must exceed
  * the life span of the constructed object, otherwise, Subscriptor
  * will throw an exception complaining that a used object is being
  * deleted.
   */
  DoFVector (const DH &dh, VECTOR &v);


  /**
   * Construct an object with a DoFHandler but without constraints and
   * use the vector provided as data vector without copying it. This
   * vector is not owned by the object and read-access only.
   *
   * @note An object created this way is effectively constant. In
   * order to avoid access errors to the data vector, it should also
   * be declared `const`. Also, vectors are expected to be in
   * distributed format (see @ref constraints on compressed and
   * distributed vectors), since the information on constraints is
   * missing.
   *
   * The life span of the arguments of this constructor must exceed
   * the life span of the constructed object, otherwise, Subscriptor
   * will throw an exception complaining that a used object is being
   * deleted.
   */
  DoFVector (const DH &dh, const VECTOR &v);

  /**
   * Destructor, deleting own data.
   */
  virtual ~DoFVector ();

  /**
   * Assignment operator. Copies the information from the other
   * object and resizes the data vector accordingly. Will only be
   * possible with an object that owns its data vector.
   */
  DoFVector<DH, VECTOR> &operator= (const DoFVector<DH, VECTOR> &);

  /**
   * Modify the vector dimensions to match the current state of the
   * dof handler and set it to zero. Will only be possible with an
   * object that owns its data vector.
   *
   * This function calls Vector::reinit() with the total number of
   * degrees of freedom and BlockVector::reinit() with the global
   * block structure of the BlockInfo object of the dof handler. If
   * the vector is an MGLevelObject, it is resized to the number of
   * levels of this object and then, the appropriate initialization is
   * run on each level.
   *
   * @todo Implementation for TrilinosWrappers::Vector and
   * PETScWrappers::Vector as well as their block versions. Check
   * whether the implementation for parallel::distributed::BlockVector
   * is as intended.
   */
  void sync ();

  ///@}

  /// @name Access to members
  ///@{
  /**
   * Access to the DoFHandler
   */
  const DH &get_dof_handler () const;

  /**
   * Access ot the constraint matrix.
   */
  const ConstraintMatrix &get_constraint_matrix () const;

  /**
   * Read-only access to the data vector
   */
  const VECTOR &get_data () const;

  /**
   * Read-write access to the data vector.
   *
   * Read-write access is only possible, if the vector is owned by the
   * object.
   */
  VECTOR &get_data ();
  ///@}

  /// @name Information on the object
  ///@{
  /**
   * Tell if the vector is mutable or not. This depends on the way
   * this object is constructed, whether the vector reference given
   * to it was a const reference or not.
   */
  bool is_const ();

  /**
   * The object is not mutable, so this function invariably returns `true`.
   */
  bool is_const () const;

  /**
   * Tell whether this object owns the vector stored in it. The
   * result is `false` if a vector was supplied to the constructor.
   */
  bool is_owner () const;

  ///@}
  /**
   * Exception thrown when a function tries to modify data of an
   * object that does not own it. Check out the different
   * constructors to fix this problem.
   */
  DeclExceptionMsg(ExcNotOwner,
                   "You are trying to modify the data vector of a DoFVector"
                   "object that does not own its data vector. Only DoFVector"
                   " objects that have not received the data vector in their"
                   " constructor can perform this action.");
private:

  /**
   * Disable copy constructor.
   */
  DoFVector (const DoFVector<DH, VECTOR> &);

  /**
   * The pool for allocating the vector.
   */
  GrowingVectorMemory<VECTOR> mem;

  /**
   * Empty constrainty matrix for constructors without one.
   */
  static const ConstraintMatrix no_constraints;

  /**
   * Pointer to the dof handler.
   */
  SmartPointer<const DH, DoFVector<DH, VECTOR> > dh;

  /**
   * Pointer to the dof handler.
   */
  SmartPointer<const ConstraintMatrix, DoFVector<DH, VECTOR> > constraints;

  /**
   * The actual data vector, if owned by this object. This pointer
   * being nonzero implies #other_data is zero.
   */
  VECTOR *my_data;

  /**
   * The data vector if referencing a vector not owned by this
   * object. This pointer being nonzero implies #my_data is zero. If
   * this pointer and #my_data are zero, it is implied that this object
   * was created with a reference to a const vector.
   */
  SmartPointer<VECTOR, DoFVector<DH, VECTOR> > other_data;

  /**
   * The data vector if referencing a vector from a const object. Is
   * equal to either #my_data or #other_data if the object is
   * writeable. If this object was constructed with a const vector
   * reference, both other pointers are zero, but this one is still
   * accessible.
   */
  SmartPointer<const VECTOR, DoFVector<DH, VECTOR> > const_data;
};


template <class DH, class VECTOR>
inline
const DH &
DoFVector<DH, VECTOR>::get_dof_handler () const
{
  return *dh;
}


template <class DH, class VECTOR>
inline
const ConstraintMatrix &
DoFVector<DH, VECTOR>::get_constraint_matrix () const
{
  return *constraints;
}


template <class DH, class VECTOR>
inline
const VECTOR &
DoFVector<DH, VECTOR>::get_data () const
{
  return *const_data;
}


template <class DH, class VECTOR>
inline
VECTOR &
DoFVector<DH, VECTOR>::get_data ()
{
  if (my_data)
    return *my_data;
  Assert (other_data != 0, ExcNotInitialized());
  return *other_data;
}


template <class DH, class VECTOR>
inline
bool
DoFVector<DH, VECTOR>::is_const ()
{
  return (my_data == 0 && other_data == 0) ? true : false;
}


template <class DH, class VECTOR>
inline
bool
DoFVector<DH, VECTOR>::is_const () const
{
  return true;
}


template <class DH, class VECTOR>
inline
bool
DoFVector<DH, VECTOR>::is_owner () const
{
  return my_data != 0;
}


DEAL_II_NAMESPACE_CLOSE

#endif

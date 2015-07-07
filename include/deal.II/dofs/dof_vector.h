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
 * Data vectors for finite element discretizations do not make sense
 * without their DoFHandler. Nevertheless, deal.II has assumed for
 * decades, that there is a single DoFHandler in a program and that
 * all vectors refer to it. This class encapsulates both, to begin
 * with in a minimally intrusive way, giving access to both
 * structures.
 *
 * @todo Ensure that both objects are synchronized, that is, the
 * DoFHandler is not changes without the vector being resized.
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
   */
  DoFVector (const DH &dh, const ConstraintMatrix &constraints);

  /**
   * Construct an object with a DoFHandler but without constraints and
   * borrow a vector from GrowingVectorMemory for its data. This
   * vector is then owned by this object and released when the object
   * itself expires.
   *
   * @note This constructor uses the ConstraintMatrix
   * #no_constraints. The created object thus can output vectors in
   * distributed format, but writing into vectors will most likely
   * fail without using a constraint matrix from somewhere else.
   */
  DoFVector (const DH &dh);

  /**
   * Construct an object with a DoFHandler and borrow a vector from
   * GrowingVectorMemory for its data. This vector is not owned by the
   * object and read-access only.
   *
   * @note An object created this way is effectively constant. In
   * order to avoid obscure access errors to the data vector, it
   * should also be declared `const`.
   */
  DoFVector (const DH &dh, const ConstraintMatrix &constraints, const VECTOR &v);

  /**
   * Construct an object with a DoFHandler but without constraints and
   * borrow a vector from GrowingVectorMemory for its data. This
   * vector is not owned by the object and read-access only.
   *
   * @note An object created this way is effectively constant. In
   * order to avoid access errors to the data vector, it should also
   * be declared `const`. Also, vectors are expected to be in
   * distributed format, since the information on constraints is
   * missing.
   */
  DoFVector (const DH &dh, const VECTOR &v);

  /**
   * Destructor, deleting own data.
   */
  ~DoFVector ();

  /**
   * Modify the vector dimensions to match the current state of the
   * dof handler.
   */
  void sync ();

  ///@}

  /// @name Access to members
  ///@{
  /**
   * Access to the DoFHandler
   */
  const DH &dof_handler () const;

  /**
   * Read-only access to the data vector
   */
  const VECTOR &data () const;

  /**
   * Read-write access to the data vector.
   *
   * Read-write access is only possible, if the vector is owned by the
   * object.
   */
  VECTOR &data ();
  ///@}

private:
  /**
   * The pool for allocating the vector.
   */
  GrowingVectorMemory<VECTOR> mem;

  /**
   * Empty constrainty matrix for constructors without one.
   */
  ConstraintMatrix no_constraints;

  /**
   * Pointer to the dof handler.
   */
  SmartPointer<const DH, DoFVector<DH, VECTOR> > dh;

  /**
   * Pointer to the dof handler.
   */
  SmartPointer<const ConstraintMatrix, DoFVector<DH, VECTOR> > constraints;

  /**
   * The actual data vector, if owned by this object.
   */
  VECTOR *my_data;

  /**
   * The data vector if referencing a vector not owned by this object.
   */
  SmartPointer<const VECTOR, DoFVector<DH, VECTOR> > other_data;
};


template <class DH, class VECTOR>
inline
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const ConstraintMatrix &co)
  : dh(&dh),
    constraints(&co),
    my_data(0),
    other_data(0)
{
  my_data = mem.alloc();
}


template <class DH, class VECTOR>
inline
DoFVector<DH, VECTOR>::DoFVector (const DH &dh)
  : dh(&dh),
    constraints(&no_constraints),
    my_data(0),
    other_data(0)
{
  my_data = mem.alloc();
}


template <class DH, class VECTOR>
inline
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const ConstraintMatrix &co, const VECTOR &v)
  : dh(&dh),
    constraints(&co),
    my_data(0),
    other_data(&v)
{}


template <class DH, class VECTOR>
inline
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const VECTOR &v)
  : dh(&dh),
    constraints(&no_constraints),
    my_data(0),
    other_data(&v)
{}


template <class DH, class VECTOR>
inline
DoFVector<DH, VECTOR>::~DoFVector ()
{
  if (my_data != 0)
    mem.free(my_data);
}


template <class DH, class VECTOR>
inline
const DH &
DoFVector<DH, VECTOR>::dof_handler () const
{
  return *dh;
}


template <class DH, class VECTOR>
inline
const VECTOR &
DoFVector<DH, VECTOR>::data () const
{
  if (my_data != 0)
    return *my_data;
  Assert(other_data != 0, ExcNotInitialized());
  return *other_data;
}


template <class DH, class VECTOR>
inline
VECTOR &
DoFVector<DH, VECTOR>::data ()
{
  Assert(my_data != 0, ExcNotInitialized());
  return *my_data;
}


DEAL_II_NAMESPACE_CLOSE

#endif

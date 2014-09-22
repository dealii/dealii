// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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

#ifndef __deal2__identity_matrix_h
#define __deal2__identity_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Matrix1
 *@{
 */


/**
 * Implementation of a simple class representing the identity matrix
 * of a given size, i.e. a matrix with entries
 * $A_{ij}=\delta_{ij}$. While it has the most important ingredients
 * of a matrix, in particular that one can ask for its size and
 * perform matrix-vector products with it, a matrix of this type is
 * really only useful in two contexts: preconditioning and
 * initializing other matrices.
 *

 * <h4>Initialization</h4>
 *
 * The main usefulness of this class lies in its ability to initialize
 * other matrix, like this:
 @code
   FullMatrix<double> identity (IdentityMatrix(10));
 @endcode
 *
 * This creates a $10\times 10$ matrix with ones on the diagonal and
 * zeros everywhere else. Most matrix types, in particular FullMatrix
 * and SparseMatrix, have conversion constructors and assignment
 * operators for IdentityMatrix, and can therefore be filled rather
 * easily with identity matrices.
 *
 *
 * <h4>Preconditioning</h4>
 *
 * No preconditioning at all is equivalent to preconditioning with
 * preconditioning with the identity matrix. deal.II has a specialized
 * class for this purpose, PreconditionIdentity, than can be used in a
 * context as shown in the documentation of that class. The present
 * class can be used in much the same way, although without any
 * additional benefit:
 @code
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);
  cg.solve (system_matrix, solution, system_rhs,
            IdentityMatrix(solution.size()));
 @endcode
 *
 *
 * @author Wolfgang Bangerth, 2006
 */
class IdentityMatrix
{
public:
  /**
   * Declate type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Default constructor. Creates a
   * zero-sized matrix that should
   * be resized later on using the
   * reinit() function.
   */
  IdentityMatrix ();

  /**
   * Constructor. Creates a
   * identity matrix of size #n.
   */
  IdentityMatrix (const size_type n);

  /**
   * Resize the matrix to be of
   * size #n by #n.
   */
  void reinit (const size_type n);

  /**
   * Number of rows of this
   * matrix. For the present
   * matrix, the number of rows and
   * columns are equal, of course.
   */
  size_type m () const;

  /**
   * Number of columns of this
   * matrix. For the present
   * matrix, the number of rows and
   * columns are equal, of course.
   */
  size_type n () const;

  /**
   * Matrix-vector
   * multiplication. For the
   * present case, this of course
   * amounts to simply copying the
   * input vector to the output
   * vector.
   */
  template <class VECTOR1, class VECTOR2>
  void vmult (VECTOR1       &out,
              const VECTOR2 &in) const;

  /**
   * Matrix-vector multiplication
   * with addition to the output
   * vector. For the present case,
   * this of course amounts to
   * simply adding the input
   * vector to the output vector.
   */
  template <class VECTOR1, class VECTOR2>
  void vmult_add (VECTOR1       &out,
                  const VECTOR2 &in) const;

  /**
   * Matrix-vector multiplication
   * with the transpose matrix. For
   * the present case, this of
   * course amounts to simply
   * copying the input vector to
   * the output vector.
   */
  template <class VECTOR1, class VECTOR2>
  void Tvmult (VECTOR1       &out,
               const VECTOR2 &in) const;


  /**
   * Matrix-vector multiplication
   * with the transpose matrix,
   * with addition to the output
   * vector. For the present case,
   * this of course amounts to
   * simply adding the input vector
   * to the output vector.
   */
  template <class VECTOR1, class VECTOR2>
  void Tvmult_add (VECTOR1       &out,
                   const VECTOR2 &in) const;
private:

  /**
   * Number of rows and columns of
   * this matrix.
   */
  size_type size;
};




// ------------------------- inline and template functions -------------
#ifndef DOXYGEN


inline
IdentityMatrix::IdentityMatrix ()
  :
  size (0)
{}



inline
IdentityMatrix::IdentityMatrix (const size_type n)
  :
  size (n)
{}



inline
void
IdentityMatrix::reinit (const size_type n)
{
  size = n;
}



inline
IdentityMatrix::size_type
IdentityMatrix::m () const
{
  return size;
}



inline
IdentityMatrix::size_type
IdentityMatrix::n () const
{
  return size;
}



template <class VECTOR1, class VECTOR2>
inline
void
IdentityMatrix::vmult (VECTOR1       &out,
                       const VECTOR2 &in) const
{
  Assert (out.size() == size, ExcDimensionMismatch (out.size(), size));
  Assert (in.size() == size, ExcDimensionMismatch (in.size(), size));

  out = in;
}



template <class VECTOR1, class VECTOR2>
inline
void
IdentityMatrix::vmult_add (VECTOR1       &out,
                           const VECTOR2 &in) const
{
  Assert (out.size() == size, ExcDimensionMismatch (out.size(), size));
  Assert (in.size() == size, ExcDimensionMismatch (in.size(), size));

  out += in;
}



template <class VECTOR1, class VECTOR2>
inline
void
IdentityMatrix::Tvmult (VECTOR1       &out,
                        const VECTOR2 &in) const
{
  Assert (out.size() == size, ExcDimensionMismatch (out.size(), size));
  Assert (in.size() == size, ExcDimensionMismatch (in.size(), size));

  out = in;
}



template <class VECTOR1, class VECTOR2>
inline
void
IdentityMatrix::Tvmult_add (VECTOR1       &out,
                            const VECTOR2 &in) const
{
  Assert (out.size() == size, ExcDimensionMismatch (out.size(), size));
  Assert (in.size() == size, ExcDimensionMismatch (in.size(), size));

  out += in;
}


#endif

/**@}*/

DEAL_II_NAMESPACE_CLOSE

#endif


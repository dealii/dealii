// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013, 2015 by the deal.II authors
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



/**
 * @defgroup Matrices Matrix classes
 *
 * deal.II comes with a number of different matrix classes, tailored to the
 * various purposes for which matrices are used. For example, there are full
 * matrices, sparse matrices using different storage schemes, matrices
 * composed of individual blocks, and matrices implemented as interfaces to
 * other linear algebra classes. As far as possible, all these implementations
 * share a common interface that contains at least the operations necessary to
 * write iterative linear solvers (see @ref Solvers), but also element-wise
 * access to read from and write to a matrix.
 *
 * This module is split into different parts. @ref Matrix1 "Basic matrices"
 * contains all the matrix classes actually storing entries. @ref Matrix2
 * "Derived matrices", on the other hand, only use basic matrices, but
 * implement certain operations on them. For example, TransposeMatrix provides
 * a matrix-vector multiplication that acts as if the underlying matrix had
 * been transposed, without actually ever storing a transposed matrix.
 *
 * @ref Preconditioners are matrix classes as well, since they perform linear
 * operations on vectors.
 *
 * @ingroup LAC
 */

/**
 * Template for matrix classes.
 *
 * @note This is a description of the expectations on the MATRIX
 * template argument. It is not a class by itself.
 *
 * Depending on where the MATRIX is used, its interface is expected to
 * conform to one of the following groups, which can also be found in
 * the overview.
 *
 * <h3>Solver interface</h3>
 *
 * Functions in this group are the minimal interface to use the MATRIX
 * in a linear Solver. Solvers use a matrix only as a linear operator,
 * that is, they map a vector to another. To this end, we either
 * multiply with the matrix itself or its transpose.
 *
 * The function vmult() is the bare necessity in this group. Some
 * solvers use Tvmult() as well, in which case it needs to be
 * implemented. Some derived matrices like PointerMatrix require its
 * existence, in which case it can be implemented empty with an
 * assertion <code>Assert(false, ExcNotImplemented())</code>.
 *
 * If vmult_add() and Tvmult_add() are missing, PointerMatrixAux can
 * be used to provide the missing functionality without implementing
 * it by hand.
 *
 * @ingroup LAC
 */
template <class VECTOR>
class MATRIX
{
  public:
                                     /**
                                      * @name Solver interface
                                      */
                                     /*@{*/
                                     /**
                                      * The matrix vector product $u = Av$.
                                      */
    void vmult(VECTOR& u, const VECTOR& v) const;
                                     /**
                                      * The matrix vector product $u = A^Tv$.
                                      */
    void Tvmult(VECTOR& u, const VECTOR& v) const;
                                     /**
                                      * The matrix vector product $u += Av$.
                                      */
    void vmult_add(VECTOR& u, const VECTOR& v) const;
                                     /**
                                      * The matrix vector product $u += A^Tv$.
                                      */
    void Tvmult_add(VECTOR& u, const VECTOR& v) const;
                                     /*@}*/
};

/**
 * @defgroup Matrix1 Basic matrices
 *
 * These are the actual matrix classes provided by deal.II. It is possible to
 * store values in them and retrieve them. Furthermore, they provide the full
 * interface required by linear solvers (see @ref Solvers).
 *
 * Among the matrices in this group are full matrices, different sparse
 * matrices, and block matrices. In addition, some of the classes in the
 * interfaces to other linear algebra libraries (for example the
 * PETScWrappers) are matrices.
 *
 * Most of the deal.II sparse matrix classes are separated from their sparsity
 * patterns, to make storing several matrices with the same sparsity pattern
 * more efficient. See @ref Sparsity for more information.
 *
 * @ingroup Matrices
 */


/**
 * @defgroup Matrix2 Derived matrices
 *
 * These matrices are built on top of the basic matrices. They perform special
 * operations using the interface defined in @ref Solvers.
 *
 * @ingroup Matrices
 */

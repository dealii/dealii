// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
 * This documentation topic is split into different parts. @ref Matrix1 "Basic matrices"
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
 * operations using the interface defined by
 * @ref ConceptMatrixType "the MatrixType concept".
 *
 * @ingroup Matrices
 */

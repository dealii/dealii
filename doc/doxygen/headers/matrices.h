// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

/**
 * @defgroup Sparsity Sparsity patterns
 *
 * In deal.II, sparsity patterns are typically separated from the actual
 * sparse matrices (with the exception of the SparseMatrixEZ class and some
 * classes from interfaces to external libraries such as PETSc). The reason is
 * that one often has several matrices that share the same sparsity pattern;
 * examples include the stiffness and mass matrices necessary for time
 * stepping schemes, or the left and right hand side matrix of generalized
 * eigenvalue problems. It would therefore be wasteful if each of them had to
 * store their sparsity pattern separately.
 *
 * Consequently, deal.II has sparsity pattern classes that matrix classes
 * build on. There are two main groups of sparsity pattern classes, as
 * discussed below:
 *
 *
 * <h4>"Static" sparsity patterns</h4>
 *
 * The main sparse matrix class in deal.II, SparseMatrix, only stores a value
 * for each matrix entry, but not where these entries are located. For this,
 * it relies on the information it gets from a sparsity pattern object
 * associated with this matrix. This sparsity pattern object must be of type
 * SparsityPattern.
 *
 * Because matrices are large objects and because it is comparatively
 * expensive to change them, SparsityPattern objects are built in two phases:
 * first, in a "dynamic" phase, one allocates positions where one expects
 * matrices built on it to have nonzero entries; in a second "static" phase,
 * the representation of these nonzero locations is "compressed" into the
 * usual Compressed Row Storage (CSR) format. After this, no new nonzero
 * locations can be added any more. Only after compression can a sparsity
 * pattern be associated to a matrix, since the latter requires the efficient
 * compressed data format of the former. Building a sparsity pattern during
 * the dynamic phase often happens with the DoFTools:make_sparsity_pattern()
 * function. Although this may appear a restriction, it is typically not a
 * significant problem to first build a sparsity pattern and then to write
 * into the matrix only in the previously allocated locations, since in finite
 * element codes it is normally quite clear which elements of a matrix can
 * possibly be nonzero and which are definitely zero.
 *
 * The advantage of this two-phase generation of a sparsity pattern is that
 * when it is actually used with a matrix, a very efficient format is
 * available. In particular, the locations of entries are stored in a linear
 * array that allows for rapid access friendly to modern CPU types with deep
 * hierarchies of caches. Consequently, the static SparsityPattern class is
 * the only one on which deal.II's main SparseMatrix class can work.
 *
 * The main drawback of static sparsity patterns is that their efficient
 * construction requires a reasonably good guess how many entries each of the
 * rows may maximally have. During the actual construction, for example in the
 * DoFTools:make_sparsity_pattern() function, only at most as many entries can
 * be allocated as previously stated. This is a problem because it is often
 * difficult to estimate the maximal number of entries per row. Consequently,
 * a common strategy is to first build and intermediate sparsity pattern that
 * uses a less efficient storage scheme during construction of the sparsity
 * pattern and later copy it directly into the static, compressed form. Most
 * tutorial programs do this, starting at step-2 (see also, for example the
 * step-11, step-18, and step-27 tutorial programs).
 *
 *
 * <h4>"Dynamic" or "compressed" sparsity patterns</h4>
 *
 * As explained above, it is often complicated to obtain good estimates for
 * the maximal number of entries in each row of a sparsity
 * pattern. Consequently, any attempts to allocate a regular SparsityPattern
 * with bad estimates requires huge amounts of memory, almost all of which
 * will not be used and be de-allocated upon compression.
 *
 * To avoid this, deal.II contains a number of "dynamic" or "compressed"
 * sparsity patterns that only allocate as much memory as necessary to hold
 * the currently added entries. While this saves much memory compared to the
 * worst-case behavior mentioned above, it requires the use of less efficient
 * storage schemes for insertion of elements, and the frequent allocation of
 * memory often also takes significant compute time. The tradeoff to avoid
 * excessive memory allocation cannot be avoided, however.
 *
 * The following classes implement this "dynamic" memory scheme in deal.II:
 * - CompressedSparsityPattern
 * - CompressedSimpleSparsityPattern
 * - CompressedSetSparsityPattern
 *
 * These classes have different performance characteristics and memory
 * requirements. Which one is the "best" changes from case to case because
 * it is dependent on the number of dofs, the number of couplings per dof,
 * the strategy used to insert and the amount of memory available.
 *
 * CompressedSparsityPattern and CompressedSimpleSparsityPattern are very
 * similar, where CompressedSimpleSparsityPattern trades some memory (requires
 * up to twice the memory in the worst case) for additional speed which is
 * noticeable in cases with many nonzero entries. CompressedSetSparsityPattern
 * on the other hand uses a lot more memory but may perform better in cases with
 * many nonzero entries per row. See for example the step-27
 * and step-22 tutorial programs.
 *
 * As a guideline you should start using CompressedSparsityPattern and try the
 * other variants if you run into problems. Switching between them should be as
 * simple as changing the class name because all classes have the same interface
 * for adding entries.
 * In either case, these classes are typically used in the following way
 * (replace the class CompressedSparsityPattern with a different one from above):
 * @verbatim
 * CompressedSparsityPattern compressed_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  compressed_pattern);
 * constraints.condense (compressed_pattern);
 *
 * SparsityPattern final_sparsity_pattern;
 * final_sparsity_pattern.copy_from (compressed_pattern);
 * @endverbatim
 *
 * The intermediate, compressed sparsity pattern is directly copied into the
 * "compressed" form of the final static pattern.
 *
 * <h4>Dynamic block sparsity patterns</h4>
 *
 * The following classes implement an array of dynamic sparsity patterns:
 * - BlockCompressedSparsityPattern
 * - BlockCompressedSetSparsityPattern
 * - BlockCompressedSimpleSparsityPattern
 *
 * These classes inherit the same tradeoffs regarding their efficiency from
 * their non-block classes (see above). See their documentation and
 * step-22 for more information.
 *
 * @ingroup Matrices
 */

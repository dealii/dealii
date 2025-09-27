// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2015 by the deal.II authors
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
 * @defgroup Sparsity Sparsity patterns
 *
 * Almost all finite element formulations lead to matrices that are
 * "sparse", i.e., for which the number of nonzero elements per row is
 * (i) relatively small compared to the overall size of the matrix,
 * and (ii) bounded by a fixed number that does not grow if the mesh
 * is refined. For such cases, it is more efficient to not store
 * <i>all</i> elements of the matrix, but only those that are actually
 * (or may be) nonzero. This requires storing, for each row, the
 * column indices of the nonzero entries (we call this the "sparsity
 * pattern") as well as the actual values of these nonzero
 * entries. (In practice, it sometimes happens that some of the
 * nonzero values are, in fact, zero. Sparsity patterns and sparse
 * matrices only intend to provision space for entries that <i>may</i>
 * be nonzero, and do so at a time when we don't know yet what values
 * these entries will ultimately have; they may have a zero value if a
 * coefficient or cell happens to have particular values.)
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
 * usual Compressed Sparse Row (CSR) format. After this, no new nonzero
 * locations may be added. Only after compression can a sparsity pattern be
 * associated to a matrix, since the latter requires the efficient compressed
 * data format of the former. Building a sparsity pattern during the dynamic
 * phase often happens with the DoFTools::make_sparsity_pattern()
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
 * DoFTools::make_sparsity_pattern() function, only at most as many entries can
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
 * To avoid this, deal.II contains a "dynamic" or "compressed" sparsity
 * pattern called DynamicSparsityPattern that only allocates as much memory as
 * necessary to hold the currently added entries. While this saves much memory
 * compared to the worst-case behavior mentioned above, it requires the use of
 * less efficient storage schemes for insertion of elements, and the frequent
 * allocation of memory often also takes significant compute time. The
 * tradeoff to avoid excessive memory allocation cannot be avoided, however.
 *
 * The class is typically used in the following way
 * @verbatim
 * DynamicSparsityPattern dsp (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  dsp);
 * constraints.condense (dsp);
 *
 * SparsityPattern final_sparsity_pattern;
 * final_sparsity_pattern.copy_from (dsp);
 * @endverbatim
 *
 * The intermediate, compressed sparsity pattern is directly copied into the
 * "compressed" form of the final static pattern.
 *
 * <h4>Dynamic block sparsity patterns</h4>
 *
 * The class BlockDynamicSparsityPattern implements an array of dynamic
 * sparsity patterns for constructing block matrices. See the documentation and
 * step-22 for more information.
 *
 * @ingroup LAC
 */

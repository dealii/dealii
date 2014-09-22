// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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
 * @page GlobalDoFIndex When to use types::global_dof_index instead of unsigned int
 *
 * deal.II can be configured to use 64-bit indices for degrees of freedom,
 * rather than the usual unsigned integers that default to 32-bit on most
 * current systems. This is necessary since we want to be able to solve
 * problems with more than 4 billion unknowns (the limit of what can be
 * represented with 32-bit unsigned integers). At the same time, we do not
 * want to indiscriminately replace all integers in deal.II with 64-bit
 * versions, since this would increase memory use in many places where we
 * represent quantities that will most definitely not be larger than 4 billion.
 *
 * The data type we define for these indices to keep the bulk
 * of the code base free of <code>\#ifdef</code>s is types::global_dof_index.
 * If deal.II is configured as normal, this type is <code>unsigned int</code>,
 * but can be switched to <code>unsigned long long int</code> if the right
 * flag is provided (see the ReadMe file). This page is intended to clarify
 * when types::global_dof_index must be used and when one can use a regular
 * unsigned integer:
 *
 * <dl>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexBlockIndices
 * <b>BlockIndices</b></dt>
 * <dd>
 * The number of blocks is an unsigned int because the number is expected to
 * be low, i.e less than four billions. However, the size of the block is a
 * types::global_dof_index because each block can be arbitrary large.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexCell <b>Cell</b></dt>
 * <dd>
 * The ID of cell is not unique: Cells with different levels of refinement
 * and/or on different processors can have the same ID. Thus, all the data
 * associated to cells can be unsigned int because on a single processor,
 * one one mesh level, there will definitely not be more than 4 billion
 * cells.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexDoFHandler
 * <b>DoFHandler</b></dt>
 * <dd>
 * The ID of each degree of freedom is unique in a parallel computation.
 * Therefore, degrees of freedom
 * are types::global_dof_index.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexFullMatrix
 * <b>FullMatrix</b></dt>
 * <dd>
 * The numbers of row and column are types::global_dof_index even if it is not
 * expected that someone will create a FullMatrix with so many entries.
 * However, some functions of ConstraintMatrix are templated on the matrix
 * type and thus, the
 * size of a FullMatrix has to be of the same type than the size of
 * SparseMatrix.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexSparseMatrix
 * <b>SparseMatrix</b></dt>
 * <dd>
 * The size of SparseMatrix can be arbitrary large and it is conceivable that
 * with sufficient memory on a single node, one may generate a matrix with
 * more than 4 billion rows or columns. Therefore, types::global_dof_index is
 * used. However, even for a large complex problem we can solve now, it is not
 * reasonable to expect the number of non-zero entries in a sparse matrix to
 * go over four billion. Thus, we still use unsigned int for, e.g.,
 * SparsityPattern::row_lengths and similar functions.
 * </dd>
 *
 * </dl>
 */

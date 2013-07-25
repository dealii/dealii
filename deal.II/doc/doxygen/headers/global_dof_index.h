// ---------------------------------------------------------------------
// $Id$
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
 * When the 64-bits version of deal.II is used, it becomes necessary to
 * declare as types::global_dof_index (unsigned long long int) instead of
 * unsigned int. Here, we want to clarify when types::global_dof_index must be
 * used.
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
 * The ID of cell is not unique. Cells with different levels of refinement
 * and/or on different processors can have the same ID. Thus, all the data
 * associated to cells can be unsigned int.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexDoFHandler
 * <b>DoFHandler</b></dt>
 * <dd>
 * The ID of each degree of freedom is unique. Therefore, degrees of freedom
 * are types::global_dof_index.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexFullMatrix
 * <b>FullMatrix</b></dt>
 * <dd>
 * The numbers of row and column are types::global_dof_index even if it is not
 * expected that someone will create a FullMatrix with so many entries.
 * However, ConstraintMatrix is templated on the matrix type and thus, the
 * size of a FullMatrix has to be of the same type than the size of
 * SparseMatrix.
 * </dd>
 *
 * <dt class="glossary">@anchor GlobalDoFIndexSparseMatrix
 * <b>SparseMatrix</b></dt>
 * <dd>
 * The size of SparseMatrix can be arbitrary large therefore,
 * types::global_do_index is used. However, even for a large complex problem we
 * can solve now, there is no reason for the number of non-zero entries in a 
 * sparse matrix to go over four billions. Thus, we still use unsigned int
 * for, e.g., row_lengths in the object.
 * </dd>
 *
 * </dl>
 */

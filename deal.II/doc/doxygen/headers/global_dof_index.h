//-------------------------------------------------------------------------
//    $Id: coding_conventions.h 29312 2013-04-17 13:00:07Z turcksin $
//    Version: $Name$
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @page GlobalDoFIndex When to use types::global_dof_index instead of
 * unsigned int
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
 * The size SparseMatrix can be arbitrary large therefore,
 * types::global_do_index is used.
 */

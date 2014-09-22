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
 * @defgroup mg Multilevel support
 *
 * Classes that have to do with multigrid algorithms.
 *
 * The main class with implementation of the multigrid scheme is
 * Multigrid with its function Multigrid::cycle(). It uses the
 * following abstract classes in order to perform the multigrid cycle:
 *
 * <ol>
 * <li> MGMatrixBase contains the level matrices with a fairly general
 * implementation in MGMatrix
 * <li> MGCoarseGridBase is the solver on the coarsest level.
 * <li> MGSmootherBase performs smoothing on each level.
 * <li> MGTransferBase organizes the transfer between levels.
 * </ol>
 *
 * Additionally, there is a class PreconditionMG, which is a wrapper
 * around Multigrid with the standard interface of deal.II @ref
 * Preconditioners. PreconditionMG also uses the classes inheriting
 * from MGTransferBase, for instance MGTransferPrebuilt, where it uses
 * MGTransferPrebuilt::copy_to_mg() and
 * MGTransferPrebuilt::copy_from_mg_add(), which transfer between the
 * global vector and the level vectors.
 *
 * Finally, we have several auxiliary classes, namely MGLevelObject,
 * which stores an object on each level*
 * 
 * See the step-16 and step-39 example programs on how to use this
 * functionality.

 * <h3>Multigrid and hanging nodes</h3>
 *
 * Using multigrid methods on adaptively refined meshes involves
 * more infrastructure than with regular refinement. First, in order
 * to keep the complexity optimal, we need to decide how to do the
 * smoothing on each level. And to this end, we have to define what a
 * level is in the sense of multilevel decomposition.
 *
 * First, we define that a level in the multigrid sense is constituted
 * by all cells of a certain level in the mesh hierarchy. Thus,
 * smoothing on a certain level is restricted to the subdomain which
 * consists of cells of this level or finer. This is usually referred
 * to as local smoothing. The advantage of this definition is, that
 * level matrices for the multigrid scheme can be assembled easily by
 * traversing to all cells of a certain level, and that these level
 * matrices do not contain hanging nodes.
 *
 * The disadvantage of this decomposition is, that we need additional
 * matrices to handle the issues that arise at refinement
 * edges. Furthermore, the treatment is different, depending on
 * whether the method is continuous (thus having degrees of freedom on
 * the refinement edge) or discontinuous (employs flux matrices at the
 * refinement edge). While these matrices are small, we have to
 * assemble them and notify the multigrid method of them.
 */

/**
 * This namespace contains the reimplementation of multilevel support
 * after we know what is needed in the context of local refinement and
 * block systems.
 *
 * @ingroup mg
 */
namespace mg
{
}

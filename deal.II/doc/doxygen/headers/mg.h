//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

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
 * which stores an object on each level
 * 
 * See the step-16 example program on how to use this
 * functionality.
 */


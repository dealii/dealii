//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------
/**
 * @defgroup Solvers Linear solver classes
 *
 * This module groups iterative and direct solvers, eigenvalue solvers, and
 * some control classes. All iterative solvers inherit from the class template
 * Solver, which provides some basic maintenance methods.
 *
 * The number of iteration steps of iterative solvers is controlled by
 * objects of class SolverControl or its derived class
 * ReductionControl.
 *
 * All solvers receive the matrix and vector classes as template
 * arguments. Therefore, any objects defining the interface described in the
 * documentation of Solver are admissible. These requirements are listed in
 * the documentation page of the Solver class.
 *
 * If detected during configuration (see the ReadMe file), some sparse direct
 * solvers are also supported.
 *
 * @ingroup LAC
 */

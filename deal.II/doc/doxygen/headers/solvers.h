//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006, 2007 by the deal.II authors
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
 * some control classes. All these classes operate on objects of the
 * @ref Matrices "matrix" and @ref Vectors "vector classes" defined in deal.II.
 *
 * In order to work properly, solvers that take matrix and vector classes as
 * template arguments require that these classes satisfy a certain minimal
 * interface that can be used from inside the solver. For iterative solvers,
 * this interface is defined in the Solver class. In addition, solvers are
 * controlled using objects of classes that are derived from the SolverControl
 * class (for example its derived class ReductionControl), in order to
 * determine the maximal number of iterations or a desired tolerance.
 *
 * If detected during configuration (see the ReadMe file), some sparse direct
 * solvers are also supported.
 *
 * @ingroup LAC
 */

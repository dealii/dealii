//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup Solvers Linear Solver classes
 *
 * Here we collect iterative solvers and some control classes. All
 * solvers inherit from the class template Solver, which provides some
 * basic maintenance methods.
 *
 * The number of iteration steps of all solvers is controlled by
 * objects of class SolverControl or its derived class
 * ReductionControl.
 */

/**
 * @defgroup Matrices Matrix classes
 * @{
 *
 * All matrices in this library have a common minimal interface,
 * defined through MATRIX. This interface consists of functions for
 * multiplication with appropriate vectors.
 *
 */

/**
 * @defgroup Matrix1 Basic matrices
 *
 * These are the actual matrix classes provided by deal.II. They all
 * provide some interface for allocating memory and entering values.
*/


/**
 * @defgroup Matrix2 Derived matrices
 */

/**
 * @defgroup Preconditioners Preconditioners
 * @}
 */

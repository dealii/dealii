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
 * @defgroup Matrices Matrix classes
 *
 * All matrices in this library have a common minimal interface, defined
 * through MATRIX (see Solver documentation). This interface consists of
 * functions for multiplication with appropriate vectors.
 *
 * We split this module into several parts. Basic matrices are all the matrix
 * classes actually storing their entries. Derived matrices use basic
 * matrices, but change the meaning of the matrix-vector multiplication.
 *
 * Preconditioners are matrix classes as well, since they perform linear
 * operations on vectors.
 *
 * @author Guido Kanschat, 2003
 *
 * @ingroup LAC
 */

/**
 * @defgroup Matrix1 Basic matrices
 *
 * These are the actual matrix classes provided by deal.II. It is possible to
 * store values in them and retrieve them. Furthermore, they provide the full
 * interface required by linear solvers (see Solver).
 *
 * @author Guido Kanschat, 2003
 *
 * @ingroup Matrices
 */


/**
 * @defgroup Matrix2 Derived matrices
 *
 * These matrices are built on top of the basic matrices. They perform special
 * operations using the interface defined in Solver.
 *
 * @author Guido Kanschat, 2003
 *
 * @ingroup Matrices
 */

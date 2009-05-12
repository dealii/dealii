//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------


/**
 * @defgroup Preconditioners Preconditioners
 *
 * Preconditioners are used to accelerate the iterative solution of linear
 * systems. Typical preconditioners are Jacobi, Gauss-Seidel, or SSOR, but the
 * library also supports more complex ones such as Vanka or incomplete LU
 * decompositions (ILU). In addition, sparse direct solvers can be used as
 * preconditioners when available.
 * 
 * In order to be used by deal.II solvers, preconditioners must
 * conform to the standard matrix interface. Solvers use the function
 * <tt>vmult()</tt> of the preconditioner. Some solvers may also use
 * <tt>Tvmult()</tt>.
 *
 * @ingroup LAC
 * @ingroup Matrices
 */

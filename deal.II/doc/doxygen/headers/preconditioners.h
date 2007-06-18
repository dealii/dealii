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
 * @defgroup Preconditioners Preconditioners
 *
 * Preconditioners are used to accelerate the iterative solution of linear
 * systems. Typical preconditioners are Jacobi, Gauss-Seidel, or SSOR, but the
 * library also supports more complex ones such as Vanka or incomplete LU
 * decompositions (ILU). In addition, sparse direct solvers can be used as
 * preconditioners when available.
 * 
 * In principle, and in the mathematical literature, preconditioners are
 * treated as matrices in the sense that one can do a matrix-vector
 * multiplication with them. On the other hand, one doesn't usually have an
 * element-by-element representation of these matrices, only their action on a
 * vector. The preconditioner classes therefore often have a <tt>vmult()</tt>
 * function that symbolizes the ability to perform matrix-vector
 * multiplications, just like the real matrix classes.
 *
 * @ingroup LAC
 * @ingroup Matrices
 */

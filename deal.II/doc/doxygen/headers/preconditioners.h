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
 * When talking of preconditioners, we usually expect them to be used
 * in Krylov-space methods. Nevertheless, the concept becomes clearer
 * in the standard linear defect correction
 * @f[
 *  x^{k+1} = x^k - P^{-1} \bigl(A x^k - b\bigr),
 * @f]
 * where <i>P<sup>-1</sup></i> is the preconditioner. Thus,
 * preconditioning amounts to applying a linear operator to the
 * residual. For this reason, the interface of preconditioners equals
 * the one for matrices.
 *
 * <h3>The interface</h3>
 *
 * In this section, we discuss the interface preconditioners usually
 * have to provide to work inside the deal.II library.
 *
 * <h4>Initialization</h4>
 *
 * In order to be able to be stored in containers, all preconditioners
 * have a constructor with no arguments. Since this will typically
 * produce a useless object, all preconditioners have a function
 * @code
 *   void initialize (...)
 * @endcode
 *
 * This function receives the matrix to be preconditioned as well as
 * additional required parameters and sets up the internal structures
 * of the preconditioner.
 *
 * <h4>Application of the preconditioner</h4>
 *
 * Preconditioners in deal.II are just considered linear operators.
 * Therefore, in order to be used by deal.II solvers, preconditioners must
 * conform to the standard matrix interface, namely the functions
 * @code
 * void  vmult (VECTOR& dst, const VECTOR& src) const;
 * void  Tvmult (VECTOR& dst, const VECTOR& src) const;
 * @endcode
 * Solvers use the function
 * <tt>vmult()</tt> of the preconditioner. Some solvers may also use
 * <tt>Tvmult()</tt>.
 *
 * <h4>Relaxation methods</h4>
 *
 * Additional to the interface described below, some preconditioners
 * like SOR and Jacobi have benn known as iterative methods
 * themselves. For them, an additional interface exists, consisting of
 * the functions
 * @code
 * void  step (VECTOR& dst, const VECTOR& src) const;
 * void  Tstep (VECTOR& dst, const VECTOR& src) const;
 * @endcode
 *
 * @ingroup LAC
 * @ingroup Matrices
 */

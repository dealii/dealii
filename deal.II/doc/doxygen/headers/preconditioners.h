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
 * in Krylov-space methods. In that case, the act as operators: given
 * a vector $x$, produce the result $y=P^{-1}x$ of the multiplication
 * with the preconditioning operator $P^{-1}$.
 *
 * However, some preconditioners can also be used
 * in the standard linear defect correction iteration,
 * @f[
 *  x^{k+1} = x^k - P^{-1} \bigl(A x^k - b\bigr),
 * @f]
 * where <i>P<sup>-1</sup></i> is again the preconditioner. Thus,
 * preconditioning amounts to applying a linear operator to the
 * residual.
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
 * These functions apply the preconditioning operator to the source
 * vector $src$ and return the result in $dst$ as $dst=P^{-1}src$ or
 * $dst=P^{-T}src$. Preconditioned iterative
 * dolvers use these <tt>vmult()</tt> functions of the preconditioner.
 * Some solvers may also use <tt>Tvmult()</tt>.
 *
 * <h4>Relaxation methods</h4>
 *
 * Additional to the interface described above, some preconditioners
 * like SOR and Jacobi have been known as iterative methods
 * themselves. For them, an additional interface exists, consisting of
 * the functions
 * @code
 * void  step (VECTOR& dst, const VECTOR& rhs) const;
 * void  Tstep (VECTOR& dst, const VECTOR& rhs) const;
 * @endcode
 *
 * Here, $src$ is a residual vector and $dst$ is the iterate that is
 * supposed to be updated. In other words, the operation performed by
 * these functions is
 * $dst = dst - P^{-1} (A dst - rhs)$ and $dst = dst - P^{-T} (A dst - rhs)$. The
 * functions are called this way because they perform <i>one step</i>
 * of a fixed point (Richardson) iteration. Note that preconditioners
 * store a reference to the original matrix $A$ during initialization.
 *
 * @ingroup LAC
 * @ingroup Matrices
 */

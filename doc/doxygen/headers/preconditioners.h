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
 * @defgroup Preconditioners Preconditioners and Relaxation Operators
 *
 * <h3>Preconditioners</h3>
 *
 * Preconditioners are used to accelerate the iterative solution of linear
 * systems. Typical preconditioners are Jacobi, Gauss-Seidel, or SSOR, but the
 * library also supports more complex ones such as Vanka or incomplete LU
 * decompositions (ILU). In addition, sparse direct solvers can be used as
 * preconditioners when available.
 *
 * Broadly speaking, preconditioners are operators, which are
 * multiplied with a matrix to improve conditioning. The idea is, that
 * the preconditioned system <i>P<sup>-1</sup>Ax = P<sup>-1</sup>b</i>
 * is much easier to solve than the original system <i>Ax = b</i>.What
 * this means exactly depends on the structure of  the matrix and
 * cannot be discussed here in generality. For symmetric, positive
 * definite matrices <i>A</i> and <i>P</i>, it means that the spectral
 * condition number (the quotient of greatest and smallest eigenvalue)
 * of <i>P<sup>-1</sup>A</i> is much smaller than the one of <i>A</i>.
 *
 * At hand of the simplest example, Richardson iteration, implemented
 * in SolverRichardson, the preconditioned iteration looks like
 * @f[
 *  x^{k+1} = x^k - P^{-1} \bigl(A x^k - b\bigr).
 * @f]
 * Accordingly, preconditioning amounts to applying a linear operator to the
 * residual, and consequently, the action of the preconditioner
 * <i>P<sup>-1</sup></i> is implemented as <tt>vmult()</tt>. The
 * generic interface is like for matrices
 * @code
 * class PRECONDITIONER
 * {
 *   template <class VECTOR>
 *   void vmult(VECTOR& dst, const VECTOR& src) const;
 *
 *   template <class VECTOR>
 *   void Tvmult(VECTOR& dst, const VECTOR& src) const;
 * }
 * @endcode
 * It is implemented in all the preconditioner classes in this module.
 * 
 * When used
 * in Krylov space methods, it is up to the method, whether it simply
 * replaces multiplications with <i>A</i> by those with
 * <i>P<sup>-1</sup>A</i> (for instance SolverBicgstab), or does more
 * sophisticated things. SolverCG for instance uses
 * <i>P<sup>-1</sup></i> to define an inner product, which is the
 * reason why it requires a symmetric, positive definite operator <i>P</i>.
 *
 * <h3>Relaxation methods</h3>
 *
 * Many preconditioners rely on an additive splitting <i>A = P - N</i>
 * into two matrices. In this case, the iteration step of the
 * Richardson method above can be simplified to
 * @f[
 *  x^{k+1} = P^{-1} \bigl(N x^k + b\bigr),
 * @f]
 * thus avoiding multiplication with <i>A</i> completely. We call
 * operators mapping the previous iterate <i>x<sup>k</sup></i> to the
 * next iterate in this way relaxation operators. Their generic
 * interface is
 * @code
 * class RELAXATION
 * {
 *   template <class VECTOR>
 *   void step(VECTOR& newstep, const VECTOR& oldstep, const VECTOR& rhs) const;
 *
 *   template <class VECTOR>
 *   void Tstep(VECTOR& newstep, const VECTOR& oldstep, const VECTOR& rhs) const;
 * }
 * @endcode
 * The classes with names starting with <tt>Relaxation</tt> in this
 *   module implement this interface, as well as the preconditioners
 *   PreconditionJacobi, PreconditionSOR, PreconditionSSORP,
 *   reconditionBlockJacobi, PreconditionBlockSOR, and
 *   PreconditionBlockSSOR.
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

//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup Vectors Vector classes
 *
 * Here, we list all the classes that can be used as VECTOR in linear solvers
 * (see Solver) and for matrix operations.
 */

/**
 * @defgroup Matrices Matrix classes
 * @{
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
 */

/**
 * @defgroup Matrix1 Basic matrices
 *
 * These are the actual matrix classes provided by deal.II. It is possible to
 * store values in them and retrieve them. Furthermore, they provide the full
 * interface required by linear solvers (see Solver).
 *
 * @author Guido Kanschat, 2003
 */


/**
 * @defgroup Matrix2 Derived matrices
 *
 * These matrices are built on top of the basic matrices. They perform special
 * operations using the interface defined in Solver.
 *
 * @author Guido Kanschat, 2003
 */

/**
 * @defgroup Preconditioners Preconditioners
 * In principle, preconditioners are matrices themselves. But, since they are
 * only used in a special context, we put them in a separate group.
 *
 * @author Guido Kanschat, 2003
 */

 /*@}*/

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
 *
 * All solvers receive the matrix and vector classes as template
 * arguments. Therefore, any objects defining the interface described in the
 * documentation of Solver are admissible.
 *
 * @author Guido Kanschat, 2003
 */

/**
 * @defgroup VMemory Vector memory management
 *
 * A few classes that are used to avoid allocating and deallocating vectors in
 * iterative procedures. These methods all use an object of the base class
 * VectorMemory to get their auxiliary vectors. Unfortunately, the
 * intelligence put into both implementations provided right now is not
 * overwhelming.
 *
 * @author Guido Kanschat, 2003
 */

/**
 * @defgroup PETSc PETSc
 * The classes in this module are wrappers around functionality provided by
 * the PETSc library. They provide a modern object-oriented interface that is
 * compatible with the interfaces of the other linear algebra classes in
 * deal.II.
 *
 * These classes are only available if a PETSc installation was detected
 * during configuration of deal.II. Refer to the README file for more details
 * about this.
 *
 * @author Wolfgang Bangerth, 2004
 */

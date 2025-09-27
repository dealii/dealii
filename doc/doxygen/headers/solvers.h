// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
 * @defgroup Solvers Linear solver classes
 *
 * This page groups iterative and direct solvers, eigenvalue solvers, and
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

// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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
 * @defgroup GenericLinearAlgebra Generic Linear Algebra
 *
 * This module encapsulates several namespaces that aggregate the various
 * vectors, sparse matrices, solvers and preconditioners for sparse linear
 * systems. So far it includes a subset of the native deal.II data types
 * and, if deal.II is built with them enabled, PETSc and Trilinos data for
 * both serial and MPI types as well.
 *
 * The deal.II type definitions can be found in the LinearAlgebraDealII
 * namespace (supporting real-valued problems) and LinearAlgebraDealII::Complex
 * (for complex-valued problems). Similarly, those for Trilinos are in
 * LinearAlgebraTrilinos (serial) and LinearAlgebraTrilinos::MPI (parallel),
 * and those for PETSc are in LinearAlgebraPETSc (serial) and
 * LinearAlgebraPETSc::MPI (parallel). Note that the underlying number type
 * for Trilinos and PETSc linear algebra is determined at the time of
 * compilation of the deal.II library.
 *
 * The main idea behind offering such a collection of namespaces is the
 * facilitation of easy switching between the variety options with
 * which one can solve a linear system. An example that takes advantage
 * of this is shown below (this uses Trilinos and PETSc types only):
 *
 * @code
 * #include <deal.II/lac/generic_linear_algebra.h>
 *
 * #define USE_TRILINOS_LA
 * #define USE_MPI_LA
 * namespace LA
 * {
 * #ifdef USE_TRILINOS_LA
 *
 * #if defined(DEAL_II_WITH_MPI) && defined(USE_MPI_LA)
 *   using namespace dealii::LinearAlgebraTrilinos::MPI;
 * #else
 *   using namespace dealii::LinearAlgebraTrilinos;
 * #endif
 *
 * #else
 *
 * #if defined(DEAL_II_WITH_MPI) && defined(USE_MPI_LA)
 *   using namespace dealii::LinearAlgebraPETSc::MPI;
 * #else
 *   using namespace dealii::LinearAlgebraPETSc;
 * #endif
 *
 * #endif
 * }
 * @endcode
 *
 * This set of preprocessor directives select the appropriate
 * namespace based on the existence some predefined identifier,
 * and import them into a local namespace @p LA .
 * Now one can create and use, for example, a generic
 * vector of switchable type by
 *
 * @code
 * const LA::BlockVector solution = ...
 * @endcode
 *
 * @note For this to be used most effectively, the common
 * elements of the vector/sparse matrix/solver/preconditioner interface
 * must be called during initialisation and operation.
 * For examples of this, see step-40, step-50 and step-55.
 *
 * @author Timo Heister 2008, Jean-Paul Pelteret, 2016
 *
 * @ingroup LAC
 */

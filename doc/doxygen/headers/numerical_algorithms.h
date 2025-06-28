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
 * @defgroup numerics Numerical algorithms
 *
 * This topic groups a diverse set of classes that generally implement some
 * sort of numerical algorithm on top all the basic triangulation, DoFHandler,
 * and finite element classes in the library. They are generally unconnected
 * to each other.
 *
 * Some of the classes, like DerivativeApproximation, KellyErrorEstimator and
 * SolutionTransfer, act on solutions already obtained, and compute derived
 * quantities in the first two cases, or help transferring a set of vectors
 * from one mesh to another.
 *
 * The namespaces MatrixCreator, MatrixTools, and VectorTools provide an
 * assortment of services, such as creating a Laplace matrix, projecting or
 * interpolating a function onto the present finite element space, etc.  The
 * difference to the functions in the DoFTools and FETools functions is that
 * they work on vectors (i.e. members of a finite element function space on a
 * given triangulation) or help in the creation of it. On the other hand, the
 * DoFTools functions only act on a given DoFHandler object without reference
 * to a data vector, and the FETools objects generally work with finite
 * element classes but again without any associated data vectors.
 */

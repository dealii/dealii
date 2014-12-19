// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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
 * @defgroup numerics Numerical algorithms
 *
 * This module groups a diverse set of classes that generally implement some
 * sort of numerical algorithm on top all the basic triangulation, DoFHandler,
 * and finite element classes in the library. They are generally unconnected
 * to each other.
 *
 * Some of the classes, like DerivativeApproximation, KellyErrorEstimator and
 * SolutionTransfer, act on solutions already obtained, and compute derived
 * quantities in the first two cases, or help transferring a set of vectors
 * from one mesh to another.
 *
 * The remaining classes MatrixCreator, MatrixTools, and VectorTools provide
 * an assortment of services, such as creating a Laplac matrix, projecting or
 * interpolating a function onto the present finite element space, etc.  The
 * difference to the functions in the DoFTools and FETools functions is that
 * they work on vectors (i.e. members of a finite element function space on a
 * given triangulation) or help in the creation of it. On the other hand, the
 * DoFTools functions only act on a given DoFHandler object without reference
 * to a data vector, and the FETools objects generally work with finite
 * element classes but again without any associated data vectors.
 */

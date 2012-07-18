//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

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

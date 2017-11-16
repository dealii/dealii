// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
 * @defgroup auto_symb_diff Automatic and symbolic differentiation
 *
 * @brief A module dedicated to the implementation of functions and classes that relate
 * to automatic and symbolic differentiation.
 *
 * @todo Hyper-summarize the following list of topics:
 * - Automatic differentiation
 * - Symbolic differentiation
 *
 * @section auto_diff_1 Automatic differentiation
 *
 * @todo Write a short introduction into AD. As a temporary entry, the following links
 * may be enlightening:
 * - <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">Wikipedia article</a>
 * - <a href="https://projects.coin-or.org/ADOL-C/browser/stable/2.6/ADOL-C/doc/adolc-manual.pdf?format=raw#page=1>Adol-C manual</a>
 *
 * @todo Hyper-summarize the following list of topics:
 * - Forward and reverse mode AD
 * - Taped and tapeless AD; expression templates
 *
 * @subsection auto_diff_1_1 Supported automatic differentiation libraries
 *
 * We currently have validated implementations for the following number types
 * and combinations:
 *
 *  - Taped Adol-C (n-differentiable, in theory, but internal drivers for up to second-order
 *    derivatives have been implemented)
 *  - Tapeless Adol-C (once differentiable)
 *  - Tapeless forward-mode Sacado with dynamic memory allocation (once differentiable)
 *  - Tapeless nested forward-mode Sacado (twice differentiable)
 *  - Tapeless reverse-mode Sacado (once differentiable)
 *  - Tapeless nested reverse and forward-mode Sacado (twice differentiable)
 *
 * @subsection auto_diff_1_2 How automatic differentiation is integrated into deal.II
 *
 * Since the interface to each automatic differentiation library is so vastly different,
 * a uniform internal interface to each number has been established. This allows the
 * driver classes (that provide the core functionality, and are introduced in the next
 * section) a consistent mechanism to interact with different auto-differentiation
 * libraries. Specifically, they need to be able to correctly initialize and finalize data
 * that is to be interpreted as the dependent and independent variables of a formula.
 *
 * A summary of the files that implement the interface to the supported auto-differentiable
 * numbers is as follows:
 *
 * - ad_helpers.h: Provides a set of classes to help perform automatic differentiation in a
 *   number of different contexts. These are detailed in \ref auto_diff_1_3.
 * - ad_number_types.h: Introduces an enumeration (called a type code) for the
 *   auto-differentiable number combinations that will be supported by the driver classes.
 *   The rationale behind the use of this somewhat restrictive mechanism is discussed below.
 * - ad_number_traits.h: Declare some internal classes that are to be specialized for
 *   each auto-differentiation library and/or number type. These are subsequently used to
 *   provide a uniform interface to the classes through the NumberTraits and ADNumberTraits
 *   classes which are extensively used throughout of drivers. We also provide some mechanisms
 *   to easily query select properties of these numbers, i.e. some type traits.
 * - adolc_math.h: Extension of the Adol-C math operations that allow these numbers to be used
 *   consistently throughout the library.
 * - adolc_number_types.h: Implementation of the internal classes that define how we
 *   use Adol-C numbers.
 * - adolc_product_types.h: Defines some product and scalar types that allow the use of
 *   Adol-C numbers in conjunction with the Tensor and SymmetricTensor classes.
 * - sacado_math.h: Extension of the sacado math operations that allow these numbers to be used
 *   consistently throughout the library.
 * - sacado_number_types.h: Implementation of the internal classes that define how we
 *   use the supported Sacado numbers.
 * - sacado_product_types.h: Defines some product and scalar types that allow the use of
 *   the supported Sacado numbers in conjunction with the Tensor and SymmetricTensor
 *   classes.
 *
 * By using type codes for each supported number type, we artificially limit the type
 * of auto-differentiable numbers that can be used within the library. This design choice
 * is due to the fact that its not trivial to ensure that each number type is correctly
 * initialized and that all combinations of nested (templated) types remain valid for all
 * operations performed by the library.
 * Furthermore, there are some lengthy functions within the library that are instantiated
 * for the supported number types and have internal checks that are only satisfied when a
 * auto-differentiable number, of which the library has knowledge, is used. This again
 * ensures that the integrity of all computations is maintained.
 * Finally, using a simple enumeration as a class template parameter ultimately makes it
 * really easy to switch between the type used in production code with little to no further
 * amendments required to user code.
 *
 * @subsubsection auto_diff_1_3 User interface to the automatic differentiation libraries
 *
 * @todo Summarize driver classes 
 * - %Quadrature point level
 *   - Scalar mode
 *   - %Vector mode
 * - Cell level functions
 *   - Variational formulations
 *   - Residual linearisation
 *
 */

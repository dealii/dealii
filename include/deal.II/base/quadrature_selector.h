// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_quadrature_selector_h
#define dealii_quadrature_selector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the quadrature rule passed to its constructor as a
 * string. Supported quadratures are QGauss (of all orders), QMidpoint,
 * QMilne, QSimpson, QTrapezoid and QWeddle.
 *
 * This class is useful if you want to use flexible quadrature rules, that are
 * read from a parameter file (see ParameterHandler for this).
 *
 * @ingroup Quadrature
 */
template <int dim>
class QuadratureSelector : public Quadrature<dim>
{
public:
  /**
   * Constructor. Takes the name of the quadrature rule (one of "gauss",
   * "milne", "weddle", etc) and, if it is "gauss", the number of quadrature
   * points in each coordinate direction.
   */
  QuadratureSelector(const std::string &s, const unsigned int order = 0);

  /**
   * This function returns all possible names for quadratures as a list
   * separated by <tt>|</tt>, so that you can use it for the definition of
   * parameter files (see ParameterHandler for details).
   */
  static std::string
  get_quadrature_names();

  /**
   * @addtogroup Exceptions
   * @{
   */


  /**
   * Exception
   */
  DeclException1(ExcInvalidQGaussOrder,
                 int,
                 << "You tried to generate a QGauss object with an invalid "
                 << "number " << arg1
                 << " of quadrature points in each coordinate "
                 << "direction. This number must be greater than or equal "
                 << "to 1.");
  /**
   * Exception
   */
  DeclException2(ExcInvalidOrder,
                 std::string,
                 unsigned int,
                 << "You tried to generate a " << arg1
                 << " object; no order is needed for objects of this kind, but "
                 << arg2 << " was given as argument.");
  /**
   * Exception
   */
  DeclException1(ExcInvalidQuadrature,
                 std::string,
                 << arg1 << " is not a valid name for a quadrature rule.");
  /** @} */
private:
  /**
   * This static function creates a quadrature object according to the name
   * given as a string, and the appropriate order (if the name is "gauss"). It
   * is called from the constructor.
   */
  static Quadrature<dim>
  create_quadrature(const std::string &s, const unsigned int order);
};
DEAL_II_NAMESPACE_CLOSE

#endif

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


#ifndef __deal2__quadrature_selector_h
#define __deal2__quadrature_selector_h


#include <deal.II/base/quadrature.h>
#include <deal.II/base/exceptions.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the quadrature rule passed to its constructor
 * as a string. Supported quadratures are QGauss (of all
 * orders), QMidpoint, QMilne, QSimpson,
 * QTrapez and QWeddle.
 *
 * This class is useful if you want to use flexible quadrature rules,
 * that are read from a parameter file (see ParameterHandler for
 * this).
 *
 * @ingroup Quadrature
 * @author Ralf Schulz, 2003
 */
template<int dim>
class QuadratureSelector : public Quadrature<dim>
{
public:
  /**
   * Constructor. Takes the name of
   * the quadrature rule (one of
   * "gauss", "milne", "weddle",
   * etc) and, if it iss "gauss",
   * the order of the quadrature
   * rule as argument.
   */
  QuadratureSelector (const std::string &s,
                      const unsigned int order=0);

  /**
   * This function returns all
   * possible names for quadratures
   * as a list separated by <tt>|</tt>,
   * so that you can use it for the
   * definition of parameter files
   * (see ParameterHandler for
   * details).
   */
  static std::string get_quadrature_names();

  /** @addtogroup Exceptions
   * @{ */


  /**
   * Exception
   */
  DeclException1 (ExcInvalidQGaussOrder,
                  int,
                  << "You tried to generate QGauss with an invalid order of "
                  << arg1 << " (must be >= 2)");
  /**
   * Exception
   */
  DeclException2 (ExcInvalidOrder,
                  std::string,
                  unsigned int,
                  << "You tried to generate a " << arg1
                  << " object; no order is needed (" << arg2
                  << " was given as parameter)");
  /**
   * Exception
   */
  DeclException1 (ExcInvalidQuadrature,
                  std::string,
                  << arg1
                  << " is not a valid quadrature name for a quadrature rule");
  //@}
private:
  /**
   * This static function creates a
   * quadrature object according to
   * the name given as a string,
   * and the appropriate order (if
   * the name is "gauss"). It is
   * called from the constructor.
   */
  static
  Quadrature<dim>
  create_quadrature (const std::string &s,
                     const unsigned int order);
};
DEAL_II_NAMESPACE_CLOSE

#endif

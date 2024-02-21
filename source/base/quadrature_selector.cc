// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
Quadrature<dim>
QuadratureSelector<dim>::create_quadrature(const std::string &s,
                                           const unsigned int order)
{
  if (s == "gauss")
    return QGauss<dim>(order);
  else if (s == "gauss_lobatto")
    return QGaussLobatto<dim>(order);
  else if (s == "gauss_chebyshev")
    return QGaussChebyshev<dim>(order);
  else if (s == "gauss_radau_chebyshev")
    return QGaussRadauChebyshev<dim>(order);
  else if (s == "gauss_lobatto_chebyshev")
    return QGaussLobattoChebyshev<dim>(order);
  else
    {
      AssertThrow(order == 0, ExcInvalidOrder(s, order));

      if (s == "midpoint")
        return QMidpoint<dim>();
      else if (s == "milne")
        return QMilne<dim>();
      else if (s == "simpson")
        return QSimpson<dim>();
      // The following name is DEAL_II_DEPRECATED and should be
      // removed once the 'using' declaration is removed.
      else if (s == "trapez")
        return QTrapezoid<dim>();
      else if (s == "trapezoid")
        return QTrapezoid<dim>();
      else if (s == "weddle")
        return QWeddle<dim>();
    }

  // we didn't find this name
  AssertThrow(false, ExcInvalidQuadrature(s));
  // return something to suppress
  // stupid warnings by some
  // compilers
  return Quadrature<dim>();
}



template <int dim>
QuadratureSelector<dim>::QuadratureSelector(const std::string &s,
                                            const unsigned int order)
  : Quadrature<dim>(create_quadrature(s, order).get_points(),
                    create_quadrature(s, order).get_weights())
{}



template <int dim>
std::string
QuadratureSelector<dim>::get_quadrature_names()
{
  return std::string(
    "gauss|gauss_lobatto|gauss_chebyshev|gauss_radau_chebyshev|gauss_lobatto_chebyshev|midpoint|milne|simpson|trapez|weddle");
}



// explicit instantiations
template class QuadratureSelector<1>;
template class QuadratureSelector<2>;
template class QuadratureSelector<3>;

DEAL_II_NAMESPACE_CLOSE

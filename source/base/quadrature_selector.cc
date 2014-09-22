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

#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
Quadrature<dim>
QuadratureSelector<dim>::
create_quadrature (const std::string &s,
                   const unsigned int order)
{
  if (s == "gauss")
    {
      AssertThrow(order >= 2, ExcInvalidQGaussOrder(order));
      return QGauss<dim>(order);
    }
  else
    {
      AssertThrow(order == 0, ExcInvalidOrder(s, order));

      if (s == "midpoint")        return QMidpoint<dim>();
      else if (s == "milne")      return QMilne<dim>();
      else if (s == "simpson")    return QSimpson<dim>();
      else if (s == "trapez")     return QTrapez<dim>();
      else if (s == "weddle")     return QWeddle<dim>();
    }

  // we didn't find this name
  AssertThrow (false, ExcInvalidQuadrature(s));
  // return something to suppress
  // stupid warnings by some
  // compilers
  return Quadrature<dim>();
}



template <int dim>
QuadratureSelector<dim>::QuadratureSelector (const std::string &s,
                                             const unsigned int order)
  :
  Quadrature<dim> (create_quadrature(s, order).get_points(),
                   create_quadrature(s, order).get_weights())
{
}



template <int dim>
std::string
QuadratureSelector<dim>::get_quadrature_names()
{
  return std::string("gauss|midpoint|milne|simpson|trapez|weddle");
}



// explicit instantiations
template class QuadratureSelector<1>;
template class QuadratureSelector<2>;
template class QuadratureSelector<3>;

DEAL_II_NAMESPACE_CLOSE

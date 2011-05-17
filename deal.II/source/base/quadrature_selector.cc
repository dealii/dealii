//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
Quadrature<dim>
QuadratureSelector<dim>::
create_quadrature (const std::string &s,
		   const unsigned int order)
{
  if(s == "gauss")
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

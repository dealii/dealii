//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_selector.h>
#include <base/quadrature_lib.h>


template <int dim>
Quadrature<dim>
QuadratureSelector<dim>::
create_quadrature (const std::string &s,
		   const unsigned int order)
{
  if(s == "gauss")
    {
      AssertThrow(order >= 2, ExcInvalidQGaussOrder(order));
      switch(order)
	{
	  case  2: return QGauss2<dim>();
	  case  3: return QGauss3<dim>();
	  case  4: return QGauss4<dim>();
	  case  5: return QGauss5<dim>();
	  case  6: return QGauss6<dim>();
	  case  7: return QGauss7<dim>();
	  default: return QGauss<dim>(order);
	}
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
  return QGauss2<dim>();
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

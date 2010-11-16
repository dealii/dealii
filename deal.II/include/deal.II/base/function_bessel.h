//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__function_bessel_h
#define __deal2__function_bessel_h


#include <base/config.h>
#include <base/function.h>
#include <base/point.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
/**
 * The Bessel functions of first kind or positive integer order.
 *
 * @author Guido Kanschat
 * @date 2010
 */
  template <int dim>
  class Bessel1 : public Function<dim>
  {
    public:
      Bessel1(const unsigned int order,
	      const double wave_number,
	      const Point<dim> center = Point<dim>());
      virtual double value (const Point<dim>& points, const unsigned int component) const;
      virtual void value_list (const std::vector<Point<dim> > &points,
			       std::vector<double>            &values,
			       const unsigned int              component = 0) const;
      virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				      const unsigned int  component = 0) const;
      virtual void gradient_list (const std::vector<Point<dim> > &points,
				  std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
    private:
      unsigned int order;
      double wave_number;
      Point<dim> center;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif


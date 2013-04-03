//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2010, 2011, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Test VectorFunctionTensorFunction

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>


template <int dim>
class X : public TensorFunction<1,dim>
{
public:
  virtual Tensor<1,dim> value (const Point<dim> &p) const
    {
      return p;
    }
};


template <int dim>
void check1 ()
{
  X<dim> x;
  VectorFunctionFromTensorFunction<dim>
    object (x, 1, dim+2);

  Assert (object.n_components == dim+2, ExcInternalError());

  for (unsigned int i=0; i<10; ++i)
    {
      Point<dim> p;
      for (unsigned int d=0; d<dim; ++d)
	p[d] = i+d;

      for (unsigned int c=0; c<dim+2; ++c)
	if (c==0 || c==dim+1)
	  Assert (object.value(p,c) == 0,
		ExcInternalError())
	else
	  Assert (object.value(p,c) == p[c-1],
		  ExcInternalError());

      Vector<double> v(dim+2);
      object.vector_value (p, v);
      for (unsigned int c=0; c<dim+2; ++c)
	if (c==0 || c==dim+1)
	  Assert (v(c) == 0,
		ExcInternalError())
	else
	  Assert (v(c) == p[c-1],
		  ExcInternalError());
    }

  deallog << "OK" << std::endl;
}




int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check1<1> ();
  check1<2> ();
  check1<3> ();
}



//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Generic routines to check consistency of function classes


// Check, whether the various implementations of function values are
// consistent. Arguments include the function, the amount of
// quadrature points in each direction, and the threshold above which
// we consider values unequal

template<int dim>
void
check_function_value_consistency(
  const Function<dim>& f,
  unsigned int sub,
  double threshold = 1.e-15)
{
  QMidpoint<1> mid;
  QIterated<dim> quadrature(mid, sub);

  std::vector<double> f1(quadrature.size());
  std::vector<Vector<double> > f2(quadrature.size(), Vector<double>(f.n_components));
  
  f.vector_value_list(quadrature.get_points(), f2);

  deallog << "value vs vector value list";
  for (unsigned int d=0;d<f.n_components;++d)
    for (unsigned int i=0;i<f1.size();++i)
      {
	const double v = f.value(quadrature.point(i), d);
	if (std::fabs(v-f2[i](d)) > threshold)
	  deallog << "v-vl " << d << ':' << i << ':' << v-f2[i](d);
      }
  deallog << std::endl << "value list vs vector value list";
  for (unsigned int d=0;d<f.n_components;++d)
    {
      f.value_list(quadrature.get_points(), f1, d);
      for (unsigned int i=0;i<f1.size();++i)
	{
	  if (std::fabs(f1[i]-f2[i](d)) > threshold)
	    deallog << ' ' << d << ':' << i << ':' << f1[i]-f2[i](d);
	}
    }
  deallog << std::endl;
}

// Same for gradients
template<int dim>
void
check_function_gradient_consistency(
  const Function<dim>& f,
  unsigned int sub,
  double threshold = 1.e-15)
{
  QMidpoint<1> mid;
  QIterated<dim> quadrature(mid, sub);

  std::vector<Tensor<1,dim> > f1(quadrature.size());
  std::vector<std::vector<Tensor<1,dim> > > f2(quadrature.size(),
					       std::vector<Tensor<1,dim> >(f.n_components));
  
  f.vector_gradient_list(quadrature.get_points(), f2);

  deallog << "gradient vs vector gradient list";
  for (unsigned int d=0;d<f.n_components;++d)
    for (unsigned int i=0;i<f1.size();++i)
      {
	const Tensor<1,dim> v = f.gradient(quadrature.point(i), d)-f2[i][d];
	
	if (std::sqrt(v*v) > threshold)
	  deallog << "v-vl " << d << ':' << i << ':' << v;
      }
  deallog << std::endl << "gradient list vs vector gradient list";
  for (unsigned int d=0;d<f.n_components;++d)
    {
      f.gradient_list(quadrature.get_points(), f1, d);
      for (unsigned int i=0;i<f1.size();++i)
	{
	  const Tensor<1,dim> v = f1[i]-f2[i][d];
	  if (std::sqrt(v*v) > threshold)
	    deallog << ' ' << d << ':' << i << ':' << v;
	}
    }
  deallog << std::endl;
}


//----------------------------  tensor_25.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2008, 2009, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_25.cc  ---------------------------

// check Tensor<2,dim>::operator[] (TableIndices)

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>


template <int dim>
void check ()
{
  typedef Tensor<3,dim> T;
  T t;
  for (unsigned int i=0; i<T::n_independent_components; ++i)
    t[T::unrolled_to_component_indices (i)] = (i+1)*(i+2);

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	{
	  Assert (TableIndices<3>(i,j,k)
		  ==
		  T::unrolled_to_component_indices (i*dim*dim+j*dim+k),
		  ExcInternalError());
	  Assert (T::component_to_unrolled_index(TableIndices<3>(i,j,k))
		  ==
		  i*dim*dim+j*dim+k,
		  ExcInternalError());
	  Assert (t[TableIndices<3>(i,j,k)] == t[i][j][k],
		  ExcInternalError());
	}
  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile("tensor_25/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();

  deallog << "OK" << std::endl;
}

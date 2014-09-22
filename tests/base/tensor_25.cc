// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();

  deallog << "OK" << std::endl;
}

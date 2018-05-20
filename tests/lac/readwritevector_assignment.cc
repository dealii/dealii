// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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

// Check constructor and operator=

#include "../tests.h"
#include <deal.II/base/index_set.h>
#include <deal.II/lac/read_write_vector.h>

void
test()
{
  unsigned int double_size = 2;
  unsigned int float_size  = 10;
  IndexSet     is(50);
  is.add_range(0, 2);
  is.add_index(46);
  is.add_range(10, 15);
  LinearAlgebra::ReadWriteVector<double> double_vector(is);
  LinearAlgebra::ReadWriteVector<float>  float_vector(float_size);
  deallog << "double_size " << double_vector.n_elements() << std::endl;
  deallog << "float_size " << float_vector.n_elements() << std::endl;

  double_vector = 0.;
  for(unsigned int i = 0; i < double_vector.n_elements(); ++i)
    double_vector.local_element(i) += i;
  for(unsigned int i = 0; i < float_vector.n_elements(); ++i)
    float_vector[i] = i;

  double_vector.print(deallog.get_file_stream());
  float_vector.print(deallog.get_file_stream());

  float_vector = double_vector;
  float_vector.print(deallog.get_file_stream());

  LinearAlgebra::ReadWriteVector<double> double_vector2(double_vector);
  double_vector2.print(deallog.get_file_stream());

  for(unsigned int i = 0; i < double_vector.n_elements(); ++i)
    double_vector2.local_element(i) += i;
  double_vector = double_vector2;
  double_vector.print(deallog.get_file_stream());
}

int
main()
{
  initlog();
  test();

  return 0;
}

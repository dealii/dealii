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

// Check reinit and add

#include "../tests.h"
#include <deal.II/base/index_set.h>
#include <deal.II/lac/read_write_vector.h>

#include <vector>

void
test()
{
  unsigned int size(10);
  std::vector<types::global_dof_index> indices(size);
  std::vector<float> float_values(size);
  std::vector<double> double_values(size);
  for (unsigned int i=0; i<size; ++i)
    {
      indices[i] = i;
      float_values[i] = 2.*i;
      double_values[i] = 3.*i;
    }

  LinearAlgebra::ReadWriteVector<float> float_vector;
  LinearAlgebra::ReadWriteVector<double> double_vector_1;
  LinearAlgebra::ReadWriteVector<double> double_vector_2;
  IndexSet is(50);
  is.add_range(0,10);
  is.add_index(46);
  is.add_range(11,25);
  is.compress();

  float_vector.reinit(size);
  double_vector_1.reinit(float_vector);
  double_vector_2.reinit(is);

  IndexSet is_copy(double_vector_2.get_stored_elements());
  is.print(std::cout);
  is_copy.print(std::cout);

  float_vector.add(indices,float_values);
  double_vector_1.add(indices,float_vector);
  double_vector_2.add(size,&indices[0],&double_values[0]);
  double_vector_1.swap(double_vector_2);
  LinearAlgebra::ReadWriteVector<double>::iterator val_1 = double_vector_1.begin();
  LinearAlgebra::ReadWriteVector<double>::iterator end_1 = double_vector_1.end();
  LinearAlgebra::ReadWriteVector<double>::iterator val_2 = double_vector_2.begin();
  LinearAlgebra::ReadWriteVector<double>::iterator end_2 = double_vector_2.end();
  deallog << "double_vector_1" << std::endl;
  for (; val_1<end_1; ++val_1)
    deallog << *val_1 <<std::endl;
  deallog << "double_vector_2" << std::endl;
  for (; val_2<end_2; ++val_2)
    deallog << *val_2 <<std::endl;
  double_vector_1.extract_subvector_to(indices,float_values);
  deallog << "subvector" << std::endl;
  for (unsigned int i=0; i<indices.size(); ++i)
    deallog << float_values[i] << std::endl;
}

int
main()
{
  initlog();
  test();

  return 0;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check that print_as_numpy_array works

#include <deal.II/lac/la_vector.h>

#include "../tests.h"


template <typename number>
void
print_vector()
{
  const unsigned int            size = 20;
  LinearAlgebra::Vector<number> vec(size);
  for (unsigned int i = 0; i < size; ++i)
    vec(i) = static_cast<number>(i);

  vec.print_as_numpy_array(deallog.get_file_stream());
}


int
main()
{
  initlog();

  print_vector<float>();
  print_vector<double>();
}

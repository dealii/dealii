// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// Check the functor interface


#include "../tests.h"
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/read_write_vector.templates.h>

#include <fstream>

struct Functor
{
  void operator() (double &value) const
  {
    value *= 2.;
  }
};

void test()
{
  const unsigned int size = 25;
  LinearAlgebra::ReadWriteVector<double> vector(size);
  for (unsigned int i=0; i<size; ++i)
    vector[i] = i;

  Functor functor;
  vector.apply(functor);
  for (unsigned int i=0; i<size; ++i)
    deallog << vector[i] << std::endl;
}

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  test();

  return 0;
}

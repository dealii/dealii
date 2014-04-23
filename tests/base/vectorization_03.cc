// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2012 - 2014 by the deal.II authors
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


// expression from step-48 that was computed in a wrong way for gcc-4.6.3 with
// vectorization


#include "../tests.h"

#include <deal.II/base/vectorization.h>


void test()
{
  VectorizedArray<double> array[3];
  for (unsigned int i=0; i<4; ++i)
    for (unsigned int v=0; v<VectorizedArray<double>::n_array_elements; ++v)
      array[i][v] = static_cast<double>(rand())/RAND_MAX;

  const VectorizedArray<double> tmp = 2. * array[0] - array[1] - array[2] * std::sin(array[0]);

  for (unsigned int v=0; v<VectorizedArray<double>::n_array_elements; ++v)
    deallog << tmp[v]-(2.*array[0][v]-array[1][v]-array[2][v]*std::sin(array[0][v])) << " ";
  deallog << std::endl;
}


int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();
}


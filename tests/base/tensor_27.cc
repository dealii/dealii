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


// Check the inverse of a rank-2 Tensor

// Equivalent matlab script:
/*
printf ("Tensor dim 1\n")
t1 = [2]
inv(t1)

printf ("Tensor dim 2\n")
t2 = [2 1; -1 1.5]
inv(t2)

printf ("Tensor dim 3\n")
t3 = [2 1 0.5; -1 1.5 0.25; 1.5 -0.75 1.25]
inv(t3)
*/

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(5);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  deallog << "Tensor dim 1" << std::endl;
  Tensor<2,1> t1;
  t1[0][0] = 2.0;
  deallog << invert(t1) << std::endl;
  Assert((invert(t1)*t1 - unit_symmetric_tensor<1>()).norm() < 1e-12,
         ExcMessage("Dim 1 inverse tensor definition is incorrect"));

  deallog << "Tensor dim 2" << std::endl;
  Tensor<2,2> t2;
  t2[0][0] = 2.0;
  t2[0][1] = 1.0;
  t2[1][0] = -1.0;
  t2[1][1] = 1.5;
  deallog << invert(t2) << std::endl;
  Assert((invert(t2)*t2 - unit_symmetric_tensor<2>()).norm() < 1e-12,
         ExcMessage("Dim 2 inverse tensor definition is incorrect"));

  deallog << "Tensor dim 3" << std::endl;
  Tensor<2,3> t3;
  t3[0][0] = 2.0;
  t3[0][1] = 1.0;
  t3[0][2] = 0.5;
  t3[1][0] = -1.0;
  t3[1][1] = 1.5;
  t3[1][2] = 0.25;
  t3[2][0] = 1.5;
  t3[2][1] = -0.75;
  t3[2][2] = 1.25;
  deallog << invert(t3) << std::endl;
  Assert((invert(t3)*t3 - unit_symmetric_tensor<3>()).norm() < 1e-12,
         ExcMessage("Dim 3 inverse tensor definition is incorrect"));

  deallog << "OK" << std::endl;
}

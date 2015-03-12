// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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


// check operator<< for SymmetricTensor<2,dim> and SymmetricTensor<4,dim>

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    const unsigned int dim=2;
    SymmetricTensor<2,dim> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[0][1] = 3;

    deallog << t << std::endl;
  }

  {
    const unsigned int dim=3;
    SymmetricTensor<4,dim> t;
    t[0][0][0][0] = t[1][0][1][0] = t[1][1][1][1]
                                    = t[2][2][2][2] = t[2][0][2][0] = 3;

    deallog << t << std::endl;
  }

  deallog << "OK" << std::endl;
}

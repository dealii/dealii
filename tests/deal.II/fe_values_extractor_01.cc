// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



// test that the FEValuesExtractors are copyable

#include "../tests.h"
#include <deal.II/fe/fe_values.h>

#include <fstream>




int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-8);

  {
    std::vector<FEValuesExtractors::Scalar> x;
    x.push_back (FEValuesExtractors::Scalar(42));
  }

  {
    std::vector<FEValuesExtractors::Vector> x;
    x.push_back (FEValuesExtractors::Vector(42));
  }

  deallog << "OK" << std::endl;
}

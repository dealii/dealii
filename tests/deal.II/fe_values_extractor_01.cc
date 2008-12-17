//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// test that the FEValuesExtractors are copyable

#include "../tests.h"
#include <fe/fe_values.h>

#include <fstream>




int main()
{
  std::ofstream logfile ("fe_values_extractor_01/output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-10);

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

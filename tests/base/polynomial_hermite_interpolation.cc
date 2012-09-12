//-----------------------------------------------------------------------------
//    $Id: polynomial_lagrange_product.cc 23959 2011-07-26 12:58:28Z kronbichler $
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomial.h>

using namespace Polynomials;

void plot(const std::vector<Polynomial<double> >& polynomials)
{
  LogStream::Prefix("plot");
  const unsigned int n=8;
  for (unsigned int i=0;i<=n;++i)
    {
      const double x = 1.*i/n;
      deallog << x;
      for (unsigned int p=0;p<polynomials.size();++p)
	deallog << '\t' << polynomials[p].value(x);
      deallog << std::endl;
    }
}

void interpolation_conditions(const std::vector<Polynomial<double> >& polynomials)
{
  std::vector<double> values(2);
  for (unsigned int i=0;i<polynomials.size();++i)
    {
      polynomials[i].value(0., values);
      deallog << i << "\t 0:\t" << values[0] << '\t' << values[1] << std::endl;
      polynomials[i].value(1., values);
      deallog << i << "\t 1:\t" << values[0] << '\t' << values[1] << std::endl;
    }
}


int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-12);
  
  interpolation_conditions(HermiteInterpolation::generate_complete_basis(6));
  plot(HermiteInterpolation::generate_complete_basis(6));
}

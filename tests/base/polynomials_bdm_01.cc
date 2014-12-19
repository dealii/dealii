// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// plot PolynomialsBDM on the reference cell

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/quadrature_lib.h>

#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

template <int dim>
void plot(const PolynomialsBDM<dim> &poly)
{
  QTrapez<1> base_quadrature;
  QIterated<dim> quadrature(base_quadrature, poly.degree()+3);
  std::vector<Tensor<1,dim> > values(poly.n());
  std::vector<Tensor<2,dim> > grads;
  std::vector<Tensor<3,dim> > grads2;


  for (unsigned int k=0; k<quadrature.size(); ++k)
    {
      if (k%(poly.degree()+4) == 0)
        deallog << "BDM" << poly.degree() << '<' << dim << '>' << std::endl;

      deallog << "BDM" << poly.degree() << '<' << dim << '>'
              << '\t' << quadrature.point(k);
      poly.compute(quadrature.point(k), values, grads, grads2);

      for (unsigned int i=0; i<poly.n(); ++i)
        for (unsigned int d=0; d<dim; ++d)
          deallog << '\t' << values[i][d];
      deallog << std::endl;
    }
}

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  PolynomialsBDM<2> p20(0);
  PolynomialsBDM<2> p21(1);
  PolynomialsBDM<2> p22(2);

  plot(p20);
  plot(p21);
  plot(p22);
}

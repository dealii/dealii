// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
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


// plot PolynomialsRaviartThomas on the reference cell

#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"

using namespace std;

template <int dim>
void
plot(const PolynomialsRaviartThomas<dim> &poly)
{
  QTrapezoid<1>               base_quadrature;
  QIterated<dim>              quadrature(base_quadrature, poly.degree() + 3);
  std::vector<Tensor<1, dim>> values(poly.n());
  std::vector<Tensor<2, dim>> grads;
  std::vector<Tensor<3, dim>> grads2;
  std::vector<Tensor<4, dim>> thirds;
  std::vector<Tensor<5, dim>> fourths;


  for (unsigned int k = 0; k < quadrature.size(); ++k)
    {
      if (k % (poly.degree() + 4) == 0)
        deallog << "RT" << poly.degree() << '<' << dim << '>' << std::endl;

      deallog << "RT" << poly.degree() << '<' << dim << '>' << '\t'
              << quadrature.point(k);
      poly.evaluate(
        quadrature.point(k), values, grads, grads2, thirds, fourths);

      for (unsigned int i = 0; i < poly.n(); ++i)
        for (unsigned int d = 0; d < dim; ++d)
          deallog << '\t' << values[i][d];
      deallog << std::endl;
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  PolynomialsRaviartThomas<2> p20(0);
  PolynomialsRaviartThomas<2> p21(1);
  PolynomialsRaviartThomas<2> p22(2);

  plot(p20);
  plot(p21);
  plot(p22);

  PolynomialsRaviartThomas<3> p30(0);
  PolynomialsRaviartThomas<3> p31(1);
  PolynomialsRaviartThomas<3> p32(2);

  plot(p30);
  plot(p31);
  plot(p32);
}

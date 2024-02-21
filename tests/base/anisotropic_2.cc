// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AnisotropicPolynomials


#include <deal.II/base/tensor_product_polynomials.h>

#include "../tests.h"


using namespace Polynomials;

using PolVector = std::vector<Polynomial<double>>;


void
print_2d(const AnisotropicPolynomials<2> &aniso)
{
  const unsigned int N = 10, M = 13;
  for (unsigned int i = 0; i <= N; ++i)
    {
      deallog << std::endl;
      for (unsigned int j = 0; j <= M; ++j)
        {
          deallog << 1. * i / N << ' ' << 1. * j / M << ' ';
          for (unsigned int k = 0; k < aniso.n(); ++k)
            deallog << aniso.compute_value(k, Point<2>(1. * i / N, 1. * j / M))
                    << ' ';
          deallog << std::endl;
          for (unsigned int k = 0; k < aniso.n(); ++k)
            deallog << aniso.compute_grad(k, Point<2>(1. * i / N, 1. * j / M))
                    << ' ';
          deallog << std::endl;
        }
    }
}



template <class Pol>
void
check_2d()
{
  // two checks with higher degree in
  // x or y direction
  {
    PolVector                 pols[2] = {Pol::generate_complete_basis(3),
                                         Pol::generate_complete_basis(1)};
    std::vector<PolVector>    p(&pols[0], &pols[2]);
    AnisotropicPolynomials<2> aniso(p);
    print_2d(aniso);
  }
  {
    PolVector                 pols[2] = {Pol::generate_complete_basis(2),
                                         Pol::generate_complete_basis(3)};
    std::vector<PolVector>    p(&pols[0], &pols[2]);
    AnisotropicPolynomials<2> aniso(p);

    print_2d(aniso);
  }
}


void
print_3d(const AnisotropicPolynomials<3> &aniso)
{
  const unsigned int N = 4, M = 3, P = 5;
  for (unsigned int i = 0; i <= N; ++i)
    {
      deallog << std::endl;
      for (unsigned int j = 0; j <= M; ++j)
        {
          deallog << std::endl;
          for (unsigned int k = 0; k <= P; ++k)
            {
              deallog << 1. * i / N << ' ' << 1. * j / M << ' ' << 1. * k / P;
              for (unsigned int k = 0; k < aniso.n(); ++k)
                deallog << aniso.compute_value(
                             k, Point<3>(1. * i / N, 1. * j / M, 1. * k / P))
                        << ' ';
              deallog << std::endl;
              for (unsigned int k = 0; k < aniso.n(); ++k)
                deallog << aniso.compute_grad(
                             k, Point<3>(1. * i / N, 1. * j / M, 1. * k / P))
                        << ' ';
              deallog << std::endl;
            }
        }
    }
}



template <class Pol>
void
check_3d()
{
  // three checks with higher degree
  // in x, y or z direction
  {
    PolVector                 pols[3] = {Pol::generate_complete_basis(3),
                                         Pol::generate_complete_basis(1),
                                         Pol::generate_complete_basis(1)};
    std::vector<PolVector>    p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso(p);
    print_3d(aniso);
  }
  {
    PolVector                 pols[3] = {Pol::generate_complete_basis(1),
                                         Pol::generate_complete_basis(3),
                                         Pol::generate_complete_basis(1)};
    std::vector<PolVector>    p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso(p);
    print_3d(aniso);
  }
  {
    PolVector                 pols[3] = {Pol::generate_complete_basis(1),
                                         Pol::generate_complete_basis(2),
                                         Pol::generate_complete_basis(3)};
    std::vector<PolVector>    p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso(p);
    print_3d(aniso);
  }
}



template <class Pol>
void
check()
{
  check_2d<Pol>();
  check_3d<Pol>();
}



int
main()
{
  initlog();

  deallog.push("Hierarchical");
  check<Hierarchical>();
  deallog.pop();
}

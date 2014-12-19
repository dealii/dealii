// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// check AnisotropicPolynomials


#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_product_polynomials.h>


using namespace Polynomials;

typedef std::vector<Polynomial<double> > PolVector;


void print_2d (const AnisotropicPolynomials<2> &aniso)
{
  const unsigned int N=10, M=13;
  for (unsigned int i=0; i<=N; ++i)
    {
      deallog << std::endl;
      for (unsigned int j=0; j<=M; ++j)
        {
          deallog << 1.*i/N << " "
                  << 1.*j/M << " ";
          for (unsigned int k=0; k<aniso.n(); ++k)
            deallog << aniso.compute_value (k, Point<2>(1.*i/N, 1.*j/M))
                    << " ";
          deallog << std::endl;
          for (unsigned int k=0; k<aniso.n(); ++k)
            deallog << aniso.compute_grad (k, Point<2>(1.*i/N, 1.*j/M))
                    << " ";
          deallog << std::endl;
        }
    }
}



template <class Pol>
void check_2d ()
{
  // two checks with higher degree in
  // x or y direction
  {
    PolVector pols[2] = { Pol::generate_complete_basis (3),
                          Pol::generate_complete_basis (1)
                        };
    std::vector<PolVector> p(&pols[0], &pols[2]);
    AnisotropicPolynomials<2> aniso (p);
    print_2d (aniso);
  }
  {
    PolVector pols[2] = { Pol::generate_complete_basis (2),
                          Pol::generate_complete_basis (3)
                        };
    std::vector<PolVector> p(&pols[0], &pols[2]);
    AnisotropicPolynomials<2> aniso (p);

    print_2d (aniso);
  }
}


void print_3d (const AnisotropicPolynomials<3> &aniso)
{
  const unsigned int N=4, M=3, P=5;
  for (unsigned int i=0; i<=N; ++i)
    {
      deallog << std::endl;
      for (unsigned int j=0; j<=M; ++j)
        {
          deallog << std::endl;
          for (unsigned int k=0; k<=P; ++k)
            {
              deallog << 1.*i/N << " "
                      << 1.*j/M << " "
                      << 1.*k/P;
              for (unsigned int k=0; k<aniso.n(); ++k)
                deallog << aniso.compute_value (k, Point<3>(1.*i/N, 1.*j/M, 1.*k/P))
                        << " ";
              deallog << std::endl;
              for (unsigned int k=0; k<aniso.n(); ++k)
                deallog << aniso.compute_grad (k, Point<3>(1.*i/N, 1.*j/M, 1.*k/P))
                        << " ";
              deallog << std::endl;
            }
        }
    }
}



template <class Pol>
void check_3d ()
{
  // three checks with higher degree
  // in x, y or z direction
  {
    PolVector pols[3] = { Pol::generate_complete_basis (3),
                          Pol::generate_complete_basis (1),
                          Pol::generate_complete_basis (1)
                        };
    std::vector<PolVector> p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso (p);
    print_3d (aniso);
  }
  {
    PolVector pols[3] = { Pol::generate_complete_basis (1),
                          Pol::generate_complete_basis (3),
                          Pol::generate_complete_basis (1)
                        };
    std::vector<PolVector> p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso (p);
    print_3d (aniso);
  }
  {
    PolVector pols[3] = { Pol::generate_complete_basis (1),
                          Pol::generate_complete_basis (2),
                          Pol::generate_complete_basis (3)
                        };
    std::vector<PolVector> p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso (p);
    print_3d (aniso);
  }
}



template <class Pol>
void check ()
{
  check_2d<Pol> ();
  check_3d<Pol> ();
}



int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

//   deallog.push("Lagrange");
//   check<LagrangeEquidistant> ();
//   deallog.pop();

//   deallog.push("Legendre");
//   check<Legendre> ();
//   deallog.pop();

  deallog.push("Hierarchical");
  check<Hierarchical> ();
  deallog.pop();
}

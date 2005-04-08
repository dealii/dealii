//----------------------------  anisotropic_2.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  anisotropic_2.cc  ---------------------------


// check AnisotropicPolynomials


#include "../tests.h"
#include <fstream>
#include <cmath>

#include <base/logstream.h>
#include <base/tensor_product_polynomials.h>


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
                          Pol::generate_complete_basis (1) };
    std::vector<PolVector> p(&pols[0], &pols[2]);
    AnisotropicPolynomials<2> aniso (p);
    print_2d (aniso);
  }
  {
    PolVector pols[2] = { Pol::generate_complete_basis (2),
                          Pol::generate_complete_basis (3) };
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
                          Pol::generate_complete_basis (1) };
    std::vector<PolVector> p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso (p);
    print_3d (aniso);
  }
  {
    PolVector pols[3] = { Pol::generate_complete_basis (1),
                          Pol::generate_complete_basis (3),
                          Pol::generate_complete_basis (1) };
    std::vector<PolVector> p(&pols[0], &pols[3]);
    AnisotropicPolynomials<3> aniso (p);
    print_3d (aniso);
  }
  {
    PolVector pols[3] = { Pol::generate_complete_basis (1),
                          Pol::generate_complete_basis (2),
                          Pol::generate_complete_basis (3) };
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
  std::ofstream logfile("anisotropic_2.output");
  logfile.precision(2);
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

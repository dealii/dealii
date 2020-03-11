// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// Test Chebyshev filter on Diagonal matrix with equidistance eigenvalues in
// (-1,1). Example is taken from Section 2.2. in Pieper et al, Journal of
// Computational Physics 325 (2016), 226-243
//
// The test exploits the fact that:
// x_0 =: a_i \psi_i
// H \psi_i = \lambda_i \psi_i
// p(H)x_0 = a_i p(\lambda_i) \psi_i

// Note that linfty_norm() norm will show non-zero difference for non-scaled
// case of high degree. This is ok, as the difference is small compared
// to the absolute value on the order 10^12

/* MWE in Maxima for degree 3 with scaling. The calculations are
 * done for the mode which has the maximum difference in the resulting value

Cheb2(n, x) :=
block ( [],
    if n = 0
       then 1
       else
        if n = 1
           then x
           else expand(2*x*Cheb2 (n - 1, x)
                        - Cheb2 (n - 2, x))
);
min:-0.01;
max:0.01;
c:(max+min)/2;
h:(max-min)/2;
L(x):= (x-c)/h;
ev: -0.0989011;
x:  0.904187;
aL: -0.1;
Cheb2(3,L(aL));
Cheb2(3,L(ev));
x * Cheb2(3,L(ev)) / Cheb2(3,L(aL));

// the last three produce:

-3970.0
-3839.905459408033
0.8745573293767684

 */

//#define EXTRA_OUTPUT


#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/utilities.h>

#include "../tests.h"

double
cheb2(const unsigned int d, const double x)
{
  if (d == 0)
    {
      return 1.;
    }
  else if (d == 1)
    {
      return x;
    }

  return 2. * x * cheb2(d - 1, x) - cheb2(d - 2, x);
}



void
check(const int          degree,
      const bool         scale = false,
      const double       a_L   = -0.1,
      const double       a     = -0.01,
      const double       b     = 0.01,
      const unsigned int size  = 1000)
{
  deallog << "Degree " << degree << std::endl;
  LinearAlgebra::distributed::Vector<double> ev(size), x(size), y(size),
    exact(size), diff(size);
  GrowingVectorMemory<LinearAlgebra::distributed::Vector<double>> vector_memory;

  for (unsigned int i = 0; i < size; ++i)
    ev(i) = -1. + 2. * (1. + i) / (1. + size);
  DiagonalMatrix<LinearAlgebra::distributed::Vector<double>> mat;
  mat.reinit(ev);

  x = 0.;
  // prevent overflow by not perturbing modes far away from the region
  // to be filtered
  unsigned int n_in  = 0;
  unsigned int n_out = 0;
  for (unsigned int i = 0; i < size; ++i)
    if (std::abs(ev(i)) <= std::abs(a_L))
      {
        if (ev(i) >= a && ev(i) <= b)
          n_in++;
        else
          n_out++;
        x(i) = random_value<double>();
      }

  deallog << " Modes inside/outside: " << n_in << " " << n_out << std::endl;

  // for x = x_i v_i , where v_i are eigenvectors
  // p[H]x = \sum_i x_i p(\lambda_i) v_i
  const double c = (a + b) / 2.;
  const double e = (b - a) / 2.;
  auto         L = [&](const double &x) { return (x - c) / e; };

  const double scaling = scale ? cheb2(degree, L(a_L)) : 1.; // p(L(a_L))
  deallog << " Scaling: " << scaling << " @ " << a_L << std::endl;
  exact = 0.;
  for (unsigned int i = 0; i < size; ++i)
    exact(i) = x(i) * cheb2(degree, L(ev(i))) / scaling;
  deallog << " Input norm: " << x.l2_norm() << std::endl;
  deallog << " Exact norm: " << exact.l2_norm() << std::endl;

  const double g_ = (scale ? a_L : std::numeric_limits<double>::infinity());
  y               = x;
  Utilities::LinearAlgebra::chebyshev_filter(
    y, mat, degree, std::make_pair(a, b), g_, vector_memory);
  diff = y;
  diff -= exact;

  deallog << " Filter [" << a << "," << b << "]" << std::endl;
  deallog << " Error: " << diff.linfty_norm() / exact.linfty_norm()
          << std::endl;

#ifdef EXTRA_OUTPUT
  // extra output for debugging:
  unsigned int max_i = 0;
  for (unsigned int i = 1; i < size; ++i)
    if (std::abs(diff(i)) > std::abs(diff(max_i)))
      max_i = i;

  deallog << " i =" << max_i << std::endl
          << " d =" << diff(max_i) << std::endl
          << " ev=" << ev(max_i) << std::endl
          << " x =" << x(max_i) << std::endl
          << " y =" << y(max_i) << std::endl
          << " ex=" << exact(max_i) << std::endl;
#endif
}


int
main()
{
  initlog();
  deallog << std::setprecision(6);

  deallog << "No scaling:" << std::endl;
  // no scaling:
  check(1);
  check(2);
  check(3);
  check(4);
  check(10);

  deallog << "Lower scaling:" << std::endl;
  // scaling at the lower end
  check(1, true);
  check(2, true);
  check(3, true);
  check(4, true);
  check(10, true);
  check(30, true);

  deallog << "Upper scaling:" << std::endl;
  // scaling at the upper end
  check(1, true, 0.1);
  check(2, true, 0.1);
  check(3, true, 0.1);
  check(4, true, 0.1);
  check(10, true, 0.1);
  check(30, true, 0.1);

  return 0;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


// test eigenvalue approximation by GMRES algorithm

#include "../tests.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>



template <typename number>
void test (unsigned int variant)
{
  const unsigned int n = variant % 2 == 0 ? 64 : 16;
  Vector<number> rhs(n), sol(n);
  rhs = 1.;

  LAPACKFullMatrix<number> matrix(n, n);

  // put diagonal entries of different strengths. these are very challenging
  // for GMRES and will usually take a lot of iterations until the Krylov
  // subspace is complete enough
  if (variant == 0)
    for (unsigned int i=0; i<n; ++i)
      matrix(i,i) = (i+1);
  else if (variant == 1)
    for (unsigned int i=0; i<n; ++i)
      matrix(i,i) = (i+1) * (i+1) * (i+1) * (i+1) * 1.001;
  else if (variant == 2)
    for (unsigned int i=0; i<n; ++i)
      matrix(i,i) = (i%2?1.:-1.)*(i+1);
  else if (variant == 3)
    for (unsigned int i=0; i<n; ++i)
      {
        matrix(i,i) = (i+1);
        if (i<n-1)
          matrix(i,i+1) = 1.5+i;
        if (i<n-2)
          matrix(i,i+2) = -1.65;
	matrix(i,n-1) = 2.;
	matrix(n-1,i) = -2.;
      }
  else
    Assert(false, ExcMessage("Invalid variant"));
  if (types_are_equal<number,float>::value == true)
    Assert(variant < 4, ExcMessage("Invalid_variant"));

  deallog.push(Utilities::int_to_string(variant,1));

  SolverControl control(1000, variant==1?1e-4:1e-13);
  typename SolverGMRES<Vector<number> >::AdditionalData data;
  data.max_n_tmp_vectors = 80;
  data.compute_eigenvalues = true;

  SolverGMRES<Vector<number> > solver(control, data);
  solver.solve(matrix, sol, rhs, PreconditionIdentity());

  if (variant == 0)
    {
      typename SolverCG<Vector<number> >::AdditionalData cg_data;
      cg_data.compute_eigenvalues = true;
      SolverCG<Vector<number> > solver_cg(control, cg_data);
      sol = 0;
      solver_cg.solve(matrix, sol, rhs, PreconditionIdentity());
    }

  if (variant == 3)
    {
      matrix.compute_eigenvalues();
      std::vector<std::complex<double> > eigenvalues(n);
      for (unsigned int i=0; i<n; ++i)
        eigenvalues[i] = matrix.eigenvalue(i);

      std::sort(eigenvalues.begin(), eigenvalues.end(),
                internal::SolverGMRES::complex_less_pred);

      deallog << "Actual eigenvalues:        ";
      for (unsigned int i=0; i<n; ++i)
        deallog << ' ' << eigenvalues[i];
      deallog << std::endl;
    }
  deallog.pop();
}

int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-8);

  deallog.push("double");
  test<double>(0);
  test<double>(1);
  test<double>(2);
  test<double>(3);
  deallog.pop();
}


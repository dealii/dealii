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



// we used to have constraint matrices for the FE_Q elements precomputed and
// stored, but later moved to compute them on the fly. make sure these
// correspond to the previously available precomputed ones for 3d and q=1,2


#include "../tests.h"
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>


// matrices taken from the old file deal.II/source/fe/fe_q_3d.cc
namespace FE_Q_3d
{
  static const double constraint_q1[] =
  {
    0.25,0.25,0.25,0.25,
    0.5,0,0.5,0,
    0,0.5,0,0.5,
    0.5,0.5,0,0,
    0,0,0.5,0.5
  };

  static const double constraint_q2[] =
  {
    0,0,0,0,0,0,0,0,1,
    0,0,0,0,1,0,0,0,0,
    0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,0,1,0,0,
    0,0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,0.375,-0.125,0.75,
    0,0,0,0,0,0,-0.125,0.375,0.75,
    0,0,0,0,0.375,-0.125,0,0,0.75,
    0,0,0,0,-0.125,0.375,0,0,0.75,
    0.375,0,-0.125,0,0.75,0,0,0,0,
    -0.125,0,0.375,0,0.75,0,0,0,0,
    0,0.375,0,-0.125,0,0.75,0,0,0,
    0,-0.125,0,0.375,0,0.75,0,0,0,
    0.375,-0.125,0,0,0,0,0.75,0,0,
    -0.125,0.375,0,0,0,0,0.75,0,0,
    0,0,0.375,-0.125,0,0,0,0.75,0,
    0,0,-0.125,0.375,0,0,0,0.75,0,
    0.140625,-0.046875,-0.046875,0.015625,0.28125,-0.09375,0.28125,-0.09375,0.5625,
    -0.046875,0.140625,0.015625,-0.046875,-0.09375,0.28125,0.28125,-0.09375,0.5625,
    -0.046875,0.015625,0.140625,-0.046875,0.28125,-0.09375,-0.09375,0.28125,0.5625,
    0.015625,-0.046875,-0.046875,0.140625,-0.09375,0.28125,-0.09375,0.28125,0.5625
  };
}


namespace Matrices
{
  const double *const
  constraint_matrices[] =
  {
    FE_Q_3d::constraint_q1,
    FE_Q_3d::constraint_q2
  };

  const unsigned int
  n_constraint_matrices
    = sizeof(constraint_matrices) /
      sizeof(constraint_matrices[0]);
}



void check ()
{
  // check for q=1,2
  for (unsigned int q=1; q<=2; ++q)
    {
      deallog << "q=" << q << std::endl;

      FE_Q<3> fe(q);

      FullMatrix<double> x(fe.constraints().m(),
                           fe.constraints().n());

      Assert (q<=Matrices::n_constraint_matrices,
              ExcInternalError());
      x.fill (Matrices::constraint_matrices[q-1]);

      for (unsigned int i=0; i<x.m(); ++i)
        for (unsigned int j=0; j<x.n(); ++j)
          {
            deallog << i << ' ' << j << ' '
                    << x(i,j) << ' '
                    << fe.constraints()(i,j)
                    << std::endl;
            Assert (std::fabs (x(i,j) - fe.constraints()(i,j))
                    <
                    1e-14,
                    ExcInternalError());
          }
    }
  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ();
}


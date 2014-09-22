// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// See documentation of ProductMatrix for documentation of this example

#include <deal.II/base/logstream.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

double Adata[] =
{
  .5, .1,
  .4, .2
};

double Bdata[] =
{
  .866, .5,
  -.5, .866
};


int main()
{
  FullMatrix<float> A(2,2);
  FullMatrix<double> B(2,2);

  A.fill(Adata);
  B.fill(Bdata);

  GrowingVectorMemory<Vector<double> > mem;

  ProductMatrix<Vector<double> > AB(A,B,mem);

  Vector<double> u(2);
  Vector<double> v(2);

  u(0) = 1.;
  u(1) = 2.;

  AB.vmult(v,u);

  deallog << v(0) << '\t' << v(1) << std::endl;

  AB.Tvmult(v,u);

  deallog << v(0) << '\t' << v(1) << std::endl;
}

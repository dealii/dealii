// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// Tests DiagonalMatrix


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/vector.h>



void
check()
{
  const unsigned int size = 10;
  Vector<double> vec(size);
  for (unsigned int i=0; i<size; ++i)
    vec(i) = 2;
  DiagonalMatrix<Vector<double> > mat;
  mat.reinit(vec);

  Vector<double> in (size), out(size), exact(size);
  for (unsigned int i=0; i<size; ++i)
    in(i) = 0.5;
  exact = 1;

  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult: " << out.linfty_norm() << std::endl;

  mat.Tvmult(out, in);
  out -= exact;
  deallog << "Error Tvmult: " << out.linfty_norm() << std::endl;

  out = exact;
  mat.vmult_add(out, in);
  out -= exact;
  out -= exact;
  deallog << "Error vmult_add: " << out.linfty_norm() << std::endl;

  out = exact;
  mat.Tvmult_add(out, in);
  out -= exact;
  out -= exact;
  deallog << "Error Tvmult_add: " << out.linfty_norm() << std::endl;

  for (unsigned int i=0; i<size; ++i)
    {
      mat(i,i) = i+1;
      in(i) = 1./mat(i,i);
    }

  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult set 1: " << out.linfty_norm() << std::endl;

  for (unsigned int i=0; i<size; ++i)
    {
      mat.get_vector()(i) = 1./(i+1);
      in(i) = i+1;
    }
  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult set 2: " << out.linfty_norm() << std::endl;

  for (unsigned int i=0; i<size; ++i)
    {
      const double mat_entry =
        const_cast<const DiagonalMatrix<Vector<double> > &>(mat)(i,i);
      mat(i,i) = in(i);
      in(i) = mat_entry;
    }
  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult set 3: " << out.linfty_norm() << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  check();

  return 0;
}

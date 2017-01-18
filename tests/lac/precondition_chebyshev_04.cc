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


// Similar to precondition_chebyshev_02 but checking that the Chebyshev
// initialization also works when the first vmult operation is done on a zero
// vector


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cmath>


class DiagonalMatrixManual
{
public:
  DiagonalMatrixManual() {}

  void set_vector_one (const unsigned int size)
  {
    diagonal.reinit(size);
    diagonal = 1;
  }

  void reinit (const FullMatrix<double> f)
  {
    diagonal.reinit(f.m());
    for (unsigned int i=0; i<f.m(); ++i)
      diagonal(i) = 1./f(i,i);
    diagonal.print(deallog);
  }

  typename Vector<double>::size_type m() const
  {
    return diagonal.size();
  }

  void vmult(Vector<double> &dst, const Vector<double> &src) const
  {
    dst = src;
    dst.scale(diagonal);
  }

private:
  Vector<double> diagonal;
};


void
check()
{
  const unsigned int size = 10;
  FullMatrix<double> m(size,size);
  for (unsigned int i=0; i<size; ++i)
    m(i,i) = i+1;

  Vector<double> in (size), out(size), ref(size), zero(size);
  for (unsigned int i=0; i<size; ++i)
    in(i) = (double)Testing::rand()/RAND_MAX;

  PreconditionChebyshev<FullMatrix<double>,Vector<double>,DiagonalMatrixManual> prec;
  PreconditionChebyshev<FullMatrix<double>,Vector<double>,DiagonalMatrixManual>::AdditionalData
  data;
  data.smoothing_range = 2 * size;
  data.degree = 3;
  data.preconditioner.reset(new DiagonalMatrixManual());
  data.preconditioner->set_vector_one(size);
  prec.initialize(m, data);

  deallog << "Exact inverse:     ";
  for (unsigned int i=0; i<size; ++i)
    deallog << in(i)/m(i,i) << " ";
  deallog << std::endl;

  deallog << "Check  vmult usual init: ";
  prec.vmult(ref, in);
  for (unsigned int i=0; i<size; ++i)
    deallog << ref(i) << " ";
  deallog << std::endl;

  prec.initialize(m, data);

  prec.vmult(out, zero);
  prec.vmult(out, in);
  deallog << "Check  vmult zero init:  ";
  for (unsigned int i=0; i<size; ++i)
    deallog << out(i) << " ";
  out -= ref;
  deallog << " difference norm = " << out.linfty_norm() << std::endl;

  PreconditionChebyshev<FullMatrix<double>,Vector<double>,DiagonalMatrixManual> prec2;
  prec2.initialize(m, data);
  prec2.vmult(out, zero);
  prec2.vmult(out, in);
  deallog << "Check  vmult zero init:  ";
  for (unsigned int i=0; i<size; ++i)
    deallog << out(i) << " ";
  out -= ref;
  deallog << " difference norm = " << out.linfty_norm() << std::endl;

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

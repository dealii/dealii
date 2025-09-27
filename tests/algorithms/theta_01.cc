// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// See documentation of ThetaTimestepping for documentation of this example

#include <deal.II/algorithms/operator.h>
#include <deal.II/algorithms/theta_timestepping.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

using namespace Algorithms;


class Explicit : public OperatorBase
{
public:
  Explicit(const FullMatrix<double> &matrix);
  void
  operator()(AnyData &out, const AnyData &in);

private:
  ObserverPointer<const FullMatrix<double>, Explicit> matrix;
  FullMatrix<double>                                  m;
};


class Implicit : public OperatorBase
{
public:
  Implicit(const FullMatrix<double> &matrix);
  void
  operator()(AnyData &out, const AnyData &in);

private:
  ObserverPointer<const FullMatrix<double>, Implicit> matrix;
  FullMatrix<double>                                  m;
};

// End of declarations

int
main()
{
  initlog();

  FullMatrix<double> matrix(2);
  matrix(0, 0) = 0.;
  matrix(1, 1) = 0.;
  matrix(0, 1) = numbers::PI;
  matrix(1, 0) = -numbers::PI;

  OutputOperator<Vector<double>> out;
  out.initialize_stream(deallog.get_file_stream());

  Explicit                          op_explicit(matrix);
  Implicit                          op_implicit(matrix);
  ThetaTimestepping<Vector<double>> solver(op_explicit, op_implicit);
  solver.timestep_control().start_step(0.1);
  // solver.set_output(out);

  Vector<double> value(2);
  value(0) = 1.;
  AnyData         indata;
  AnyData         outdata;
  Vector<double> *p = &value;
  outdata.add(p, "value");
  deallog << "Initial: " << value(0) << ' ' << value(1) << std::endl;
  solver.notify(Events::initial);
  solver(outdata, indata);
  deallog << "Result: " << value(0) << ' ' << value(1) << " Norm "
          << value.l2_norm() << std::endl;
}


Explicit::Explicit(const FullMatrix<double> &M)
  : matrix(&M)
{
  m.reinit(M.m(), M.n());
}


void
Explicit::operator()(AnyData &out, const AnyData &in)
{
  const double *step = in.read_ptr<double>("Timestep");

  if (this->notifications.test(Events::initial) ||
      this->notifications.test(Events::new_timestep_size))
    {
      m.equ(-*step, *matrix);
      for (unsigned int i = 0; i < m.m(); ++i)
        m(i, i) += 1.;
    }
  this->notifications.clear();
  unsigned int i = in.find("Previous iterate");
  m.vmult(*out.entry<Vector<double> *>(0), *in.read_ptr<Vector<double>>(i));
}


Implicit::Implicit(const FullMatrix<double> &M)
  : matrix(&M)
{
  m.reinit(M.m(), M.n());
}


void
Implicit::operator()(AnyData &out, const AnyData &in)
{
  const double *step = in.read_ptr<double>("Timestep");

  if (this->notifications.test(Events::initial) ||
      this->notifications.test(Events::new_timestep_size))
    {
      m.equ(*step, *matrix);
      for (unsigned int i = 0; i < m.m(); ++i)
        m(i, i) += 1.;
      m.gauss_jordan();
    }
  this->notifications.clear();

  unsigned int i = in.find("Previous time");
  m.vmult(*out.entry<Vector<double> *>(0), *in.read_ptr<Vector<double>>(i));
}

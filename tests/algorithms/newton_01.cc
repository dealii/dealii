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


#include <deal.II/algorithms/newton.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// test computing square root of 2 with newton's method


class SquareRoot : public EnableObserverPointer
{
public:
  void
  start_vector(Vector<double> &start) const;
  void
  residual(AnyData &out, const AnyData &in);
  void
  solve(AnyData &out, const AnyData &in);
};

class SquareRootResidual : public Algorithms::OperatorBase
{
  ObserverPointer<SquareRoot, SquareRootResidual> discretization;

public:
  SquareRootResidual(SquareRoot &problem)
    : discretization(&problem)
  {}

  virtual void
  operator()(AnyData &out, const AnyData &in)
  {
    discretization->residual(out, in);
  }
};

class SquareRootSolver : public Algorithms::OperatorBase
{
  ObserverPointer<SquareRoot, SquareRootSolver> solver;

public:
  SquareRootSolver(SquareRoot &problem)
    : solver(&problem)
  {}

  virtual void
  operator()(AnyData &out, const AnyData &in)
  {
    solver->solve(out, in);
  }
};

void
SquareRoot::residual(AnyData &out, const AnyData &in)
{
  Vector<double> &v = *out.entry<Vector<double> *>(0);
  // residuum = 0
  v(0)                    = 0.;
  const Vector<double> &x = *in.entry<const Vector<double> *>("Newton iterate");
  v(0)                    = x * x - 2.;
}

void
SquareRoot::solve(AnyData &out, const AnyData &in)
{
  Vector<double> &v       = *out.entry<Vector<double> *>(0);
  v(0)                    = 0;
  const Vector<double> &x = *in.entry<const Vector<double> *>("Newton iterate");
  const Vector<double> &r =
    *in.entry<const Vector<double> *>("Newton residual");

  v(0) = 1. / 2. / x(0) * r(0);
}



void
test()
{
  SquareRoot         square_root;
  SquareRootSolver   sq_solver(square_root);
  SquareRootResidual sq_residual(square_root);

  Algorithms::OutputOperator<Vector<double>> output;

  Algorithms::Newton<Vector<double>> newton(sq_residual, sq_solver);
  newton.initialize(output);

  AnyData in_data;
  AnyData out_data;

  Vector<double> solution(1);
  solution = 10.;
  out_data.add<Vector<double> *>(&solution, "solution");

  newton.control.set_reduction(1.e-20);
  newton.control.log_history(true);
  newton.debug_vectors = true;
  newton(out_data, in_data);
  deallog << " square root " << (*out_data.read<Vector<double> *>(0))(0)
          << std::endl;
}



int
main()
{
  initlog();

  test();
}

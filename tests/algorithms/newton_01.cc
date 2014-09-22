// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/algorithms/operator.h>
#include <deal.II/algorithms/newton.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/named_data.h>


//test computing square root of 2 with newton's method

using namespace dealii;

class SquareRoot : public Subscriptor
{
public:
  void start_vector (Vector<double> &start) const;
  void residual (AnyData &out, const AnyData &in);
  void solve (AnyData &out, const AnyData &in);
};

class SquareRootResidual : public
  Algorithms::Operator<Vector<double> >
{
  SmartPointer<SquareRoot, SquareRootResidual>
  discretization;
public:

  SquareRootResidual(SquareRoot &problem)
    : discretization(&problem)
  {}

  virtual void operator ()(AnyData &out, const AnyData &in)
  {
    discretization->residual(out,in);
  }
};

class SquareRootSolver : public
  Algorithms::Operator<Vector<double> >
{
  SmartPointer<SquareRoot, SquareRootSolver>
  solver;
public:

  SquareRootSolver(SquareRoot &problem)
    : solver(&problem)
  {}

  virtual void operator ()(AnyData &out, const AnyData &in)
  {
    solver->solve(out,in);
  }
};

void SquareRoot::residual (AnyData &out, const AnyData &in)
{
  Vector<double>& v = *out.entry<Vector<double>*>(0);
  //residuum = 0
  v(0) = 0.;
  const Vector<double> &x = *in.entry<const Vector<double>*>("Newton iterate");
  v(0) = x*x - 2.;
}

void SquareRoot::solve (AnyData &out, const AnyData &in)
{
  Vector<double>& v = *out.entry<Vector<double>*>(0);
  v(0) = 0;
  const Vector<double> &x = *in.entry<const Vector<double>*>("Newton iterate");
  const Vector<double> &r = *in.entry<const Vector<double>*>("Newton residual");

  v(0) = 1./2./x(0)* r(0);
}



void test ()
{
  SquareRoot square_root;
  SquareRootSolver sq_solver (square_root);
  SquareRootResidual sq_residual (square_root);

  Algorithms::OutputOperator<Vector<double> > output;

  Algorithms::Newton<Vector<double> > newton(
    sq_residual,
    sq_solver);
  newton.initialize(output);

  NamedData<Vector<double>*> in_data;
  NamedData<Vector<double>*> out_data;

  Vector<double> solution (1);
  solution = 10.;
  Vector<double> *p = &solution;
  out_data.add(p, "solution");

  newton.control.set_reduction(1.e-20);
  newton.control.log_history(true);
  newton.debug_vectors = true;
  newton(out_data, in_data);
  deallog << " square root " << (*out_data(0))(0) << std::endl;
}




int main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}

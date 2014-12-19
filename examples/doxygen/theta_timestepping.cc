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


// See documentation of ThetaTimestepping for documentation of this example

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/algorithms/operator.h>
#include <deal.II/algorithms/theta_timestepping.h>

#include <iostream>

using namespace dealii;
using namespace Algorithms;


class Explicit : public Operator<Vector<double> >
{
public:
  Explicit(const FullMatrix<double> &matrix);
  void operator() (AnyData &out, const AnyData &in);

private:
  SmartPointer<const FullMatrix<double>, Explicit> matrix;
  FullMatrix<double> m;
};


class Implicit : public Operator<Vector<double> >
{
public:
  Implicit(const FullMatrix<double> &matrix);
  void operator() (AnyData &out, const AnyData &in);

private:
  SmartPointer<const FullMatrix<double>, Implicit> matrix;
  FullMatrix<double> m;
};

// End of declarations

int main()
{
  FullMatrix<double> matrix(2);
  matrix(0,0) = 0.;
  matrix(1,1) = 0.;
  matrix(0,1) = 3.14;
  matrix(1,0) = -3.14;

  OutputOperator<Vector<double> > out;
  out.initialize_stream(std::cout);

  Explicit op_explicit(matrix);
  Implicit op_implicit(matrix);
  ThetaTimestepping<Vector<double> > solver(op_explicit, op_implicit);
  solver.set_output(out);

  Vector<double> value(2);
  value(0) = 1.;
  AnyData indata;
  AnyData outdata;
  outdata.add(&value, "value");

  solver.notify(Events::initial);
  solver(outdata, indata);
}


Explicit::Explicit(const FullMatrix<double> &M)
  :
  matrix(&M)
{
  m.reinit(M.m(), M.n());
}


void
Explicit::operator() (AnyData &out, const AnyData &in)
{
  const double timestep = *in.read_ptr<double>("Timestep");
  if (this->notifications.test(Events::initial) || this->notifications.test(Events::new_timestep_size))
    {
      m.equ(-timestep, *matrix);
      for (unsigned int i=0; i<m.m(); ++i)
        m(i,i) += 1.;
    }
  this->notifications.clear();
  m.vmult(*out.entry<Vector<double>*>(0),
          *in.read_ptr<Vector<double> >("Previous iterate"));
}


Implicit::Implicit(const FullMatrix<double> &M)
  :
  matrix(&M)
{
  m.reinit(M.m(), M.n());
}


void
Implicit::operator() (AnyData &out, const AnyData &in)
{
  const double timestep = *in.read_ptr<double>("Timestep");
  if (this->notifications.test(Events::initial) || this->notifications.test(Events::new_timestep_size))
    {
      m.equ(timestep, *matrix);
      for (unsigned int i=0; i<m.m(); ++i)
        m(i,i) += 1.;
      m.gauss_jordan();
    }
  this->notifications.clear();
  m.vmult(*out.entry<Vector<double>*>(0),
          *in.read_ptr<Vector<double> >("Previous time"));
}



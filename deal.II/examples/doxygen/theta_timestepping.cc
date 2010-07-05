//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

// See documentation of ThetaTimestepping for documentation of this example

#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>

#include <numerics/operator.h>
#include <numerics/theta_timestepping.h>

using namespace dealii;
using namespace Algorithms;


class Explicit
  : public Operator<Vector<double> >
{
  public:
    Explicit(const FullMatrix<double>& matrix);
    void operator() (NamedData<Vector<double>*>& out,
		     const NamedData<Vector<double>*>& in);

    void initialize_timestep_data(const TimestepData&);
  private:
    const TimestepData* timestep_data;
    SmartPointer<const FullMatrix<double>, Explicit> matrix;
    FullMatrix<double> m;
};

  
class Implicit
  : public Operator<Vector<double> >
{
  public:
    Implicit(const FullMatrix<double>& matrix);
    void operator() (NamedData<Vector<double>*>& out,
		     const NamedData<Vector<double>*>& in);

    void initialize_timestep_data(const TimestepData&);
  private:
    const TimestepData* timestep_data;
    SmartPointer<const FullMatrix<double>, Implicit> matrix;
    FullMatrix<double> m;
};

  
int main()
{
  FullMatrix<double> matrix(2);
  matrix(0,0) = -.1;
  matrix(1,1) = -.1;
  matrix(0,1) = 31.4;
  matrix(1,0) = -31.4;

  Explicit op_explicit(matrix);
  Implicit op_implicit(matrix);
  ThetaTimestepping<Vector<double> > solver(op_explicit, op_implicit);
  op_explicit.initialize_timestep_data(solver.explicit_data());
  op_implicit.initialize_timestep_data(solver.implicit_data());
  solver.notify(Events::initial);
  
  Vector<double> value(2);
  NamedData<Vector<double>*> indata;
  NamedData<Vector<double>*> outdata;
  Vector<double>* p = &value;
  outdata.add(p, "value");
  solver(outdata, indata);
}


Explicit::Explicit(const FullMatrix<double>& M)
		:
		matrix(&M)
{}


void
Explicit::initialize_timestep_data(const TimestepData& t)
{
  timestep_data = &t;
}


void
Explicit::operator() (NamedData<Vector<double>*>& out, const NamedData<Vector<double>*>& in)
{
  this->notifications.print(deallog);
  deallog << std::endl;
  unsigned int i = in.find("Previous time");
  m.vmult(*out(0), *in(i));
}


Implicit::Implicit(const FullMatrix<double>& M)
		:
		matrix(&M)
{}


void
Implicit::initialize_timestep_data(const TimestepData& t)
{
  timestep_data = &t;
}


void
Implicit::operator() (NamedData<Vector<double>*>& out, const NamedData<Vector<double>*>& in)
{
  this->notifications.print(deallog);
  deallog << std::endl;
  unsigned int i = in.find("Previous time");
  m.vmult(*out(0), *in(i));
}



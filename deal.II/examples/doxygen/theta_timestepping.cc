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

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/operator.h>
#include <deal.II/numerics/theta_timestepping.h>
#include <deal.II/numerics/operator.h>

#include <iostream>

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

// End of declarations

int main()
{
  FullMatrix<double> matrix(2);
  matrix(0,0) = 1.;
  matrix(1,1) = 1.;
  matrix(0,1) = 31.4;
  matrix(1,0) = -31.4;

  OutputOperator<Vector<double> > out;
  out.initialize_stream(std::cout);
  
  Explicit op_explicit(matrix);
  Implicit op_implicit(matrix);
  ThetaTimestepping<Vector<double> > solver(op_explicit, op_implicit);
  op_explicit.initialize_timestep_data(solver.explicit_data());
  op_implicit.initialize_timestep_data(solver.implicit_data());
  solver.set_output(out);
  
  Vector<double> value(2);
  value(0) = 1.;
  NamedData<Vector<double>*> indata;
  NamedData<Vector<double>*> outdata;
  Vector<double>* p = &value;
  outdata.add(p, "value");
  
  solver.notify(Events::initial);
  solver(outdata, indata);
}


Explicit::Explicit(const FullMatrix<double>& M)
		:
		matrix(&M)
{
  m.reinit(M.m(), M.n());
}


void
Explicit::initialize_timestep_data(const TimestepData& t)
{
  timestep_data = &t;
}


void
Explicit::operator() (NamedData<Vector<double>*>& out, const NamedData<Vector<double>*>& in)
{
  if (this->notifications.test(Events::initial) || this->notifications.test(Events::new_timestep_size))
    {
      m.equ(-timestep_data->step, *matrix);
      for (unsigned int i=0;i<m.m();++i)
	m(i,i) += 1.;
    }
  this->notifications.clear();
  unsigned int i = in.find("Previous iterate");
  m.vmult(*out(0), *in(i));
}


Implicit::Implicit(const FullMatrix<double>& M)
		:
		matrix(&M)
{
  m.reinit(M.m(), M.n());
}


void
Implicit::initialize_timestep_data(const TimestepData& t)
{
  timestep_data = &t;
}


void
Implicit::operator() (NamedData<Vector<double>*>& out, const NamedData<Vector<double>*>& in)
{
  if (this->notifications.test(Events::initial) || this->notifications.test(Events::new_timestep_size))
    {
      m.equ(timestep_data->step, *matrix);
      for (unsigned int i=0;i<m.m();++i)
	m(i,i) += 1.;
      m.gauss_jordan();
    }
  this->notifications.clear();
  
  unsigned int i = in.find("Previous time");
  m.vmult(*out(0), *in(i));
}



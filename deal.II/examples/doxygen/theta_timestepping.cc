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
  : Operator<Vector<double> >
{
  public:
    void operator() (NamedData<Vector<double>*>& out,
		     const NamedData<Vector<double>*>& in);

    void initialize_timestep_data(const TimeStepData&);
  private:
    const TimeStepData* timestep_data;
};

  
class Implicit
  : Operator<Vector<double> >
{
  public:
    void operator() (NamedData<Vector<double>*>& out,
		     const NamedData<Vector<double>*>& in);

    void initialize_timestep_data(const TimeStepData&);
  private:
    const TimeStepData* timestep_data;
};

  
int main()
{
  Explicit op_explicit;
  Implicit op_implicit;
  ThetaTimestepping<Vector<double> > solver(op_explicit, op_implicit);
  op_explicit.initialize_timestep_data(solver.explicit_data());
  op_implicit.initialize_timestep_data(solver.implicit_data());

  Vector<double> value(2);
  NamedData<Vector<double>*> indata;
  NamedData<Vector<double>*> outdata;
  outdata.add(&value, "value");
  solver(outdata, indata);
}

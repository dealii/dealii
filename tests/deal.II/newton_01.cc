//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/numerics/operator.h>
#include <deal.II/numerics/newton.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/base/named_data.h>


//test computing square root of 2 with newton's method

using namespace dealii;

class SquareRoot : public Subscriptor
{
  public:
    void start_vector (Vector<double>& start) const;
    void residual (NamedData<Vector<double>*>& out, 
        const NamedData<Vector<double>*>& in);
    void solve (NamedData<Vector<double>*>& out, 
        const NamedData<Vector<double>*>& in);
};

class SquareRootResidual : public 
                       Algorithms::Operator<Vector<double> >
{
    SmartPointer<SquareRoot, SquareRootResidual>
      discretization;
  public:

    SquareRootResidual(SquareRoot& problem)
      : discretization(&problem)
    {}

    virtual void operator ()(
        NamedData<Vector<double>*>& out, 
        const NamedData<Vector<double>*>& in) 
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

    SquareRootSolver(SquareRoot& problem)
      : solver(&problem)
    {}

    virtual void operator ()(
        NamedData<Vector<double>*>& out, 
        const NamedData<Vector<double>*>& in) 
    {
      solver->solve(out,in);
    }
};

void SquareRoot::residual (
    NamedData<Vector<double>*>& out, 
    const NamedData<Vector<double>*>& in) 
{
  //residuum = 0
  *out(0) = 0.;
  const Vector<double> &x = *in(in.find("Newton iterate")); 
  *out(0) = x*x - 2.;
}

void SquareRoot::solve (
    NamedData<Vector<double>*>& out, 
    const NamedData<Vector<double>*>& in)
{
  *out(0) = 0;
  const Vector<double> &x = *in(in.find("Newton iterate")); 
  const Vector<double> &r = *in(in.find("Newton residual"));

  (*out(0))(0) = 1./2./x(0)* r(0);
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
  Vector<double>* p = &solution;
  out_data.add(p, "solution");

  newton.control.set_reduction(1.e-20);
  newton.control.log_history(true);
  newton.debug_vectors = true;
  newton(out_data, in_data);
  deallog << " square root " << (*out_data(0))(0) << std::endl;
}

  
  

int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}

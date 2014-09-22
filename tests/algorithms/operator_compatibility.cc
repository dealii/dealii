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

// Test the compatibility functionality for Operator::operator()

#include "../tests.h"
#include <deal.II/algorithms/operator.h>
#include <deal.II/lac/vector.h>

using namespace Algorithms;

class Op1 : public Operator<Vector<double> >
{
  public:
  virtual void operator() (AnyData& out, const AnyData& in)
  {
    Vector<double>* u = out.entry<Vector<double>*>("u");
    deallog << "u " << (*u)(0) << std::endl;
    
    const Vector<double>* v = in.entry<Vector<double>*>("v");
    deallog << "v " << (*v)(0) << std::endl;
  }
};

class Op2 : public Operator<Vector<double> >
{
  public:
  virtual void operator() (NamedData<Vector<double>*>& out,
			   const NamedData<Vector<double>*>& in)
  {
    Vector<double>* u = out(out.find("u"));
    deallog << "u " << (*u)(0) << std::endl;
    
    const Vector<double> * const v = in(in.find("v"));
    deallog << "v " << (*v)(0) << std::endl;
  }
};

void test(Operator<Vector<double> >& op)
{
  Vector<double> u(1);
  Vector<double>* p = &u;
  Vector<double> v(1);
  u(0) = 3.;
  v(0) = 7.;
  
  AnyData o1;
  o1.add (p, "u");
  AnyData i1;
  i1.add (&v, "v");
  
  NamedData<Vector<double>*> o2;
  o2.add (p, "u");
  NamedData<Vector<double>*> i2;
  i2.add (&v, "v");
  
  deallog << "Call with AnyData" << std::endl;
  op(o1, i1);
  deallog << "Call with NamedData" << std::endl;
  op(o2, i2);
}

int main()
{
  //  deal_II_exceptions::disable_abort_on_exception();
  
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  deallog << "Operator with AnyData" << std::endl;
  Op1 op1;
  test(op1);
  
  deallog << "Operator with NamedData" << std::endl;
  Op2 op2;
  test(op2);

  // deallog << "Operator with NamedData" << std::endl;
  // Op2 op2;
  // test(op1);
}

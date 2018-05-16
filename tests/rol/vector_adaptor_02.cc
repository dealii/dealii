// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#include <cmath>
#include <iostream>
#include <sstream>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/optimization/rol/vector_adaptor.h>

#include "ROL_Objective.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Use ROL to minimize the objective function, f(x,y) = x^2 + y^2.

using namespace dealii;

using VectorType = typename dealii::Vector<double>;

template<class Real=double, typename Xprim=Rol::VectorAdaptor<VectorType> >
class QuadraticObjective : public ROL::Objective<Real>
{

private:

  Teuchos::RCP<const VectorType>
  get_rcp_to_VectorType (const ROL::Vector<Real> &x)
  {
    return (Teuchos::dyn_cast<const Xprim>(x)).getVector();
  }

  Teuchos::RCP<dealii::Vector<Real> >
  get_rcp_to_VectorType (ROL::Vector<Real> &x)
  {
    return (Teuchos::dyn_cast<Xprim>(x)).getVector();
  }

public:

  Real
  value (const ROL::Vector<Real> &x,
         Real                    &tol)
  {
    Assert (x.dimension()==2,
            ExcInternalError());

    return x.dot(x);
  }

  void
  gradient (ROL::Vector<Real>       &g,
            const ROL::Vector<Real> &x,
            Real                    &tol)
  {
    Teuchos::RCP<const VectorType> xp = this->get_rcp_to_VectorType(x);
    Teuchos::RCP<VectorType> gp = this->get_rcp_to_VectorType(g);

    (*gp)[0] = 2. * (*xp)[0];
    (*gp)[1] = 2. * (*xp)[1];
  }

};

void
test (const double x,
      const double y)
{
  typedef double RealT;

  QuadraticObjective<RealT> quad_objective;

  Teuchos::RCP<std::ostream> outStream = Teuchos::rcp(&std::cout, false);
  Teuchos::RCP<VectorType> x_rcp       = Teuchos::rcp (new VectorType);

  x_rcp->reinit (2);

  (*x_rcp)[0] = x;
  (*x_rcp)[1] = y;

  Rol::VectorAdaptor<VectorType> x_rol(x_rcp);

  Teuchos::ParameterList parlist;
  // Set parameters.
  parlist.sublist("Secant").set("Use as Preconditioner", false);
  // Define algorithm.
  ROL::Algorithm<RealT> algo("Line Search", parlist);

  // Run Algorithm
  algo.run(x_rol, quad_objective, true, *outStream);

  Teuchos::RCP<const VectorType> xg = x_rol.getVector();
  std::cout << "The solution to minimization problem is: ";
  std::cout << (*xg)[0] << " " << (*xg)[1] << std::endl;
}

int
main (int argc, char **argv)
{
  try
    {
      test(  10,  -2);
      test(-0.1, 0.1);
      test( 9.1,-6.1);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }

  return 0;
}

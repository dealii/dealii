// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/optimization/rol/vector_adaptor.h>

#include <cmath>
#include <iostream>
#include <sstream>

#include "../tests.h"

#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Objective.hpp"
#include "ROL_StatusTest.hpp"

// Use ROL to minimize the objective function, f(x,y) = x^2 + y^2.

using namespace dealii;

using VectorType = typename dealii::Vector<double>;

template <class Real = double, typename Xprim = Rol::VectorAdaptor<VectorType>>
class QuadraticObjective : public ROL::Objective<Real>
{
private:
  ROL::Ptr<const VectorType>
  get_rolptr_to_VectorType(const ROL::Vector<Real> &x)
  {
    return (dynamic_cast<const Xprim &>(x)).getVector();
  }

  ROL::Ptr<dealii::Vector<Real>>
  get_rolptr_to_VectorType(ROL::Vector<Real> &x)
  {
    return (dynamic_cast<Xprim &>(x)).getVector();
  }

public:
  Real
  value(const ROL::Vector<Real> &x, Real & /*tol*/) override
  {
    Assert(x.dimension() == 2, ExcInternalError());

    return x.dot(x);
  }

  void
  gradient(ROL::Vector<Real>       &g,
           const ROL::Vector<Real> &x,
           Real & /*tol*/) override
  {
    ROL::Ptr<const VectorType> xp = this->get_rolptr_to_VectorType(x);
    ROL::Ptr<VectorType>       gp = this->get_rolptr_to_VectorType(g);

    (*gp)[0] = 2. * (*xp)[0];
    (*gp)[1] = 2. * (*xp)[1];
  }
};

void
test(const double x, const double y)
{
  using RealT = double;

  QuadraticObjective<RealT> quad_objective;

  ROL::Ptr<std::ostream> outStream =
    ROL::makePtrFromRef<std::ostream>(std::cout);
  ROL::Ptr<VectorType> x_rcp = ROL::makePtr<VectorType>();

  x_rcp->reinit(2);

  (*x_rcp)[0] = x;
  (*x_rcp)[1] = y;

  Rol::VectorAdaptor<VectorType> x_rol(x_rcp);

  ROL::ParameterList parlist;

#if DEAL_II_TRILINOS_VERSION_GTE(12, 18, 0)
  // Define algorithm in three intuitive and easy steps.
  //
  // For the future developer: if this ever fails again, copy the relevant
  // changes out of packages/rol/example/rosenbrock/example_01.cpp found in
  // the Trilinos git repository - the package documentation might be
  // outdated.
  ROL::Ptr<ROL::Step<RealT>> step =
    ROL::makePtr<ROL::LineSearchStep<RealT>>(parlist);
  ROL::Ptr<ROL::StatusTest<RealT>> status =
    ROL::makePtr<ROL::StatusTest<RealT>>(parlist);
  ROL::Algorithm<RealT> algo(step, status, false);
#else
  // Define algorithm.
  ROL::Algorithm<RealT> algo("Line Search", parlist);
#endif

  // Set parameters.
  parlist.sublist("Secant").set("Use as Preconditioner", false);

  // Run Algorithm.
  algo.run(x_rol, quad_objective, true, *outStream);

  ROL::Ptr<const VectorType> xg = x_rol.getVector();
  std::cout << "The solution to minimization problem is: ";
  std::cout << (*xg)[0] << ' ' << (*xg)[1] << std::endl;
}

int
main()
{
  try
    {
      test(10, -2);
      test(-0.1, 0.1);
      test(9.1, -6.1);
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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

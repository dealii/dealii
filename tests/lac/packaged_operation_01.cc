// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test PackagedOperation for
//   dealii::Vector<double>
//   dealii::SparseMatrix<double>

#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



void
test_applies(std::string                              description,
             const PackagedOperation<Vector<double>> &expr)
{
  // test apply
  Vector<double> tmp = expr;
  deallog << description << ": " << tmp << std::endl;

  // test apply_add
  for (auto &i : tmp)
    i = 100.;
  expr.apply_add(tmp);
  deallog << "100. * 1_n + " << description << ": " << tmp << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  static const int dim = 2;

  // Tests:

  Vector<double> u(25);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i + 1);
    }

  deallog << "u: " << u << std::endl;

  // creation via operator+, operator-, operator*

  test_applies("u + u", u + u);
  test_applies("u - u", u - u);
  test_applies("3. * u", 3. * u);
  test_applies("u * 3.", u * 3.);

  // creation via mixed operator+, operator-

  auto expr = 2. * u;

  test_applies("2. * u + u", expr + u);
  test_applies("2. * u - u", expr - u);

  test_applies("u + 2. * u", u + expr);
  test_applies("u - 2. * u", u - expr);

  // operator+, operator-, operator*

  PackagedOperation<Vector<double>> expr2 = u;

  test_applies("2. * u + u", expr + expr2);
  test_applies("2. * u - u", expr - expr2);
  test_applies("3. * u", 3. * expr2);
  test_applies("u * 3.", expr2 * 3.);
}

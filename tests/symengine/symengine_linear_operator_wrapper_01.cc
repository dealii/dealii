// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Check that the SymEngine wrapper for linear operators based on deal.II
// sparse linear algebra classes works as expected

#include <deal.II/differentiation/sd.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"

using namespace dealii;
namespace SD = dealii::Differentiation::SD;

template <typename number_t>
void
test_LO_Vector()
{
  typedef SD::Expression SD_number_t;

  // --- Linear operators ---

  const unsigned int rc = 10;
  SparsityPattern    sparsity_pattern(rc, rc, /*n_entries_per_row =*/1);
  for (unsigned int i = 0; i < rc; ++i)
    {
      sparsity_pattern.add(i, i);
    }
  sparsity_pattern.compress();

  // Setup a linear system
  SparseMatrix<number_t> la_A(sparsity_pattern);
  SparseMatrix<number_t> la_B(sparsity_pattern);
  Vector<number_t>       la_u(rc);
  Vector<number_t>       la_v(rc);
  for (unsigned int i = 0; i < rc; ++i)
    {
      la_A.set(i, i, 2.0);
      la_B.set(i, i, 5.0);
      la_u(i) = 2.0 * i;
      la_v(i) = 3.0 * i;
    }

  // Create linear operators
  LinearOperator<Vector<number_t>> lo_A(la_A);
  LinearOperator<Vector<number_t>> lo_B(la_B);

  // --- Values ---
  deallog.push("Initialize");
  SD_number_t A(lo_A, "A");
  SD_number_t B(lo_B, "B");
  SD_number_t u(la_u, "u");
  SD_number_t v(la_v, "v");
  deallog << "A: " << A << std::endl;
  deallog << "B: " << B << std::endl;
  deallog << "u: " << u << std::endl;
  deallog << "v: " << v << std::endl;
  deallog.pop();

  deallog.push("Relational operators (LO)");
  deallog << "A == A: " << (A == A) << std::endl;
  deallog << "A == B: " << (A == B) << std::endl;
  deallog << "A != A: " << (A != A) << std::endl;
  deallog << "A != B: " << (A != B) << std::endl;
  //   deallog << "A < B: " << (A < B) << std::endl; // Don't make sense -
  //   memory addresses could change deallog << "A > B: " << (A > B) <<
  //   std::endl;
  deallog.pop();

  deallog.push("Relational operators (Vec)");
  deallog << "u == u: " << (u == u) << std::endl;
  deallog << "u == v: " << (u == v) << std::endl;
  deallog << "u != u: " << (u != u) << std::endl;
  deallog << "u != v: " << (u != v) << std::endl;
  //   deallog << "u < v: " << (u < u) << std::endl; // Don't make sense -
  //   memory addresses could change deallog << "u > v: " << (u > v) <<
  //   std::endl;
  deallog.pop();

  //   deallog.push("Set new values");
  //   a = ADNumberType(5);
  //   b = ADNumberType(9);
  //   deallog << "a: " << a << std::endl;
  //   deallog << "b: " << b << std::endl;
  //   deallog.pop();

  //   deallog.push("Relational operators 2");
  //   deallog << "a == b: " << (a == b) << std::endl;
  //   deallog << "a != b: " << (a != b) << std::endl;
  //   deallog << "a < b: " << (a < b) << std::endl;
  //   deallog << "a > b: " << (a > b) << std::endl;
  //   deallog.pop();

  deallog.push("Math operators (LO)");
  deallog << "A+B: " << (A + B) << std::endl;
  deallog << "A-B: " << (A - B) << std::endl;
  deallog << "A*B: " << (A * B) << std::endl;

  SD_number_t C(A);
  C += B;
  deallog << "C+=B: " << C << std::endl;
  C -= B;
  deallog << "C-=B: " << C << std::endl;
  C *= B;
  deallog << "C*=B: " << C << std::endl;
  deallog.pop();

  deallog.push("Math operators (Vec)");
  deallog << "u+v: " << (u + v) << std::endl;
  deallog << "u-v: " << (u - v) << std::endl;
  deallog << "u*v: " << (u * v) << std::endl;

  SD_number_t w(u);
  w += v;
  deallog << "w+=v: " << w << std::endl;
  w -= v;
  deallog << "w-=v: " << w << std::endl;
  w *= v;
  deallog << "w*=v: " << w << std::endl;
  deallog.pop();

  deallog.push("Math operators (LO + Vec)");
  deallog << "A*u: " << (A * u) << std::endl;
  //   deallog << "v*B: " << (v*B) << std::endl; // No sensible way to implement
  //   this
  deallog << "B*A*(u-C*v): " << B * A * (u - C * v) << std::endl;
  //   deallog << "A-B: " << (A-B) << std::endl;
  //   deallog << "A*B: " << (A*B) << std::endl;
  //
  //   SD_number_t C (A);
  //   C+=B;
  //   deallog << "C+=B: " << C << std::endl;
  //   C-=B;
  //   deallog << "C-=B: " << C << std::endl;
  //   C*=B;
  //   deallog << "C*=B: " << C << std::endl;
  deallog.pop();

  // --- Normal numbers ---
  deallog.push("Math operators with normal numbers (LO)");
  const number_t s(2);
  C *= s;
  deallog << "C*=s: " << C << std::endl;
  C /= s;
  deallog << "C/=s: " << C << std::endl;
  C = std::pow(A, 3);
  deallog << "A^3: " << C << std::endl;
  deallog.pop();

  deallog.push("Math operators with normal numbers (Vec)");
  const number_t t(3);
  w = u;
  w *= t;
  deallog << "w*=s: " << w << std::endl;
  w /= t;
  deallog << "w/=s: " << w << std::endl;
  deallog.pop();

  deallog.push("Math operators with symbols (LO)");
  const SD_number_t x("x");
  deallog << "x: " << x << std::endl;
  C = A * B * x;
  deallog << "C=A*B*x: " << C << std::endl;
  //   C+=x;
  //   deallog << "C+=x: " << C << std::endl;
  C *= x;
  deallog << "C*=x: " << C << std::endl;
  //   C-=x;
  //   deallog << "C-=x: " << C << std::endl;
  C /= x;
  deallog << "C/=x: " << C << std::endl;
  deallog.pop();
}

int
main()
{
  initlog();

  deallog.push("Float");
  test_LO_Vector<float>();
  deallog.pop();

  deallog.push("Double");
  test_LO_Vector<double>();
  deallog.pop();

  deallog << "OK" << std::endl;
}

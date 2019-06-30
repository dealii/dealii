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


// Check that we can perform evaluation of a SymEngine wrapper result using
// linear operators based on deal.II sparse linear algebra classes

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

  deallog.push("Evaluation (Vec)");
  const SD_number_t w = (u + v);
  deallog << "w: " << w << std::endl;
  const Vector<number_t> w_eval = static_cast<Vector<number_t>>(w);
  deallog << "w_eval: " << w_eval << std::endl;
  deallog << "w=u+v: " << Vector<number_t>(la_u + la_v) << std::endl;
  deallog.pop();

  deallog.push("Evaluation (LO*vec)");
  const SD_number_t C = (A + B);
  deallog << "C: " << C << std::endl;
  const Vector<number_t> x_eval = static_cast<Vector<number_t>>(C * w);
  deallog << "x_eval: " << x_eval << std::endl;
  deallog << "x=(A*B)*(u+v): "
          << Vector<number_t>((lo_A + lo_B) * (la_u + la_v)) << std::endl;
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

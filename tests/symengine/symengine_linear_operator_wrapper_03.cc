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


// Check that we can perform substitution and subsequent evaluation of a
// SymEngine wrapper result using linear operators based on deal.II sparse
// linear algebra classes

#include <deal.II/differentiation/sd.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_selector.h>
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
  LinearOperator<Vector<number_t>> lo_A_T = transpose_operator(lo_A);

  PreconditionSelector<SparseMatrix<number_t>, Vector<number_t>>
    preconditioner_A_inv("jacobi");
  preconditioner_A_inv.use_matrix(la_A);
  ReductionControl solver_control_A_inv(la_A.m(), 1.0e-30, 1e-6);
  SolverSelector<Vector<number_t>> solver_A_inv;
  solver_A_inv.select("cg");
  solver_A_inv.set_control(solver_control_A_inv);
  LinearOperator<Vector<number_t>> lo_A_inv =
    inverse_operator(lo_A, solver_A_inv, preconditioner_A_inv);

  deallog.push("Symbolic operations (LO + Vec)");
  const SD_number_t symb_lo_A("lo_A");
  const SD_number_t symb_lo_B("lo_B");
  const SD_number_t symb_lo_A_T("lo_A_T");
  const SD_number_t symb_lo_A_inv("lo_A_inv");
  const SD_number_t symb_la_u("la_u");
  const SD_number_t symb_la_v("la_v");
  deallog << "symb_lo_A: " << symb_lo_A << std::endl;
  deallog << "symb_lo_B: " << symb_lo_B << std::endl;
  deallog << "symb_lo_A_T: " << symb_lo_A_T << std::endl;
  deallog << "symb_lo_A_inv: " << symb_lo_A_inv << std::endl;
  deallog << "symb_la_u: " << symb_la_u << std::endl;
  deallog << "symb_la_v: " << symb_la_v << std::endl;

  const SD_number_t symb_la_w = symb_lo_A * symb_la_u;
  deallog << "symb_la_w: " << symb_la_w << std::endl;

  const SD_number_t symb_la_x =
    (symb_lo_A + symb_lo_B) * (symb_la_u + symb_la_v);
  deallog << "symb_la_x: " << symb_la_x << std::endl;

  // Note: The internal representation of this result violates the
  // non-commutative property of these operations.
  const SD_number_t symb_la_y = symb_lo_A_inv * (symb_la_u + symb_la_v);
  deallog << "symb_la_y: " << symb_la_y << std::endl;

  // Who knows what this operation could possibly represent, but we'll try it
  // anyway...
  const SD_number_t symb_la_z =
    (2.0 * symb_lo_A_inv + symb_lo_A_T * symb_lo_A - 0.5 * symb_lo_B) *
    (symb_la_u / 2.0 + 2.0 * symb_la_v);
  deallog << "symb_la_z: " << symb_la_z << std::endl;

  // This is an ambiguous case, testing the commutative property of the
  // operations i.e. Will it return () or ()?...
  //  const SD_number_t symb_la_t = symb_la_v*(symb_lo_A*symb_la_u);
  //  deallog << "symb_la_t: " << symb_la_t << std::endl;
  deallog.pop();

  deallog.push("Substitution");
  SD_number_t::substitution_map_t sub_vals;
  sub_vals[symb_lo_A.get_RCP()]     = SD_number_t(lo_A, "A").get_RCP();
  sub_vals[symb_lo_B.get_RCP()]     = SD_number_t(lo_B, "B").get_RCP();
  sub_vals[symb_lo_A_T.get_RCP()]   = SD_number_t(lo_A_T, "A_T").get_RCP();
  sub_vals[symb_lo_A_inv.get_RCP()] = SD_number_t(lo_A_inv, "A_inv").get_RCP();
  sub_vals[symb_la_u.get_RCP()]     = SD_number_t(la_u, "u").get_RCP();
  sub_vals[symb_la_v.get_RCP()]     = SD_number_t(la_v, "v").get_RCP();

  deallog << "lo_A: " << symb_lo_A.substitute(sub_vals) << std::endl;
  deallog << "lo_B: " << symb_lo_B.substitute(sub_vals) << std::endl;
  deallog << "la_u: " << symb_la_u.substitute(sub_vals) << std::endl;
  deallog << "la_v: " << symb_la_v.substitute(sub_vals) << std::endl;
  deallog << "la_w: " << symb_la_w.substitute(sub_vals) << std::endl;
  deallog << "la_x: " << symb_la_x.substitute(sub_vals) << std::endl;
  deallog << "la_y: " << symb_la_y.substitute(sub_vals) << std::endl;
  deallog << "la_z: " << symb_la_z.substitute(sub_vals) << std::endl;
  //  deallog << "la_t: " << symb_la_t.substitute(sub_vals) << std::endl;
  deallog.pop();

  deallog.push("Evaluation");
  const Vector<number_t> la_u_result =
    static_cast<Vector<number_t>>(symb_la_u.substitute(sub_vals));
  const Vector<number_t> la_v_result =
    static_cast<Vector<number_t>>(symb_la_v.substitute(sub_vals));
  const Vector<number_t> la_w_result =
    static_cast<Vector<number_t>>(symb_la_w.substitute(sub_vals));
  const Vector<number_t> la_x_result =
    static_cast<Vector<number_t>>(symb_la_x.substitute(sub_vals));
  const Vector<number_t> la_y_result =
    static_cast<Vector<number_t>>(symb_la_y.substitute(sub_vals));
  const Vector<number_t> la_z_result =
    static_cast<Vector<number_t>>(symb_la_z.substitute(sub_vals));
  //  const Vector<number_t> la_t_result = static_cast< Vector<number_t>
  //  >(symb_la_t.substitute(sub_vals));

  deallog << "la_u: " << la_u_result << std::endl;
  deallog << "la_v: " << la_v_result << std::endl;
  deallog << "la_w: " << la_w_result << std::endl;
  deallog << "la_x: " << la_x_result << std::endl;
  deallog << "la_y: " << la_y_result << std::endl;
  deallog << "la_z: " << la_z_result << std::endl;
  //  deallog << "la_t: " << la_t_result << std::endl;
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

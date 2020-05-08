// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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


// same as gmres_reorthogonalize_04 but with forced orthogonalization

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



// reimplements things from vector.templates.h in another way which is more
// prone to roundoff, the linker will select this version instead of the one
// in the deal.II library
namespace dealii
{
  template <typename Number>
  template <typename Number2>
  Number Vector<Number>::operator*(const Vector<Number2> &v) const
  {
    Number sum = 0;
    for (unsigned int i = 0; i < size(); ++i)
      sum += values[i] * v.values[i];
    return sum;
  }
  template <typename Number>
  typename Vector<Number>::real_type
  Vector<Number>::l2_norm() const
  {
    real_type sum = 0;
    for (unsigned int i = 0; i < size(); ++i)
      sum += values[i] * values[i];
    return std::sqrt(sum);
  }
} // namespace dealii



template <typename number>
void
test()
{
  const unsigned int n = 200;
  Vector<number>     rhs(n), sol(n);
  rhs = 1.;

  // only add diagonal entries
  SparsityPattern sp(n, n);
  sp.compress();
  SparseMatrix<number> matrix(sp);

  for (unsigned int i = 0; i < n; ++i)
    matrix.diag_element(i) = (i + 1);

  SolverControl control(1000, 1e2 * std::numeric_limits<number>::epsilon());
  typename SolverGMRES<Vector<number>>::AdditionalData data;
  data.max_n_tmp_vectors          = 202;
  data.force_re_orthogonalization = true;

  SolverGMRES<Vector<number>> solver(control, data);
  auto print_re_orthogonalization = [](int accumulated_iterations) {
    deallog.get_file_stream() << "Re-orthogonalization enabled at step "
                              << accumulated_iterations << std::endl;
  };
  solver.connect_re_orthogonalization_slot(print_re_orthogonalization);
  solver.solve(matrix, sol, rhs, PreconditionIdentity());
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog.push("double");
  test<double>();
  deallog.pop();
  deallog.push("float");
  test<float>();
  deallog.pop();
}

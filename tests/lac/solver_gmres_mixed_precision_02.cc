// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Same as solver_gmres_mixed_precision_01, but for SolverFGMRES and
// SolverMPGMRES, the latter with two preconditioners and the truncated
// strategy, which applies the preconditioners also to older basis vectors
// that need to be converted from single to double precision. (The full
// strategy is not tested because the diagonal matrix and preconditioners of
// this test commute, which leads to linearly dependent search directions.)

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



template <typename VectorType>
struct MyDiagonalMatrix
{
  MyDiagonalMatrix(const VectorType &diagonal)
    : diagonal(diagonal)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = src;
    dst.scale(diagonal);
  }

  const VectorType &diagonal;
};



template <typename Number, typename MemorySpace>
void
fill_vector(LinearAlgebra::distributed::Vector<Number, MemorySpace> &vec,
            const std::function<double(unsigned int)>               &function)
{
  LinearAlgebra::ReadWriteVector<Number> rw(vec.locally_owned_elements());
  for (unsigned int i = 0; i < vec.size(); ++i)
    rw(i) = function(i);
  vec.import_elements(rw, VectorOperation::insert);
}



template <typename Number>
void
fill_vector(LinearAlgebra::distributed::BlockVector<Number> &vec,
            const std::function<double(unsigned int)>       &function)
{
  unsigned int offset = 0;
  for (unsigned int b = 0; b < vec.n_blocks(); ++b)
    {
      fill_vector(vec.block(b),
                  [&](const unsigned int i) { return function(offset + i); });
      offset += vec.block(b).size();
    }
}



template <typename Number>
void
fill_vector(Vector<Number>                            &vec,
            const std::function<double(unsigned int)> &function)
{
  for (unsigned int i = 0; i < vec.size(); ++i)
    vec(i) = function(i);
}



template <typename VectorType>
void
create_vector(VectorType &vec, const unsigned int size)
{
  if constexpr (IsBlockVector<VectorType>::value)
    {
      vec.reinit(2);
      vec.block(0).reinit(size / 2);
      vec.block(1).reinit(size - size / 2);
      vec.collect_sizes();
    }
  else
    vec.reinit(size);
}



template <typename VectorType, typename BasisVectorType>
void
test(const LinearAlgebra::OrthogonalizationStrategy strategy,
     const unsigned int                             variant)
{
  // choose a size that is not a multiple of the tile sizes in the kernels
  const unsigned int size = 20013;
  VectorType         diagonal, inverse_diagonal, rhs, sol, exact;
  create_vector(diagonal, size);
  create_vector(inverse_diagonal, size);
  create_vector(rhs, size);
  create_vector(sol, size);
  create_vector(exact, size);
  fill_vector(diagonal,
              [](const unsigned int i) { return 1. + 50. * (i % 211) / 211.; });
  // a rough approximation of the inverse diagonal as second preconditioner
  fill_vector(inverse_diagonal, [](const unsigned int i) {
    return 1. / (1. + 40. * (i % 211) / 211.);
  });
  fill_vector(rhs, [](const unsigned int i) { return 1. + 0.01 * (i % 7); });
  fill_vector(exact, [](const unsigned int i) {
    return (1. + 0.01 * (i % 7)) / (1. + 50. * (i % 211) / 211.);
  });

  MyDiagonalMatrix<VectorType> matrix(diagonal);
  MyDiagonalMatrix<VectorType> preconditioner(inverse_diagonal);

  // tolerance well below the precision of the basis vectors
  SolverControl control(1000, 1e-11 * rhs.l2_norm());

  // Variants 0-1 run FGMRES and variants 2-3 MPGMRES, each with a small
  // basis with many restarts and with a large basis without restarts
  if (variant < 2)
    {
      typename SolverFGMRES<VectorType, BasisVectorType>::AdditionalData data;
      data.max_basis_size             = variant == 0 ? 19 : 150;
      data.orthogonalization_strategy = strategy;
      SolverFGMRES<VectorType, BasisVectorType> solver(control, data);
      solver.solve(matrix, sol, rhs, PreconditionIdentity());
    }
  else
    {
      typename SolverMPGMRES<VectorType, BasisVectorType>::AdditionalData data;
      data.max_basis_size                 = variant == 2 ? 9 : 150;
      data.orthogonalization_strategy     = strategy;
      data.use_truncated_mpgmres_strategy = true;
      SolverMPGMRES<VectorType, BasisVectorType> solver(control, data);
      solver.solve(matrix, sol, rhs, PreconditionIdentity(), preconditioner);
    }

  sol -= exact;
  const double relative_error = sol.l2_norm() / exact.l2_norm();
  deallog << (variant == 0 ? "FGMRES" :
              variant == 1 ? "FGMRES, no restart" :
              variant == 2 ? "MPGMRES" :
                             "MPGMRES, no restart")
          << " converged in " << control.last_step() << " steps, error "
          << (relative_error < 1e-9 ? "ok" : std::to_string(relative_error))
          << std::endl;
}



template <typename VectorType, typename BasisVectorType>
void
test_all()
{
  deallog.push("CGS");
  for (unsigned int variant = 0; variant < 4; ++variant)
    test<VectorType, BasisVectorType>(
      LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt,
      variant);
  deallog.pop();
  deallog.push("DCGS");
  for (unsigned int variant = 0; variant < 4; ++variant)
    test<VectorType, BasisVectorType>(
      LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt,
      variant);
  deallog.pop();
  deallog.push("MGS");
  for (unsigned int variant = 0; variant < 4; ++variant)
    test<VectorType, BasisVectorType>(
      LinearAlgebra::OrthogonalizationStrategy::modified_gram_schmidt, variant);
  deallog.pop();
}



int
main()
{
  initlog();

  {
    using VectorType      = LinearAlgebra::distributed::Vector<double>;
    using BasisVectorType = LinearAlgebra::distributed::Vector<float>;
    deallog.push("double basis");
    test_all<VectorType, VectorType>();
    deallog.pop();
    deallog.push("float basis");
    test_all<VectorType, BasisVectorType>();
    deallog.pop();
  }

  {
    using VectorType =
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default>;
    using BasisVectorType =
      LinearAlgebra::distributed::Vector<float, MemorySpace::Default>;
    deallog.push("float basis default");
    test_all<VectorType, BasisVectorType>();
    deallog.pop();
  }

  {
    using VectorType      = LinearAlgebra::distributed::BlockVector<double>;
    using BasisVectorType = LinearAlgebra::distributed::BlockVector<float>;
    deallog.push("float basis block");
    test_all<VectorType, BasisVectorType>();
    deallog.pop();
  }

  {
    deallog.push("float basis serial");
    test_all<Vector<double>, Vector<float>>();
    deallog.pop();
  }
}

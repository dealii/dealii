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

// Check the optimized orthogonalization kernels in SolverGMRES and
// SolverFGMRES for LinearAlgebra::distributed::Vector with
// MemorySpace::Default. The iteration counts must coincide with the ones
// obtained with MemorySpace::Host, which uses a different implementation.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>

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



template <typename Number, typename MemorySpace>
void
test(const LinearAlgebra::OrthogonalizationStrategy strategy,
     const unsigned int                             variant)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;

  // choose a size that is not a multiple of the tile sizes in the kernels
  const unsigned int size = 20013;
  VectorType         diagonal(size), rhs(size), sol(size), exact(size);
  fill_vector(diagonal,
              [](const unsigned int i) { return 1. + 50. * (i % 211) / 211.; });
  fill_vector(rhs, [](const unsigned int i) { return 1. + 0.01 * (i % 7); });
  fill_vector(exact, [](const unsigned int i) {
    return (1. + 0.01 * (i % 7)) / (1. + 50. * (i % 211) / 211.);
  });

  MyDiagonalMatrix<VectorType> matrix(diagonal);

  SolverControl control(1000,
                        (std::is_same_v<Number, float> ? 1e-5 : 1e-10) *
                          rhs.l2_norm());

  // a small basis size gives several restarts and a basis whose size is not
  // a multiple of the block size in the kernels
  const unsigned int basis_size = 19;

  if (variant < 2)
    {
      typename SolverGMRES<VectorType>::AdditionalData data;
      data.max_basis_size             = basis_size;
      data.orthogonalization_strategy = strategy;
      data.right_preconditioning      = (variant == 1);
      SolverGMRES<VectorType> solver(control, data);
      solver.solve(matrix, sol, rhs, PreconditionIdentity());
    }
  else
    {
      typename SolverFGMRES<VectorType>::AdditionalData data;
      data.max_basis_size             = basis_size;
      data.orthogonalization_strategy = strategy;
      SolverFGMRES<VectorType> solver(control, data);
      solver.solve(matrix, sol, rhs, PreconditionIdentity());
    }

  sol -= exact;
  const double relative_error = sol.l2_norm() / exact.l2_norm();
  deallog << (variant == 0 ? "GMRES left" :
              variant == 1 ? "GMRES right" :
                             "FGMRES")
          << " converged in " << control.last_step() << " steps, error "
          << (relative_error < (std::is_same_v<Number, float> ? 1e-4 : 1e-9) ?
                "ok" :
                std::to_string(relative_error))
          << std::endl;
}



template <typename Number, typename MemorySpace>
void
test_all()
{
  for (const auto strategy :
       {LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt,
        LinearAlgebra::OrthogonalizationStrategy::
          delayed_classical_gram_schmidt,
        LinearAlgebra::OrthogonalizationStrategy::modified_gram_schmidt})
    {
      deallog.push(
        strategy ==
            LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt ?
          "CGS" :
        strategy == LinearAlgebra::OrthogonalizationStrategy::
                      delayed_classical_gram_schmidt ?
          "DCGS" :
          "MGS");
      for (unsigned int variant = 0; variant < 3; ++variant)
        test<Number, MemorySpace>(strategy, variant);
      deallog.pop();
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  deallog.push("double host");
  test_all<double, MemorySpace::Host>();
  deallog.pop();

  deallog.push("double default");
  test_all<double, MemorySpace::Default>();
  deallog.pop();

  deallog.push("float host");
  test_all<float, MemorySpace::Host>();
  deallog.pop();

  deallog.push("float default");
  test_all<float, MemorySpace::Default>();
  deallog.pop();
}

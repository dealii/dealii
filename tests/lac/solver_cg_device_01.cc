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

// Check the fused kernels in SolverCG for LinearAlgebra::distributed::Vector
// with MemorySpace::Default and the preconditioners PreconditionIdentity and
// DiagonalMatrix. The iteration counts must coincide with the ones obtained
// with MemorySpace::Host, where the same algorithm is run with the vector
// updates interleaved into the matrix-vector product, and the solution must
// satisfy the linear system up to the requested tolerance.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <Kokkos_Core.hpp>

#include "../tests.h"



// A symmetric positive definite tridiagonal matrix with variable diagonal
// entries. The coupling between the entries is only within the locally owned
// range of each MPI process, such that no ghost exchange is necessary.
template <typename VectorType>
struct MyTridiagonalMatrix
{
  using Number = typename VectorType::value_type;

  MyTridiagonalMatrix(const VectorType &diagonal)
    : diagonal(diagonal)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    const std::size_t size    = src.locally_owned_size();
    const Number     *diag    = diagonal.get_values();
    const Number     *src_ptr = src.get_values();
    Number           *dst_ptr = dst.get_values();
    if constexpr (std::is_same_v<typename VectorType::memory_space,
                                 MemorySpace::Default>)
      {
        using ExecutionSpace =
          MemorySpace::Default::kokkos_space::execution_space;
        Kokkos::parallel_for(
          "test_tridiagonal_vmult",
          Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<std::size_t>>(
            ExecutionSpace(), 0, size),
          KOKKOS_LAMBDA(const std::size_t i) {
            Number result = diag[i] * src_ptr[i];
            if (i > 0)
              result -= src_ptr[i - 1];
            if (i + 1 < size)
              result -= src_ptr[i + 1];
            dst_ptr[i] = result;
          });
      }
    else
      for (std::size_t i = 0; i < size; ++i)
        {
          Number result = diag[i] * src_ptr[i];
          if (i > 0)
            result -= src_ptr[i - 1];
          if (i + 1 < size)
            result -= src_ptr[i + 1];
          dst_ptr[i] = result;
        }
  }

  // Variant with the functions run before and after the matrix-vector
  // product, enabling the specialized algorithm on the host
  void
  vmult(
    VectorType                                                        &dst,
    const VectorType                                                  &src,
    const std::function<void(const unsigned int, const unsigned int)> &before,
    const std::function<void(const unsigned int, const unsigned int)> &after)
    const
  {
    before(0, dst.locally_owned_size());
    vmult(dst, src);
    after(0, dst.locally_owned_size());
  }

  const VectorType &diagonal;
};



template <typename Number, typename MemorySpace>
void
fill_vector(LinearAlgebra::distributed::Vector<Number, MemorySpace> &vec,
            const std::function<double(unsigned int)>               &function)
{
  LinearAlgebra::ReadWriteVector<Number> rw(vec.locally_owned_elements());
  for (const auto i : vec.locally_owned_elements())
    rw(i) = function(i);
  vec.import_elements(rw, VectorOperation::insert);
}



template <typename Number, typename MemorySpace>
void
test(const unsigned int size, const bool use_diagonal_preconditioner)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;

  const MPI_Comm     comm    = MPI_COMM_WORLD;
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(comm);
  const unsigned int my_rank = Utilities::MPI::this_mpi_process(comm);
  IndexSet           owned(size);
  owned.add_range(static_cast<std::size_t>(size) * my_rank / n_procs,
                  static_cast<std::size_t>(size) * (my_rank + 1) / n_procs);

  VectorType diagonal(owned, comm), rhs(owned, comm), sol(owned, comm);
  fill_vector(diagonal, [](const unsigned int i) {
    return 2.02 + 50. * (i % 211) / 211.;
  });
  fill_vector(rhs, [](const unsigned int i) { return 1. + 0.01 * (i % 7); });

  MyTridiagonalMatrix<VectorType> matrix(diagonal);

  DiagonalMatrix<VectorType> preconditioner;
  preconditioner.get_vector().reinit(owned, comm);
  fill_vector(preconditioner.get_vector(), [](const unsigned int i) {
    return 1. / (2.02 + 50. * (i % 211) / 211.);
  });

  const double         tolerance = std::is_same_v<Number, float> ? 1e-5 : 1e-10;
  SolverControl        control(1000, tolerance * rhs.l2_norm());
  SolverCG<VectorType> solver(control);
  if (use_diagonal_preconditioner)
    solver.solve(matrix, sol, rhs, preconditioner);
  else
    solver.solve(matrix, sol, rhs, PreconditionIdentity());

  // check the true residual
  VectorType residual(owned, comm);
  matrix.vmult(residual, sol);
  residual -= rhs;
  const double relative_residual = residual.l2_norm() / rhs.l2_norm();
  deallog << "size " << size << ", "
          << (use_diagonal_preconditioner ? "diagonal" : "identity")
          << ": converged in " << control.last_step() << " steps, residual "
          << (relative_residual < 20. * tolerance ?
                "ok" :
                std::to_string(relative_residual))
          << std::endl;
}



template <typename Number, typename MemorySpace>
void
test_all()
{
  // choose sizes so that the solver stops both after an even and an odd
  // number of iterations
  for (const unsigned int size : {20013U, 1001U, 37U})
    for (const bool use_diagonal : {false, true})
      test<Number, MemorySpace>(size, use_diagonal);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

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

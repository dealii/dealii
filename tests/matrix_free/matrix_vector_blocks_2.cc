// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test the correctness of matrix free matrix-vector product with block vectors
// consisting of many blocks with respect to the MPI data exchange
//
// @note Unlike `matrix_vector_blocks.cc` this tests the explicit case.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <iostream>
#include <limits>

#include "../tests.h"

// Simple operator that adds one (i.e., \phi^n = 1 + \phi^n-1).
template <int dim, int fe_degree, typename Number>
class MatrixFreeBlock
{
public:
  MatrixFreeBlock(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
  {
    data.initialize_dof_vector(invm);

    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> fe_eval(data);
    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        fe_eval.reinit(cell);

        for (const unsigned int q : fe_eval.quadrature_point_indices())
          fe_eval.submit_value(make_vectorized_array(Number(1.0)), q);

        fe_eval.integrate_scatter(EvaluationFlags::values, invm);
      }

    invm.compress(VectorOperation::add);
    for (unsigned int k = 0; k < invm.locally_owned_size(); ++k)
      invm.local_element(k) =
        invm.local_element(k) > 10.0 * std::numeric_limits<Number>::epsilon() ?
          1.0 / invm.local_element(k) :
          1.0;
  }

  void
  apply(LinearAlgebra::distributed::BlockVector<Number>       &dst,
        const LinearAlgebra::distributed::BlockVector<Number> &src) const
  {
    AssertDimension(src.n_blocks(), dst.n_blocks());

    data.cell_loop(&MatrixFreeBlock::local_apply, this, dst, src, true);

    for (unsigned int block = 0; block < dst.n_blocks(); ++block)
      dst.block(block).scale(invm);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number>                         &data,
              LinearAlgebra::distributed::BlockVector<Number>       &dst,
              const LinearAlgebra::distributed::BlockVector<Number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        for (unsigned int block = 0; block < src.n_blocks(); ++block)
          {
            phi.gather_evaluate(src.block(block), EvaluationFlags::values);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(Number(1) + phi.get_value(q), q);

            phi.integrate_scatter(EvaluationFlags::values, dst.block(block));
          }
      }
  }

  const MatrixFree<dim, Number>             &data;
  LinearAlgebra::distributed::Vector<Number> invm;
};



template <int dim, int fe_degree>
void
test()
{
  using number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  const IndexSet &owned_set    = dof.locally_owned_dofs();
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dof);

  AffineConstraints<double> constraints(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGaussLobatto<1>                           quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeBlock<dim, fe_degree, number> mf(mf_data);

  // make sure that the value we set here at least includes some case where we
  // need to go to the alternative case of calling the full
  // update_ghost_values()
  Assert(
    LinearAlgebra::distributed::BlockVector<number>::communication_block_size <
      80,
    ExcInternalError());
  for (unsigned int n_blocks = 5; n_blocks < 81; n_blocks *= 2)
    {
      LinearAlgebra::distributed::BlockVector<number> phi, phi_old;
      {
        std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          partitioners(n_blocks, mf_data.get_vector_partitioner());

        phi.reinit(partitioners);
        phi_old.reinit(partitioners);
      }

      // after initialization via MatrixFree we are in the write state for
      // ghosts
      AssertThrow(!phi.has_ghost_elements(), ExcInternalError());

      deallog << "Average value with " << n_blocks << " blocks:";

      // run 10 times to make a possible error more likely to show up
      for (unsigned int run = 0; run < 10; ++run)
        {
          mf.apply(phi, phi_old);
          AssertThrow(!phi.has_ghost_elements(), ExcInternalError());

          const double avg_val = phi.l1_norm() / (n_blocks * dof.n_dofs());
          deallog << ' ' << avg_val;

          phi.swap(phi_old);
        }
      deallog << std::endl;
    }
  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  mpi_initlog();
  deallog << std::setprecision(4);

  deallog.push("2d");
  test<2, 2>();
  deallog.pop();

  deallog.push("3d");
  test<3, 2>();
  deallog.pop();
}

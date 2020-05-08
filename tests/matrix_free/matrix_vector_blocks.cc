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



// test the correctness of matrix free matrix-vector product with block vectors
// consisting of many blocks with respect to the MPI data exchange

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

#include "matrix_vector_mf.h"



template <int dim, int fe_degree, typename Number>
class MatrixFreeBlock
{
public:
  MatrixFreeBlock(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
  {}

  void
  vmult(LinearAlgebra::distributed::BlockVector<Number> &      dst,
        const LinearAlgebra::distributed::BlockVector<Number> &src) const
  {
    data.cell_loop(&MatrixFreeBlock::local_apply, this, dst, src, true);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number> &                        data,
              LinearAlgebra::distributed::BlockVector<Number> &      dst,
              const LinearAlgebra::distributed::BlockVector<Number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    AssertDimension(src.n_blocks(), dst.n_blocks());
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        for (unsigned int block = 0; block < src.n_blocks(); ++block)
          {
            phi.gather_evaluate(src.block(block), true, true, false);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                phi.submit_value(Number(10) * phi.get_value(q), q);
                phi.submit_gradient(phi.get_gradient(q), q);
              }
            phi.integrate_scatter(true, true, dst.block(block));
          }
      }
  }

  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree>
void
test()
{
  typedef double number;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->is_locally_owned())
      if (cell->center().norm() < 0.2)
        cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(4 - dim);
  if (tria.begin(tria.n_levels() - 1)->is_locally_owned())
    tria.begin(tria.n_levels() - 1)->set_refine_flag();
  if (tria.last()->is_locally_owned())
    tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 11 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (cell->is_locally_owned())
          if (counter % (7 - i) == 0)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dof, relevant_set);

  AffineConstraints<double> constraints(relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.overlap_communication_computation = false;
    mf_data.reinit(dof, constraints, quad, data);
  }

  MatrixFreeTest<dim,
                 fe_degree,
                 number,
                 LinearAlgebra::distributed::Vector<number>>
                                          mf_ref(mf_data);
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
      LinearAlgebra::distributed::BlockVector<number> in, out, ref;
      in.reinit(n_blocks);
      for (unsigned int block = 0; block < n_blocks; ++block)
        mf_data.initialize_dof_vector(in.block(block));
      out.reinit(in);
      ref.reinit(in);

      // after initialization via MatrixFree we are in the write state for
      // ghosts
      AssertThrow(in.has_ghost_elements() == false, ExcInternalError());

      // fill each block with random numbers and do the block-wise
      // matrix-vector product for reference
      for (unsigned int block = 0; block < n_blocks; ++block)
        {
          for (unsigned int i = 0; i < in.block(block).local_size(); ++i)
            {
              const unsigned int glob_index = owned_set.nth_index_in_set(i);
              if (constraints.is_constrained(glob_index))
                continue;
              in.block(block).local_element(i) = random_value<double>();
            }
          mf_ref.vmult(ref.block(block), in.block(block));
        }

      // explicitly update ghosts so that we can read them:
      in.update_ghost_values();
      AssertThrow(in.has_ghost_elements(), ExcInternalError());

      deallog << "Norm of difference with " << n_blocks << " blocks:";

      // run 10 times to make a possible error more
      // likely to show up
      for (unsigned int run = 0; run < 10; ++run)
        {
          mf.vmult(out, in);

          // since we made "in" ghosted, make sure it is still ghosted
          AssertThrow(in.has_ghost_elements(), ExcInternalError());

          out -= ref;
          const double diff_norm = out.linfty_norm();
          deallog << " " << diff_norm;
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

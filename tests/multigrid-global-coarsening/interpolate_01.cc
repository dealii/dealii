// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * Test MGTransferGlobalCoarsening::interpolate_to_mg() for FE_Q and FE_DGQ.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
create_quadrant(Triangulation<dim> &tria, const unsigned int n_refinements)
{
  GridGenerator::hyper_cube(tria, -1.0, +1.0);

  if (n_refinements == 0)
    return;

  tria.refine_global(1);

  for (unsigned int i = 1; i < n_refinements; i++)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (int d = 0; d < dim; d++)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      tria.execute_coarsening_and_refinement();
    }

  AssertDimension(tria.n_global_levels() - 1, n_refinements);
}



template <int dim>
std::shared_ptr<Utilities::MPI::Partitioner>
create_partitioner(const DoFHandler<dim> &dof_handler)
{
  return std::make_shared<Utilities::MPI::Partitioner>(
    dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_active_dofs(dof_handler),
    dof_handler.get_mpi_communicator());
}



template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;

    return p[0];
  }
};



template <int dim, typename Number = double>
void
test(const unsigned int n_refinements,
     const unsigned int fe_degree_fine,
     const bool         do_dg)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_quadrant(tria, n_refinements);

  const auto trias =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(tria);

  MGLevelObject<DoFHandler<dim>> dof_handlers(min_level, max_level);
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers(min_level,
                                                               max_level);

  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &dof_handler = dof_handlers[l];

      dof_handler.reinit(*trias[l]);

      if (do_dg)
        dof_handler.distribute_dofs(FE_DGQ<dim>{fe_degree_fine});
      else
        dof_handler.distribute_dofs(FE_Q<dim>{fe_degree_fine});
    }

  for (unsigned int l = min_level; l < max_level; ++l)
    transfers[l + 1].reinit(dof_handlers[l + 1], dof_handlers[l]);

  MGTransferGlobalCoarsening<dim, VectorType> transfer(
    transfers, [&](const auto l, auto &vec) {
      vec.reinit(create_partitioner(dof_handlers[l]));
    });

  VectorType vec;
  vec.reinit(create_partitioner(dof_handlers[max_level]));

  RightHandSideFunction<dim> fu;

  VectorTools::interpolate(dof_handlers[max_level], fu, vec);

  MGLevelObject<VectorType> results(min_level, max_level);

  transfer.interpolate_to_mg(dof_handlers[max_level], results, vec);

  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      results[l].update_ghost_values();
      Vector<float> norm_per_cell(trias[l]->n_active_cells());
      VectorTools::integrate_difference(dof_handlers[l],
                                        results[l],
                                        fu,
                                        norm_per_cell,
                                        QGauss<dim>(fe_degree_fine + 2),
                                        VectorTools::L2_norm);
      const double error_L2_norm =
        VectorTools::compute_global_error(*trias[l],
                                          norm_per_cell,
                                          VectorTools::L2_norm);


      deallog << l << " " << results[l].l2_norm() << " " << error_L2_norm
              << std::endl;
    }

  deallog << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  deallog.precision(8);

  for (unsigned int n_refinements = 2; n_refinements <= 4; ++n_refinements)
    for (unsigned int degree = 1; degree <= 4; ++degree)
      {
        test<2>(n_refinements, degree, true);
        test<2>(n_refinements, degree, false);
      }
}

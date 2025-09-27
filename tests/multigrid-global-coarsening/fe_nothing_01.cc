// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
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
 * Test MGTwoLevelTransfer::interpolate() for DoFHandlers with `FE_Nothing`.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
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

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim>
class AnalyticalSolution : public Function<dim>
{
public:
  AnalyticalSolution()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    (void)component;

    return p[0];
  }

private:
};

template <int dim, typename Number = double>
void
do_test(const hp::FECollection<dim> &fe_fine,
        const hp::FECollection<dim> &fe_coarse)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  // create grid
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria);

  if (fe_fine.size() > 1)
    {
      for (const auto &cell : dof_handler_fine.active_cell_iterators())
        if (cell->center()[0] < 0.5)
          cell->set_active_fe_index(0);
        else
          cell->set_active_fe_index(1);
    }

  dof_handler_fine.distribute_dofs(fe_fine);

  DoFHandler<dim> dof_handler_coarse(tria);

  if (fe_coarse.size() > 1)
    {
      for (const auto &cell : dof_handler_coarse.active_cell_iterators())
        if (cell->center()[0] < 0.5)
          cell->set_active_fe_index(0);
        else
          cell->set_active_fe_index(1);
    }
  dof_handler_coarse.distribute_dofs(fe_coarse);

  AffineConstraints<Number> dummy;
  dummy.close();

  // setup transfer operator
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;
  transfer.reinit_polynomial_transfer(dof_handler_fine,
                                      dof_handler_coarse,
                                      dummy,
                                      dummy);

  LinearAlgebra::distributed::Vector<Number> vec_fine(
    dof_handler_fine.n_dofs());
  LinearAlgebra::distributed::Vector<Number> vec_coarse(
    dof_handler_coarse.n_dofs());

  Tensor<1, dim, Number> exp;
  exp[0] = 1.0;

  VectorTools::interpolate(dof_handler_fine,
                           Functions::Monomial<dim>(exp),
                           vec_fine);

  transfer.interpolate(vec_coarse, vec_fine);

  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler_coarse, vec_coarse, "solution_coarse");
  data_out.add_data_vector(dof_handler_fine, vec_fine, "solution_fine");

  data_out.build_patches();

  static unsigned int counter = 0;

#if 0
  data_out.write_vtu_with_pvtu_record(
    "./", "result", counter++, tria.get_mpi_communicator(), 3, 1);
#else
  deallog << std::endl;
  data_out.write_vtk(deallog.get_file_stream());
#endif
}

template <int dim, typename Number>
void
test(const FiniteElement<dim> &fe_fine, const FiniteElement<dim> &fe_coarse)
{
  deallog.push("0");
  do_test(hp::FECollection<dim>(fe_fine), hp::FECollection<dim>(fe_coarse));
  deallog.pop();
  deallog.push("1");
  do_test(hp::FECollection<dim>(fe_fine, FE_Nothing<dim>()),
          hp::FECollection<dim>(fe_coarse));
  deallog.pop();
  deallog.push("2");
  do_test(hp::FECollection<dim>(fe_fine),
          hp::FECollection<dim>(fe_coarse, FE_Nothing<dim>()));
  deallog.pop();
  deallog.push("3");
  do_test(hp::FECollection<dim>(fe_fine, FE_Nothing<dim>()),
          hp::FECollection<dim>(fe_coarse, FE_Nothing<dim>()));
  deallog.pop();
}

template <int dim, typename Number>
void
test()
{
  const unsigned int fe_degree = 2;

  deallog.push("CG<->CG");
  test<dim, Number>(FE_Q<dim>(fe_degree), FE_Q<dim>(fe_degree));
  deallog.pop();
  deallog.push("DG<->CG");
  test<dim, Number>(FE_DGQ<dim>(fe_degree), FE_Q<dim>(fe_degree));
  deallog.pop();
  deallog.push("CG<->DG");
  test<dim, Number>(FE_Q<dim>(fe_degree), FE_DGQ<dim>(fe_degree));
  deallog.pop();
  deallog.push("DG<->DG");
  test<dim, Number>(FE_DGQ<dim>(fe_degree), FE_DGQ<dim>(fe_degree));
  deallog.pop();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  test<2, double>();
}

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

// Test MappingQCache initialization with point lambda/Function

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class Solution : public Function<dim>
{
public:
  Solution(const bool is_displacement_function)
    : Function<dim>(dim)
    , is_displacement_function(is_displacement_function)
  {}

  double
  value(const Point<dim> &point, const unsigned int component) const
  {
    return std::sin(point[component] * 0.5 * numbers::PI) -
           (is_displacement_function ? point[component] : 0.0);
  }

private:
  const bool is_displacement_function;
};

static int counter = 0;

template <int dim, typename Fu>
void
do_test(const unsigned int degree,
        const unsigned int mapping_degree,
        const Fu          &fu,
        const bool         is_displacement_function)
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::MeshSmoothing::none,
    parallel::distributed::Triangulation<
      dim>::Settings::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 4);
  tria.refine_global(1);

  FESystem<dim>   fe(FE_Q<dim>(degree), dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  LinearAlgebra::distributed::Vector<double> vector(
    dof_handler.locally_owned_dofs(),
    locally_relevant_dofs,
    dof_handler.get_mpi_communicator());

  VectorTools::interpolate(dof_handler, fu, vector);

  {
    MappingQ<dim>      mapping(mapping_degree);
    MappingQCache<dim> mapping_cache(mapping_degree);
    mapping_cache.initialize(mapping,
                             dof_handler,
                             vector,
                             is_displacement_function);

    DataOut<dim> data_out;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    data_out.attach_triangulation(tria);

    data_out.build_patches(mapping_cache,
                           2,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);

#if false
    std::ofstream output("test." + std::to_string(counter++) + "." + std::to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)) + ".vtk");
    data_out.write_vtk(output);
#else
    data_out.write_vtk(deallog.get_file_stream());
#endif
  }

  {
    MGLevelObject<LinearAlgebra::distributed::Vector<double>> vectors(
      0, tria.n_global_levels() - 1);
    MGTransferMatrixFree<dim, double> transfer;
    transfer.build(dof_handler);
    transfer.interpolate_to_mg(dof_handler, vectors, vector);

    MappingQ<dim>      mapping(mapping_degree);
    MappingQCache<dim> mapping_cache(mapping_degree);
    mapping_cache.initialize(mapping,
                             dof_handler,
                             vectors,
                             is_displacement_function);

    const unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    for (unsigned int lvl = 0; lvl < tria.n_global_levels(); ++lvl)
      {
        DataOut<dim> data_out;

        DataOutBase::VtkFlags flags;
        flags.write_higher_order_cells = true;
        data_out.set_flags(flags);

        data_out.attach_triangulation(tria);

        data_out.set_cell_selection(
          [&](const typename Triangulation<dim>::cell_iterator &cell) {
            return (cell->level_subdomain_id() == rank) &&
                   (static_cast<unsigned int>(cell->level()) == lvl);
          });
        data_out.build_patches(mapping_cache,
                               2,
                               DataOut<dim>::curved_inner_cells);

#if true
        data_out.write_vtu_with_pvtu_record(
          "./", "mg_solution", lvl, MPI_COMM_WORLD, 1, 1);
#else
        data_out.write_vtk(deallog.get_file_stream());
#endif
      }
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  do_test<2>(3, 3, Solution<2>(true), true);
  do_test<2>(3, 3, Solution<2>(false), false);

  do_test<2>(3, 4, Solution<2>(true), true);
  do_test<2>(3, 4, Solution<2>(false), false);
}

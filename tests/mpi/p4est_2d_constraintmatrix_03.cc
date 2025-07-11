// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double>.distribute() for a distributed mesh
// with Trilinos
// Mesh: shell with random refinement

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test(const bool use_manifold_for_normal, const Mapping<dim> &mapping)
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  const double R0 = 6371000. - 2890000.;
  const double R1 = 6371000. - 35000.;


  GridGenerator::hyper_shell(tr, Point<dim>(), R0, R1, 12, true);
  GridTools::copy_boundary_to_manifold_id(tr);

  static SphericalManifold<dim> boundary;
  tr.set_manifold(0, boundary);
  tr.set_manifold(1, boundary);

  tr.refine_global(1);
  for (unsigned int step = 0; step < 20; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
                                                          tr.begin_active(),
                                                        endc = tr.end();

      for (; cell != endc; ++cell)
        if (Testing::rand() % 42 == 1)
          cell->set_refine_flag();

      tr.execute_coarsening_and_refinement();
    }

  DoFHandler<dim> dofh(tr);

  static FESystem<dim> fe(FE_Q<dim>(1 + 1), dim, FE_Q<dim>(1), 1);

  dofh.distribute_dofs(fe);

  const IndexSet &owned_set    = dofh.locally_owned_dofs();
  const IndexSet  dof_set      = DoFTools::extract_locally_active_dofs(dofh);
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);
  x = 2.0;

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);

  AffineConstraints<double> cm(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dofh, cm);
  ComponentMask velocity_mask(dim + 1, true);

  velocity_mask.set(dim, false);

  VectorTools::interpolate_boundary_values(
    dofh, 0, Functions::ZeroFunction<dim>(dim + 1), cm, velocity_mask);

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(1);


  VectorTools::compute_no_normal_flux_constraints(
    dofh, 0, no_normal_flux_boundaries, cm, mapping, use_manifold_for_normal);

  cm.close();

  cm.distribute(x);
  x_rel = x;

  TrilinosWrappers::MPI::Vector x_dub;
  x_dub.reinit(complete_index_set(dof_set.size()));
  x_dub.reinit(x_rel, false, true);

  if (myid == 0)
    x_dub.print(deallog.get_file_stream(), 8, true, false);

  tr.reset_manifold(0);
  tr.reset_manifold(1);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();
  deallog << "Normal with manifold:" << std::endl;
  test<2>(true, MappingQ1<2>());
  deallog << std::endl;

  deallog << "Normal with mapping:" << std::endl;
  test<2>(false, MappingQ<2>(3));
}

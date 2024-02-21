// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// create a shared tria mesh and refine it

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim, int spacedim>
void
write_mesh(const parallel::shared::Triangulation<dim, spacedim> &tria,
           const char                                           *filename_)
{
  DataOut<dim> data_out;
  data_out.attach_triangulation(tria);
  Vector<float> subdomain(tria.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = tria.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches();
  const std::string filename =
    (filename_ + Utilities::int_to_string(tria.locally_owned_subdomain(), 4));
  {
    std::ofstream output(filename + ".vtu");
    data_out.write_vtu(output);
  }
}



template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);

  AssertThrow(tr.with_artificial_cells() == false, ExcInternalError());

  const std::vector<unsigned int> &true_subdomain_ids_of_cells =
    tr.get_true_subdomain_ids_of_cells();

  AssertThrow(true_subdomain_ids_of_cells.size() == tr.n_active_cells(),
              ExcInternalError());


  GridGenerator::hyper_cube(tr);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  deallog << " locally_owned_subdomain(): " << tr.locally_owned_subdomain()
          << "\n"
          << " n_active_cells: " << tr.n_active_cells() << "\n"
          << " n_levels: " << tr.n_levels() << "\n"
          << " n_global_levels: " << tr.n_global_levels()
          << "\n"
          //<< " n_locally_owned_active_cells: " <<
          // tr.n_locally_owned_active_cells() << "\n"
          //<< " n_global_active_cells: " << tr.n_global_active_cells() << "\n"
          << std::endl;

  /*deallog << "n_locally_owned_active_cells_per_processor: ";
  std::vector<unsigned int> v = tr.n_locally_owned_active_cells_per_processor();
  for (unsigned int i=0;i<v.size();++i)
    deallog << v[i] << ' ';
    deallog << std::endl;*/

  // until parmetis is stable, do not output partitioning
  // deallog << "subdomains: ";
  typename parallel::shared::Triangulation<dim>::active_cell_iterator it =
    tr.begin_active();
  for (unsigned int index = 0; it != tr.end(); ++it, ++index)
    {
      // check that true subdomain_ids are the same as those, stored in
      // cell->subdomain_id()
      AssertThrow(true_subdomain_ids_of_cells[index] == it->subdomain_id(),
                  ExcInternalError());
      // deallog << it->subdomain_id() << ' ';
    }
  deallog << std::endl;

  // write_mesh(tr, "mesh");
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}

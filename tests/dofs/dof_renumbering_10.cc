// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test cell_wise with p::d::Tria

#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim, int fe_degree = 1>
void
test(const bool adaptive_ref = true)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_core =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);
  if (adaptive_ref)
    {
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->center().norm() < 0.5)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->center()[0] < 0.2)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
  else
    {
      tria.refine_global(1);
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  using CELL = typename DoFHandler<dim>::active_cell_iterator;
  std::vector<CELL> cell_order;
  for (auto &cell : dof.active_cell_iterators())
    if (cell->is_locally_owned())
      cell_order.push_back(cell);

  std::sort(cell_order.begin(),
            cell_order.end(),
            [](const CELL &a, const CELL &b) {
              std::vector<double> p1(dim), p2(dim);
              for (unsigned int d = 0; d < dim; ++d)
                {
                  p1[d] = a->center()[d];
                  p2[d] = b->center()[d];
                }
              return std::lexicographical_compare(p1.begin(),
                                                  p1.end(),
                                                  p2.begin(),
                                                  p2.end());
            });

  DoFRenumbering::cell_wise(dof, cell_order);

  deallog << std::endl << "cell order:" << std::endl;
  std::vector<types::global_dof_index> cell_dofs(fe.n_dofs_per_cell());
  for (const auto &c : cell_order)
    {
      c->get_active_or_mg_dof_indices(cell_dofs);
      deallog << c->center() << ':';
      for (const auto d : cell_dofs)
        deallog << ' ' << d;
      deallog << std::endl;
    }

  IndexSet owned_set = dof.locally_owned_dofs();

  // output in Gnuplot for visual inspection
  if (dim == 2)
    {
      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      DoFTools::map_dofs_to_support_points(mapping, dof, support_points);

      const std::string href = (adaptive_ref ? "" : "global_");
      const std::string base_filename =
        href + "grid" + dealii::Utilities::int_to_string(dim) + "_" +
        dealii::Utilities::int_to_string(fe_degree) + "_p" +
        dealii::Utilities::int_to_string(this_mpi_core);

      const std::string filename = base_filename + ".gp";
      std::ofstream     f(filename);

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
        << std::endl
        << "set output \"" << base_filename << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics" << std::endl
        << "unset ytics" << std::endl
        << "unset grid" << std::endl
        << "unset border" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels tc rgb 'red' nopoint notitle, '-' with labels point pt 4 offset 1,1 notitle"
        << std::endl;
      GridOut().write_gnuplot(tria, f);
      f << 'e' << std::endl;

      // output cell order
      for (unsigned int index = 0; index < cell_order.size(); ++index)
        f << cell_order[index]->center() << " \"" << index << "\"\n";

      f << std::flush;
      f << 'e' << std::endl << std::endl;

      // output only owned support points
      for (auto it = support_points.cbegin(); it != support_points.cend();)
        {
          if (owned_set.is_element(it->first))
            ++it;
          else
            support_points.erase(it++);
        }

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);

      f << 'e' << std::endl;
    }
}


int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  test<2>(false);

  return 0;
}

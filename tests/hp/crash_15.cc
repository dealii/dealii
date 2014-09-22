// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// in crash_14, we crash because we apparently don't unify all line dofs when
// we have FESystem(FE_Q(2),FE_DGQ(i)) for a bunch of different indices i
// coming together. note that all lines dofs are of the (always same) FE_Q

// apparently, what is happening is that we don't unify more than 2 finite
// elements on an edge in 3d, according to a comment at the top of
// hp::DoFHandler::compute_line_dof_identities at the time of this writing
// (and refering to a comment in the hp paper). there is now code that deals
// with the more narrow special case we have here

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (1);

  hp::FECollection<dim> fe_collection;
  for (unsigned int i=0; i<tria.n_active_cells(); ++i)
    fe_collection.push_back(FESystem<dim> (FE_Q<dim> (2), 1,
                                           FE_DGQ<dim> (i % 4), 1));

  hp::DoFHandler<dim> dof_handler(tria);

  unsigned int fe_index = 0;
  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell, ++fe_index)
    {
      deallog << "Setting fe_index=" << fe_index << " on cell " << cell
              << std::endl;
      cell->set_active_fe_index (fe_index);
    }

  dof_handler.distribute_dofs(fe_collection);

  std::vector<types::global_dof_index> line_dof_indices_1 (fe_collection[0].dofs_per_line +
                                                           2 * fe_collection[0].dofs_per_vertex);
  std::vector<types::global_dof_index> line_dof_indices_2 (fe_collection[0].dofs_per_line +
                                                           2 * fe_collection[0].dofs_per_vertex);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
      if (cell->line(line)->n_active_fe_indices() > 1)
        {
          deallog << "line=" << cell->line(line)
                  << std::endl;
          for (unsigned int i=0; i<cell->line(line)->n_active_fe_indices(); ++i)
            {
              deallog << "  active_fe_index="
                      << cell->line(line)->nth_active_fe_index(i)
                      << " ("
                      << fe_collection[cell->line(line)->nth_active_fe_index(i)].get_name()
                      << ")"
                      << std::endl;

              cell->line(line)->get_dof_indices (line_dof_indices_1,
                                                 cell->line(line)->nth_active_fe_index(i));
              deallog << "  dof indices=";
              for (unsigned int p=0; p<line_dof_indices_1.size(); ++p)
                deallog << line_dof_indices_1[p] << ' ';
              deallog << std::endl;
            }


          // if there are multiple active fe
          // indices, make sure that all their
          // fe indices were unified
          for (unsigned int i=0; i<cell->line(line)->n_active_fe_indices(); ++i)
            for (unsigned int j=i+1; j<cell->line(line)->n_active_fe_indices(); ++j)
              {
                cell->line(line)->get_dof_indices (line_dof_indices_1,
                                                   cell->line(line)->nth_active_fe_index(i));
                cell->line(line)->get_dof_indices (line_dof_indices_2,
                                                   cell->line(line)->nth_active_fe_index(j));

                Assert (line_dof_indices_1 == line_dof_indices_2,
                        ExcInternalError());
              }
        }
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// an extract of hp_constraints_q_system_x_01. something still goes wrong,
// akin to what happens in crash_15. (it turned out to be bogus index
// computations.)

char logname[] = "output";


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <vector>



template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
                                  FE_DGQ<dim>(j), 1));

  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  hp::DoFHandler<dim>        dof_handler(triangulation);

  // distribute fe_indices randomly
  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (Testing::rand() % fe.size());
  dof_handler.distribute_dofs (fe);

  // loop over all lines and make sure that
  // all the DoF indices on these lines are
  // identical
  std::vector<types::global_dof_index> indices_1;
  std::vector<types::global_dof_index> indices_2;

  std::set<unsigned int> line_already_treated;

  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
      if (line_already_treated.find (cell->line(l)->index())
          == line_already_treated.end())
        // line not yet treated
        {
          const typename hp::DoFHandler<dim>::active_line_iterator
          line = cell->line(l);
          deallog << "line=" << line << std::endl;
          line_already_treated.insert (line->index());

          for (unsigned int f=0; f<line->n_active_fe_indices(); ++f)
            {
              indices_1.resize (fe[line->nth_active_fe_index(f)].dofs_per_line +
                                2 * fe[line->nth_active_fe_index(f)].dofs_per_vertex);
              line->get_dof_indices (indices_1,
                                     line->nth_active_fe_index(f));

              deallog << "  fe index=" << line->nth_active_fe_index(f)
                      << ", indices=";
              for (unsigned int i=0; i<indices_1.size(); ++i)
                deallog << indices_1[i] << ' ';

              deallog << std::endl;
            }

          for (unsigned int f=0; f<line->n_active_fe_indices(); ++f)
            {
              indices_1.resize (fe[line->nth_active_fe_index(f)].dofs_per_line +
                                2 * fe[line->nth_active_fe_index(f)].dofs_per_vertex);
              line->get_dof_indices (indices_1,
                                     line->nth_active_fe_index(f));
              for (unsigned int g=f+1; g<line->n_active_fe_indices(); ++g)
                if (fe[line->nth_active_fe_index(f)].dofs_per_line
                    ==
                    fe[line->nth_active_fe_index(g)].dofs_per_line)
                  {
                    indices_2.resize (fe[line->nth_active_fe_index(g)].dofs_per_line +
                                      2 * fe[line->nth_active_fe_index(g)].dofs_per_vertex);
                    line->get_dof_indices (indices_2,
                                           line->nth_active_fe_index(g));
                    Assert (indices_1 == indices_2,
                            ExcInternalError());
                  }
            }
        }
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<3> ();
}


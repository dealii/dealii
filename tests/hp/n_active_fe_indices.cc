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



// record the output of n_active_fe_indices for cells, faces, and
// edges, as well as fe_index_is_active


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

#include <fstream>


template <int dim>
void check_cells (const hp::DoFHandler<dim> &dof_handler)
{
  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      deallog << "cell=" << cell << std::endl;
      deallog << "n="
              << cell->n_active_fe_indices ()
              << std::endl;
      deallog << "x=";
      for (unsigned int i=0; i<dof_handler.get_fe().size(); ++i)
        deallog << cell->fe_index_is_active(i);
      deallog << std::endl;

      Assert (cell->n_active_fe_indices () == 1, ExcInternalError());
    }
}


void check_faces (const hp::DoFHandler<1> &)
{}


template <int dim>
void check_faces (const hp::DoFHandler<dim> &dof_handler)
{
  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        deallog << "face=" << cell->face(f) << std::endl;
        deallog << "n="
                << cell->face(f)->n_active_fe_indices ()
                << std::endl;
        deallog << "x=";
        for (unsigned int i=0; i<dof_handler.get_fe().size(); ++i)
          deallog << cell->face(f)->fe_index_is_active(i);
        deallog << std::endl;

        Assert (cell->face(f)->n_active_fe_indices () >= 1, ExcInternalError());
        Assert (cell->face(f)->n_active_fe_indices () <= 2, ExcInternalError());
      }
}


void check_edges (const hp::DoFHandler<1> &)
{}


void check_edges (const hp::DoFHandler<2> &)
{}


template <int dim>
void check_edges (const hp::DoFHandler<dim> &dof_handler)
{
  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_cell; ++e)
      {
        deallog << "edge=" << cell->line(e) << std::endl;
        deallog << "n="
                << cell->line(e)->n_active_fe_indices ()
                << std::endl;
        deallog << "x=";
        for (unsigned int i=0; i<dof_handler.get_fe().size(); ++i)
          deallog << cell->line(e)->fe_index_is_active(i);
        deallog << std::endl;

        Assert (cell->line(e)->n_active_fe_indices () >= 1, ExcInternalError());
        Assert (cell->line(e)->n_active_fe_indices () <= dof_handler.get_fe().size(),
                ExcInternalError());
      }
}



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (2);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim> (1));
  fe_collection.push_back(FE_Q<dim> (2));
  fe_collection.push_back(FE_Q<dim> (3));
  fe_collection.push_back(FE_Q<dim> (4));

  hp::DoFHandler<dim> dof_handler(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (Testing::rand() % fe_collection.size());

  dof_handler.distribute_dofs(fe_collection);

  check_cells (dof_handler);
  check_faces (dof_handler);
  check_edges (dof_handler);
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}

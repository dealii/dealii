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



// when we have FESystem(FE_Q(k),FE_DGQ(1)) and FESystem(FE_Q(k),FE_DGQ(2)) on
// an interface, make sure all DoFs are unified. this is also tested in a
// slightly cleaner way by crash_15


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

  for (unsigned int k=1; k<4; ++k)
    {
      hp::FECollection<dim> fe_collection;
      for (unsigned int i=0; i<tria.n_active_cells(); ++i)
        fe_collection.push_back(FESystem<dim> (FE_Q<dim> (k), 1,
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

      std::vector<types::global_dof_index> face_dof_indices;
      std::vector<types::global_dof_index> neighbor_face_dof_indices;
      for (typename hp::DoFHandler<dim>::active_cell_iterator
           cell=dof_handler.begin_active();
           cell!=dof_handler.end(); ++cell)
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (!cell->at_boundary(face))
            {
              Assert (cell->get_fe().dofs_per_face
                      ==
                      cell->neighbor(face)->get_fe().dofs_per_face,
                      ExcInternalError());

              face_dof_indices.resize (cell->get_fe().dofs_per_face);
              neighbor_face_dof_indices.resize (cell->neighbor(face)
                                                ->get_fe().dofs_per_face);

              cell->face(face)->get_dof_indices (face_dof_indices,
                                                 cell->active_fe_index());
              cell->face(face)->get_dof_indices (neighbor_face_dof_indices,
                                                 cell->neighbor(face)
                                                 ->active_fe_index());

              deallog << "cell=" << cell << ", face=" << face << std::endl;
              deallog << "fe1=" << cell->get_fe().get_name()
                      << ", fe2=" << cell->neighbor(face)->get_fe().get_name()
                      << std::endl;

              for (unsigned int i=0; i<face_dof_indices.size(); ++i)
                {
                  deallog << face_dof_indices[i] << std::endl;

                  Assert (face_dof_indices[i] ==
                          neighbor_face_dof_indices[i],
                          ExcInternalError());
                }

              deallog << std::endl;
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

// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// check serialization for DoFHandler<1,dim>

#include "serialization.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

namespace dealii
{
  template <int dim, int spacedim>
  bool operator == (const DoFHandler<dim,spacedim> &t1,
                    const DoFHandler<dim,spacedim> &t2)
  {
    // test a few attributes, though we can't
    // test everything unfortunately...
    typename DoFHandler<dim,spacedim>::cell_iterator
    c1 = t1.begin(),
    c2 = t2.begin();
    for (; (c1 != t1.end()) && (c2 != t2.end()); ++c1, ++c2)
      {
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            if (c1->vertex(v) != c2->vertex(v))
              return false;
            if (c1->vertex_index(v) != c2->vertex_index(v))
              return false;
          }

        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          {
            if (c1->face(f)->at_boundary() != c2->face(f)->at_boundary())
              return false;

            if (c1->face(f)->at_boundary())
              {
                if (c1->face(f)->boundary_indicator() !=
                    c2->face(f)->boundary_indicator())
                  return false;
              }
            else
              {
                if (c1->neighbor(f)->level() != c2->neighbor(f)->level())
                  return false;
                if (c1->neighbor(f)->index() != c2->neighbor(f)->index())
                  return false;
              }
          }

        if (c1->active() && c2->active()
            && (c1->subdomain_id() != c2->subdomain_id()))
          return false;

        if (c1->material_id() != c2->material_id())
          return false;

        if (c1->user_index() != c2->user_index())
          return false;

        if (c1->user_flag_set() != c2->user_flag_set())
          return false;

        if (c1->get_fe().get_name() != c2->get_fe().get_name())
          return false;

        if (c1->active_fe_index() != c2->active_fe_index())
          return false;

        // compare dofs on this cell and then on the faces
        if (c1->has_children() == false)
          {
            std::vector<types::global_dof_index> local_dofs_1 (c1->get_fe().dofs_per_cell);
            std::vector<types::global_dof_index> local_dofs_2 (c2->get_fe().dofs_per_cell);

            c1->get_dof_indices (local_dofs_1);
            c2->get_dof_indices (local_dofs_2);
            if (local_dofs_1 != local_dofs_2)
              return false;

            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                std::vector<types::global_dof_index> local_dofs_1 (c1->get_fe().dofs_per_face);
                std::vector<types::global_dof_index> local_dofs_2 (c2->get_fe().dofs_per_face);

                c1->face(f)->get_dof_indices (local_dofs_1);
                c2->face(f)->get_dof_indices (local_dofs_2);
                if (local_dofs_1 != local_dofs_2)
                  return false;
              }
          }
      }

    // also check the order of raw iterators as they contain
    // something about the history of the triangulation
    typename DoFHandler<dim,spacedim>::cell_iterator
    r1 = t1.begin(),
    r2 = t2.begin();
    for (; (r1 != t1.end()) && (r2 != t2.end()); ++r1, ++r2)
      {
        if (r1->level() != r2->level())
          return false;
        if (r1->index() != r2->index())
          return false;
      }

    return true;
  }
}

template <int dim, int spacedim>
void do_boundary (Triangulation<dim,spacedim> &t1)
{
  typename Triangulation<dim,spacedim>::cell_iterator
  c1 = t1.begin();
  for (; c1 != t1.end(); ++c1)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (c1->at_boundary(f))
        c1->face(f)->set_boundary_indicator (42);
}


template <int spacedim>
void do_boundary (Triangulation<1,spacedim> &)
{}


template <int dim, int spacedim>
void test ()
{
  Triangulation<dim,spacedim> tria;

  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);
  tria.begin_active()->set_subdomain_id (1);
  tria.begin_active()->set_material_id (2);
  tria.begin_active()->set_user_index (3);
  tria.begin_active()->set_user_flag ();
  tria.begin_active()->set_refine_flag (RefinementCase<dim>::cut_x);

  do_boundary (tria);

  FESystem<dim,spacedim> fe (FE_Q<dim,spacedim>(2),dim,
                             FE_Q<dim,spacedim>(1),1);

  DoFHandler<dim,spacedim> dof_1 (tria);
  DoFHandler<dim,spacedim> dof_2 (tria);

  dof_1.distribute_dofs (fe);
  dof_2.distribute_dofs (fe);

  // right now, both DoFHandlers are the same. Renumber one of the them
  DoFRenumbering::Cuthill_McKee (dof_1);

  verify (dof_1, dof_2);
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1,1> ();
  test<1,2> ();
  test<2,2> ();
  test<2,3> ();
  test<3,3> ();

  deallog << "OK" << std::endl;
}

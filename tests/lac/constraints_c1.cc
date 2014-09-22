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

/*
 * Author: Guido Kanschat, 2010, 2012
 */

// H2-elements implemented through constraints on the degrees of
// freedom. After adding all constrained lines, the program hangs in
// constraints.close().
// This is because it contains a cycle. The test now checks that this
// is detected.


#include <deal.II/base/job_identifier.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace dealii;

template <int dim>
void
setup_constraints(const DoFHandler<dim> &dof_handler)
{
  ConstraintMatrix     constraints;
  constraints.clear();
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  // Set up derivative constraints to
  // make element C1
  std::vector<std::vector<std::vector<double> > > vertex_constraints(
    GeometryInfo<dim>::vertices_per_cell, std::vector<std::vector<double> >(
      dim+1, std::vector<double>(fe.dofs_per_cell, 0.)));

  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    {
      Point<dim> v;
      for (unsigned int d=0; d<dim; ++d)
        {
          const unsigned int ds = 1<<d;
          const unsigned int id = vertex/ds;
          v[d] = (id%2 != 0) ? 1. : 0.;
        }

      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
          for (unsigned int d=0; d<dim; ++d)
            {
              const double f = fe.shape_grad(j,v)[d];
              if (std::fabs(f) > 1.e-12)
                vertex_constraints[vertex][d][j] = f;
            }
          const double f = fe.shape_grad_grad(j,v)[0][1];
          if (std::fabs(f) > 1.e-12)
            vertex_constraints[vertex][dim][j] = f;
        }
    }

  std::vector<types::global_dof_index> cell_indices(fe.dofs_per_cell);
  std::vector<types::global_dof_index> neighbor_indices(fe.dofs_per_cell);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      deallog << "New cell" << std::endl;

      cell->get_dof_indices(cell_indices);

      // Do lower left and upper
      // right vertex
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; face+=2)
        {
          if (cell->at_boundary(face)) continue;
          for (unsigned int fvertex=0; fvertex<GeometryInfo<dim>::vertices_per_face; ++fvertex)
            {
              const unsigned int other_face = (face==0) ? (2+fvertex) : fvertex;
              if (cell->at_boundary(other_face)) continue;

              const unsigned int vertex = GeometryInfo<dim>::face_to_cell_vertices(face, fvertex);
              typename DoFHandler<dim>::active_cell_iterator neighbor = cell->neighbor(other_face);
              const unsigned int neighbor_face = cell->neighbor_of_neighbor(other_face);
              neighbor->get_dof_indices(neighbor_indices);
              unsigned int neighbor_vertex = GeometryInfo<dim>::face_to_cell_vertices(neighbor_face, fvertex);

              const unsigned int d = GeometryInfo<dim>::unit_normal_direction[other_face];

              std::vector<std::pair<types::global_dof_index, double> > rhs;
              const unsigned int constrained = fe.face_to_cell_index(
                                                 (fvertex == 0) ? 2 : fe.dofs_per_face-1, face);
              double constrained_weight = 0.;

              for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                {
                  if (i == constrained)
                    constrained_weight = vertex_constraints[vertex][d][i];
                  else if (vertex_constraints[vertex][d][i] != 0.)
                    rhs.push_back(std::make_pair(cell_indices[i], -vertex_constraints[vertex][d][i]));
                  if (vertex_constraints[neighbor_vertex][d][i] != 0.)
                    rhs.push_back(std::make_pair(neighbor_indices[i], vertex_constraints[neighbor_vertex][d][i]));
                }
              deallog << " v" << vertex
                      << " f" << face
                      << " o" << other_face
                      << " d" << d
                      << " l" << constrained
                      << " g" << cell_indices[constrained]
//          << " w " << constrained_weight
                      << " rhs ";
              for (unsigned int i=0; i<rhs.size(); ++i)
                deallog << ' ' << rhs[i].first;
              deallog << std::endl;
              for (unsigned int i=0; i<rhs.size(); ++i)
                rhs[i].second /= constrained_weight;
              constraints.add_line(cell_indices[constrained]);
              constraints.add_entries(cell_indices[constrained], rhs);
            }
        }
    }

  deallog << "Closing" << std::endl;
  try
    {
      constraints.close();
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
  deallog << "Closed" << std::endl;
}


template <int dim>
void
run(const FiniteElement<dim> &fe)
{
  deallog << "Element: " << fe.get_name() << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global(2);
  deallog << "Triangulation "
          << triangulation.n_active_cells() << " cells, "
          << triangulation.n_levels() << " levels" << std::endl;
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  setup_constraints (dof_handler);
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  FE_Q<2> fe(3);
  run(fe);
}

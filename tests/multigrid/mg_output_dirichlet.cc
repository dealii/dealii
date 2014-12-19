// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// Check the results of MGTransferPrebuilt with Dirichlet boundary conditions

//TODO:[GK] Add checks for RT again!

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim, typename number, int spacedim>
void
reinit_vector (const dealii::DoFHandler<dim,spacedim> &mg_dof,
               MGLevelObject<dealii::Vector<number> > &v)
{
  for (unsigned int level=v.min_level();
       level<=v.max_level(); ++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }
}


template<int dim>
void refine_mesh (Triangulation<dim> &triangulation)
{
  bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    {
      for (unsigned int vertex=0;
           vertex < GeometryInfo<dim>::vertices_per_cell;
           ++vertex)
        {
          const Point<dim> p = cell->vertex(vertex);
          const Point<dim> origin = (dim == 2 ?
                                     Point<dim>(0,0) :
                                     Point<dim>(0,0,0));
          const double dist = p.distance(origin);
          if (dist<0.25/numbers::PI)
            {
              cell->set_refine_flag ();
              cell_refined = true;
              break;
            }
        }
    }
  if (!cell_refined) //if no cell was selected for refinement, refine global
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement ();
}

template <int dim>
void initialize (const DoFHandler<dim> &dof,
                 Vector<double> &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dof.begin_active();
       cell != dof.end(); ++cell)
    {
      cell->get_dof_indices(dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        u(dof_indices[i]) = ++counter;
    }
}

template <int dim>
void initialize (const DoFHandler<dim> &dof,
                 MGLevelObject<Vector<double> > &u)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  const unsigned int dofs_per_face = dof.get_fe().dofs_per_face;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> face_indices(dofs_per_face);

  for (unsigned int l=0; l<dof.get_tria().n_levels(); ++l)
    {
      for (typename DoFHandler<dim>::cell_iterator
           cell = dof.begin_mg(l);
           cell != dof.end_mg(l); ++cell)
        {
          cell->get_mg_dof_indices(dof_indices);
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            {
              cell->face(f)->get_mg_dof_indices(cell->level(),face_indices);
              if (cell->face(f)->at_boundary())
                {
                  for (unsigned int i=0; i<dofs_per_face; ++i)
                    u[cell->level()](face_indices[i]) = 1;
                }
              else
                {
                  for (unsigned int i=0; i<dofs_per_face; ++i)
                    if (u[cell->level()](face_indices[i]) != 1)
                      u[cell->level()](face_indices[i]) = 0;
                }
            }
        }
    }
}

template <int dim>
void print (const DoFHandler<dim> &dof,
            MGLevelObject<Vector<double> > &u)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (unsigned int l=0; l<dof.get_tria().n_levels(); ++l)
    {
      deallog << std::endl;
      deallog << "Level " << l << std::endl;
      for (typename DoFHandler<dim>::cell_iterator
           cell = dof.begin_mg(l);
           cell != dof.end_mg(l); ++cell)
        {
          cell->get_mg_dof_indices(dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            deallog << ' ' << (int)u[l](dof_indices[i]);
        }
    }
}

template <int dim>
void print_diff (const DoFHandler<dim> &dof_1, const DoFHandler<dim> &dof_2,
                 const Vector<double> &u, const Vector<double> &v)
{
  Vector<double> diff;
  diff.reinit (u);
  const unsigned int dofs_per_cell = dof_1.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices_1(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_2(dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator
       cell_1 = dof_1.begin_active(), cell_2 = dof_2.begin_active();
       cell_1 != dof_1.end(); ++cell_1, ++cell_2)
    {
      cell_1->get_dof_indices(dof_indices_1);
      cell_2->get_dof_indices(dof_indices_2);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        diff(dof_indices_1[i]) = u(dof_indices_1[i]) - v(dof_indices_2[i]);
    }
  deallog << std::endl;
  deallog << "diff " << diff.l2_norm() << std::endl;
}

template <int dim>
void check_simple(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr,-1,1);
  typename FunctionMap<dim>::type      dirichlet_boundary_functions;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary_functions[0] = &homogeneous_dirichlet_bc;

  tr.refine_global (1);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs(fe);

  DoFHandler<dim> mgdof_renumbered(tr);
  mgdof_renumbered.distribute_dofs(fe);
  mgdof_renumbered.distribute_mg_dofs(fe);

  ConstraintMatrix     hnc, hnc_renumbered;
  hnc.clear ();
  hnc_renumbered.clear ();
  DoFTools::make_hanging_node_constraints (
    mgdof,
    hnc);

  DoFTools::make_hanging_node_constraints (
    mgdof_renumbered,
    hnc_renumbered);
  hnc.close ();
  hnc_renumbered.close ();


  std::vector<unsigned int> block_component (4,0);
  block_component[2] = 1;
  block_component[3] = 1;

  DoFRenumbering::component_wise (mgdof_renumbered, block_component);
  for (unsigned int level=0; level<tr.n_levels(); ++level)
    DoFRenumbering::component_wise (mgdof_renumbered, level, block_component);


  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mgdof, dirichlet_boundary_functions);

  MGConstrainedDoFs mg_constrained_dofs_renumbered;
  mg_constrained_dofs_renumbered.initialize(mgdof_renumbered, dirichlet_boundary_functions);

  MGTransferPrebuilt<Vector<double> > transfer(hnc, mg_constrained_dofs);
  transfer.build_matrices(mgdof);

  MGTransferPrebuilt<Vector<double> > transfer_renumbered(hnc_renumbered, mg_constrained_dofs_renumbered);
  transfer_renumbered.build_matrices(mgdof_renumbered);

  Vector<double> u(mgdof.n_dofs());
  initialize(mgdof,u);

  MGLevelObject<Vector<double> > v(0, tr.n_levels()-1);
  reinit_vector(mgdof, v);

  transfer.copy_to_mg(mgdof, v, u);
  print(mgdof, v);

  u=0;
  initialize(mgdof, v);
  for (unsigned int l=0; l<tr.n_levels()-1; ++l)
    transfer.prolongate (l+1, v[l+1], v[l]);

  transfer.copy_from_mg(mgdof, u, v);
  Vector<double> diff = u;

  initialize(mgdof_renumbered,u);
  reinit_vector(mgdof_renumbered, v);
  transfer_renumbered.copy_to_mg(mgdof_renumbered, v, u);
  deallog << "copy_to_mg" << std::endl;
  print(mgdof_renumbered, v);

  u=0;
  initialize(mgdof_renumbered, v);
  for (unsigned int l=0; l<tr.n_levels()-1; ++l)
    transfer_renumbered.prolongate (l+1, v[l+1], v[l]);

  transfer_renumbered.copy_from_mg(mgdof_renumbered, u, v);
  print_diff(mgdof_renumbered, mgdof, u, diff);
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  //check_simple (FESystem<2>(FE_Q<2>(1), 2));
  check_simple (FESystem<2>(FE_Q<2>(1), 4));
}

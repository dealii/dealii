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


// like _01 but on adaptively refined grid


//TODO:[GK] Add checks for RT again!
#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
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
#include <deal.II/multigrid/mg_dof_handler.h>
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
reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
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
void initialize (const MGDoFHandler<dim> &dof,
                 Vector<double> &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
       cell = dof.begin_active();
       cell != dof.end(); ++cell)
    {
      cell->get_dof_indices(dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        u(dof_indices[i]) = ++counter;
    }
}

template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
                 MGLevelObject<Vector<double> > &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator
  cell = dof.begin(0);
  cell->get_mg_dof_indices(dof_indices);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    u[0](dof_indices[i]) = ++counter;
}

template <int dim>
void print (const MGDoFHandler<dim> &dof,
            MGLevelObject<Vector<double> > &u)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (unsigned int l=0; l<dof.get_tria().n_levels(); ++l)
    {
      deallog << std::endl;
      deallog << "Level " << l << std::endl;
      for (typename MGDoFHandler<dim>::cell_iterator
           cell = dof.begin(l);
           cell != dof.end(l); ++cell)
        {
          cell->get_mg_dof_indices(dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            deallog << ' ' << (int)u[l](dof_indices[i]);
        }
    }
}

template <int dim>
void print_diff (const MGDoFHandler<dim> &dof_1, const MGDoFHandler<dim> &dof_2,
                 const Vector<double> &u, const Vector<double> &v)
{
  Vector<double> diff;
  diff.reinit (u);
  const unsigned int dofs_per_cell = dof_1.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices_1(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_2(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
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
  tr.refine_global (1);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);
  refine_mesh(tr);


  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);

  MGDoFHandler<dim> mgdof_renumbered(tr);
  mgdof_renumbered.distribute_dofs(fe);

  std::vector<unsigned int> block_component (4,0);
  block_component[2] = 1;
  block_component[3] = 1;

  DoFRenumbering::component_wise (mgdof_renumbered, block_component);
  for (unsigned int level=0; level<tr.n_levels(); ++level)
    DoFRenumbering::component_wise (mgdof_renumbered, level, block_component);

  MGTransferPrebuilt<Vector<double> > transfer;
  transfer.build_matrices(mgdof);

  MGTransferPrebuilt<Vector<double> > transfer_renumbered;
  transfer_renumbered.build_matrices(mgdof_renumbered);

  Vector<double> u(mgdof.n_dofs());
  initialize(mgdof,u);
  deallog << "initialize " << u.l2_norm() << std::endl;

  MGLevelObject<Vector<double> > v(0, tr.n_levels()-1);
  reinit_vector(mgdof, v);

  transfer.copy_to_mg(mgdof, v, u);
  //print(mgdof, v);

  u=0;
  initialize(mgdof, v);
  for (unsigned int l=0; l<tr.n_levels()-1; ++l)
    transfer.prolongate (l+1, v[l+1], v[l]);

  transfer.copy_from_mg(mgdof, u, v);
  Vector<double> diff = u;

  initialize(mgdof_renumbered,u);
  deallog << "initialize " << u.l2_norm() << std::endl;

  reinit_vector(mgdof_renumbered, v);
  transfer_renumbered.copy_to_mg(mgdof_renumbered, v, u);
  deallog << "copy_to_mg" << std::endl;
  //print(mgdof_renumbered, v);

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

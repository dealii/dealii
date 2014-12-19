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


// Compare MGTransferPrebuilt, MGTransferBlock and MGTransferSelect

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_block.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>
#include <numeric>

using namespace std;


template <int dim, typename number, int spacedim>
void
reinit_vector_by_blocks (
  const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
  MGLevelObject<dealii::Vector<number> > &v,
  const unsigned int selected_block,
  std::vector<std::vector<types::global_dof_index> > &ndofs)
{
  const unsigned int n_blocks = mg_dof.get_fe().n_blocks();
  Assert(selected_block < n_blocks, ExcIndexRange(selected_block, 0, n_blocks));

  std::vector<bool> selected(n_blocks, false);
  selected[selected_block] = true;

  if (ndofs.size() == 0)
    {
      std::vector<std::vector<types::global_dof_index> >
      new_dofs(mg_dof.get_tria().n_levels(),
               std::vector<types::global_dof_index>(selected.size()));
      std::swap(ndofs, new_dofs);
      MGTools::count_dofs_per_block (mg_dof, ndofs);
    }

  for (unsigned int level=v.min_level();
       level<=v.max_level(); ++level)
    {
      v[level].reinit(ndofs[level][selected_block]);
    }
}


template <int dim, typename number, int spacedim>
void
reinit_vector_by_blocks (
  const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
  MGLevelObject<BlockVector<number> > &v,
  const std::vector<bool> &sel,
  std::vector<std::vector<types::global_dof_index> > &ndofs)
{
  std::vector<bool> selected=sel;
  // Compute the number of blocks needed
  const unsigned int n_selected
    = std::accumulate(selected.begin(),
                      selected.end(),
                      0U);

  if (ndofs.size() == 0)
    {
      std::vector<std::vector<types::global_dof_index> >
      new_dofs(mg_dof.get_tria().n_levels(),
               std::vector<types::global_dof_index>(selected.size()));
      std::swap(ndofs, new_dofs);
      MGTools::count_dofs_per_block (mg_dof, ndofs);
    }

  for (unsigned int level=v.min_level();
       level<=v.max_level(); ++level)
    {
      v[level].reinit(n_selected, 0);
      unsigned int k=0;
      for (unsigned int i=0; i<selected.size() && (k<v[level].n_blocks()); ++i)
        {
          if (selected[i])
            {
              v[level].block(k++).reinit(ndofs[level][i]);
            }
          v[level].collect_sizes();
        }
    }
}




template <int dim>
void check_block(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  std::vector<bool> selected(fe.n_blocks(), true);

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MGDoFHandler<dim> mgdof(tr);
  DoFHandler<dim> &dof=mgdof;
  mgdof.distribute_dofs(fe);

  // Make sure all orderings are the
  // same
  Point<dim> direction;
  for (unsigned int d=0; d<dim; ++d)
    direction[d] = d*d*d;
  DoFRenumbering::downstream(dof, direction);
  DoFRenumbering::component_wise(dof);
  for (unsigned int l=0; l<tr.n_levels(); ++l)
    {
      DoFRenumbering::downstream(mgdof, l, direction);
      DoFRenumbering::component_wise(mgdof, l);
    }

  // Store sizes
  vector<types::global_dof_index> ndofs(fe.n_blocks());
  DoFTools::count_dofs_per_block(mgdof, ndofs);
  std::vector<std::vector<types::global_dof_index> > mg_ndofs(mgdof.get_tria().n_levels(),
                                                              std::vector<types::global_dof_index>(fe.n_blocks()));
  MGTools::count_dofs_per_block(mgdof, mg_ndofs);

  MGTransferPrebuilt<BlockVector<double> > transfer;
  MGTransferBlock<double> transfer_block;
  MGTransferBlockSelect<double> transfer_select;
  transfer.build_matrices(mgdof);
  transfer_block.build_matrices(dof, mgdof, selected);
  transfer_select.build_matrices(dof, mgdof, 0);

  BlockVector<double> u2(mg_ndofs[2]);
  BlockVector<double> u1(mg_ndofs[1]);
  BlockVector<double> u0(mg_ndofs[0]);
  BlockVector<double> v2(mg_ndofs[2]);
  BlockVector<double> v1(mg_ndofs[1]);
  BlockVector<double> v0(mg_ndofs[0]);

  // Prolongate a constant function
  // twice
  u0 = 1;
  transfer.prolongate(1,u1,u0);
  transfer.prolongate(2,u2,u1);
  v0 = 1;
  transfer_block.prolongate(1,v1,v0);
  transfer_block.prolongate(2,v2,v1);
  v0.add(-1., u0);
  v1.add(-1., u1);
  v2.add(-1., u2);
  // These outputs are just the
  // number of dofs on each level
  deallog << "Prolongate " << v0.l2_norm()
          << ' ' << v1.l2_norm()
          << ' ' << v2.l2_norm()
          << std::endl;

  v0 = 1.;
  transfer_select.prolongate(1,v1.block(0),v0.block(0));
  transfer_select.prolongate(2,v2.block(0),v1.block(0));
  v0.add(-1., u0);
  v1.add(-1., u1);
  v2.add(-1., u2);
  deallog << "Select     " << v0.block(0).l2_norm()
          << ' ' << v1.block(0).l2_norm()
          << ' ' << v2.block(0).l2_norm()
          << std::endl;

  v2 = u2;
  u1 = 0.;
  u0 = 0.;
  transfer.restrict_and_add(2,u1,u2);
  transfer.restrict_and_add(1,u0,u1);
  v1 = 0.;
  v0 = 0.;
  transfer_block.restrict_and_add(2,v1,v2);
  transfer_block.restrict_and_add(1,v0,v1);
  v0.add(-1., u0);
  v1.add(-1., u1);
  deallog << "Restrict " << v0.l2_norm()
          << ' ' << v1.l2_norm()
          << std::endl;

  v1 = 0.;
  v0 = 0.;
  transfer_select.restrict_and_add(2,v1.block(0),v2.block(0));
  transfer_select.restrict_and_add(1,v0.block(0),v1.block(0));
  v0.add(-1., u0);
  v1.add(-1., u1);
  deallog << "Select   " << v0.block(0).l2_norm()
          << ' ' << v1.block(0).l2_norm()
          << std::endl;


  // Check copy to mg and back
  BlockVector<double> u;
  u.reinit (ndofs);
  for (unsigned int i=0; i<u.size(); ++i)
    u(i) = i+1;

  std::vector<std::vector<types::global_dof_index> > cached_sizes;
  MGLevelObject<BlockVector<double> > v;
  MGLevelObject<BlockVector<double> > wb;
  MGLevelObject<Vector<double> > ws;
  v.resize(0, tr.n_levels()-1);
  wb.resize(0, tr.n_levels()-1);
  ws.resize(0, tr.n_levels()-1);
  reinit_vector_by_blocks(mgdof, v, selected, cached_sizes);
  reinit_vector_by_blocks(mgdof, wb, selected, cached_sizes);
  reinit_vector_by_blocks(mgdof, ws, 0, cached_sizes);

  deallog << "copy to mg";
  transfer.copy_to_mg(mgdof, v, u);
  transfer_block.copy_to_mg(mgdof, wb, u);
  transfer_select.copy_to_mg(mgdof, ws, u);

  for (unsigned int l=0; l<tr.n_levels(); ++l)
    {
      wb[l].add(-1., v[l]);
      ws[l].add(-1., v[l].block(0));
      deallog << ' ' << wb[l].l2_norm();
      deallog << ' ' << ws[l].l2_norm();
    }
  deallog << std::endl << "copy from mg ";
  // Now do the opposite: fill a
  // multigrid vector counting the
  // dofs and see where the numbers go
  u = 0.;
  for (unsigned int i=0; i<v[2].size(); ++i)
    {
      v[2](i) = i+1;
      wb[2](i) = i+1;
    }
  for (unsigned int i=0; i<ws[2].size(); ++i)
    ws[2](i) = i+1;
  BlockVector<double> uu;
  uu.reinit(u);
  transfer.copy_from_mg_add(mgdof, u, v);
  transfer_block.copy_from_mg_add(mgdof, uu, v);
  u.add(-1., uu);
  deallog << u.l2_norm() << std::endl;

}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  std::vector<double> factors;

  FE_DGQ<2> q0(0);
  FE_DGQ<2> q1(1);
  FE_Q<2> cq1(1);
  FE_RaviartThomasNodal<2> rt0(0);
  FE_RaviartThomasNodal<2> rt1(1);

  FESystem<2> fe0(rt1, 1, q1, 1);
  FESystem<2> fe1(rt0, 2, q0, 2);

  check_block(q0);
  check_block(cq1);
  check_block(rt1);

  check_block(fe0);
  check_block(fe1);
}

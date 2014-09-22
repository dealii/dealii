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
void check_block(const FiniteElement<dim> &fe,
                 const vector<bool> &selected,
                 const vector<double> &factors)
{
  deallog << fe.get_name() << std::endl << "selected ";
  for (unsigned int i=0; i<selected.size(); ++i)
    if (selected[i])
      deallog << ' ' << i;
  deallog << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MGDoFHandler<dim> mgdof(tr);
  DoFHandler<dim> &dof=mgdof;
  mgdof.distribute_dofs(fe);
  DoFRenumbering::component_wise(mgdof);
  vector<types::global_dof_index> ndofs(fe.n_blocks());
  DoFTools::count_dofs_per_block(mgdof, ndofs);

  for (unsigned int l=0; l<tr.n_levels(); ++l)
    DoFRenumbering::component_wise(mgdof, l);
  std::vector<std::vector<types::global_dof_index> > mg_ndofs(mgdof.get_tria().n_levels(),
                                                              std::vector<types::global_dof_index>(fe.n_blocks()));
  MGTools::count_dofs_per_block(mgdof, mg_ndofs);

  deallog << "Global  dofs:";
  for (unsigned int i=0; i<ndofs.size(); ++i)
    deallog << ' ' << ndofs[i];
  deallog << std::endl;
  for (unsigned int l=0; l<mg_ndofs.size(); ++l)
    {
      deallog << "Level " << l << " dofs:";
      for (unsigned int i=0; i<mg_ndofs[l].size(); ++i)
        deallog << ' ' << mg_ndofs[l][i];
      deallog << std::endl;
    }

  PrimitiveVectorMemory<Vector<double> > mem;
  MGTransferBlock<double> transfer;
  transfer.build_matrices(dof, mgdof, selected);
  if (factors.size()>0)
    transfer.initialize(factors, mem);

  MGLevelObject< BlockVector<double> > u(0, tr.n_levels()-1);

  reinit_vector_by_blocks(mgdof, u, selected, mg_ndofs);

  // Prolongate a constant function
  // twice
  u[0] = 1;
  transfer.prolongate(1,u[1],u[0]);
  transfer.prolongate(2,u[2],u[1]);
  // These outputs are just the
  // number of dofs on each level
  deallog << "u0";
  for (unsigned int b=0; b<u[0].n_blocks(); ++b)
    deallog << '\t' << (int) (u[0].block(b)*u[0].block(b)+.4);
  deallog << std::endl << "u1";
  for (unsigned int b=0; b<u[1].n_blocks(); ++b)
    deallog << '\t' << (int) (u[1].block(b)*u[1].block(b)+.4);
  deallog << std::endl << "u2";
  for (unsigned int b=0; b<u[2].n_blocks(); ++b)
    deallog << '\t' << (int) (u[2].block(b)*u[2].block(b)+.4);
  deallog << std::endl;

  u[1] = 0.;
  u[0] = 0.;
  transfer.restrict_and_add(2,u[1],u[2]);
  transfer.restrict_and_add(1,u[0],u[1]);
  // After adding the restrictions,
  // things get bigger.
  deallog << "u1";
  for (unsigned int b=0; b<u[1].n_blocks(); ++b)
    deallog << '\t' << (int) (u[1].block(b)*u[1].block(b)+.5);
  deallog << std::endl << "u0";
  for (unsigned int b=0; b<u[0].n_blocks(); ++b)
    deallog << '\t' << (int) (u[0].block(b)*u[0].block(b)+.5);
  deallog << std::endl;

  // Check copy to mg and back
  // Fill a global vector by counting
  // from one up
  BlockVector<double> v;
  v.reinit (ndofs);
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i+1;

  // See what part gets copied to mg
  u.resize(0, tr.n_levels()-1);
  reinit_vector_by_blocks(mgdof, u, selected, mg_ndofs);

  transfer.copy_to_mg(mgdof, u, v);
  for (unsigned int i=0; i<u[2].size(); ++i)
    deallog << ' ' << (int) u[2](i);
  deallog << std::endl;

  // Now do the opposite: fill a
  // multigrid vector counting the
  // dofs and see where the numbers go
  for (unsigned int i=0; i<u[2].size(); ++i)
    u[2](i) = i+1;
  v = 0.;
  transfer.copy_from_mg(mgdof, v, u);
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << (int) v(i);
  deallog << std::endl;
  v.equ(-1., v);
  transfer.copy_from_mg_add(mgdof, v, u);
  deallog << "diff " << v.l2_norm() << std::endl;
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
  FE_RaviartThomasNodal<2> rt0(0);
  FE_RaviartThomasNodal<2> rt1(1);

  FESystem<2> fe0(rt1, 1, q1, 1);
  FESystem<2> fe1(rt0, 2, q0, 2);

  vector<bool> s1(2, true);
  deallog << "All" << std::endl;
  check_block(fe0, s1, factors);

  s1[1] = false;
  deallog << "Velocity" << std::endl;
  check_block(fe0, s1, factors);

  s1[1] = true;
  s1[0] = false;
  deallog << "Pressure" << std::endl;
  check_block(fe0, s1, factors);

  s1.resize(4,true);
  s1[0] = false;
  s1[2] = false;
  check_block(fe1, s1, factors);
}

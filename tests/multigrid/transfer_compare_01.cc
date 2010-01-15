//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 - 2007, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// Compare MGTransferPrebuilt, MGTransferBlock and MGTransferSelect

#include "../tests.h"
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgq.h>
#include <fe/fe_q.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_transfer_block.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_level_object.h>

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
  std::vector<std::vector<unsigned int> >& ndofs)
{
  std::vector<bool> selected=sel;
				   // Compute the number of blocks needed
  const unsigned int n_selected
    = std::accumulate(selected.begin(),
		      selected.end(),
		      0U);

  if (ndofs.size() == 0)
    {
      std::vector<std::vector<unsigned int> >
	new_dofs(mg_dof.get_tria().n_levels(),
		 std::vector<unsigned int>(selected.size()));
      std::swap(ndofs, new_dofs);
      MGTools::count_dofs_per_block (mg_dof, ndofs);
    }

  for (unsigned int level=v.get_minlevel();
       level<=v.get_maxlevel();++level)
    {
      v[level].reinit(n_selected, 0);
      unsigned int k=0;
      for (unsigned int i=0;i<selected.size() && (k<v[level].n_blocks());++i)
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
void check_block(const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  std::vector<bool> selected(fe.n_blocks(), true);

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MGDoFHandler<dim> mgdof(tr);
  DoFHandler<dim>& dof=mgdof;
  mgdof.distribute_dofs(fe);

				   // Make sure all orderings are the
				   // same
  Point<dim> direction;
  for (unsigned int d=0;d<dim;++d)
    direction[d] = d*d*d;
  DoFRenumbering::downstream(dof, direction);
  DoFRenumbering::component_wise(dof);
  for (unsigned int l=0;l<tr.n_levels();++l)
    {
      DoFRenumbering::downstream(mgdof, l, direction);
      DoFRenumbering::component_wise(mgdof, l);
    }

				   // Store sizes
  vector<unsigned int> ndofs(fe.n_blocks());
  DoFTools::count_dofs_per_block(mgdof, ndofs);
  std::vector<std::vector<unsigned int> > mg_ndofs(mgdof.get_tria().n_levels());
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
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i+1;

  std::vector<std::vector<unsigned int> > cached_sizes;
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

  for (unsigned int l=0; l<tr.n_levels();++l)
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
  for (unsigned int i=0;i<v[2].size();++i)
    {
      v[2](i) = i+1;
      wb[2](i) = i+1;
    }
  for (unsigned int i=0;i<ws[2].size();++i)
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
  std::ofstream logfile("transfer_compare_01/output");
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

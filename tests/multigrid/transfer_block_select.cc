//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 - 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------


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
#include <multigrid/mg_transfer_block.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_level_object.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void check_select(const FiniteElement<dim>& fe, unsigned int selected)
{
  deallog << fe.get_name()
	  << " select " << selected << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  
  MGDoFHandler<dim> mgdof(tr);
  DoFHandler<dim>& dof=mgdof;
  mgdof.distribute_dofs(fe);
  DoFRenumbering::component_wise(static_cast<DoFHandler<dim>&>(mgdof));
  vector<unsigned int> ndofs(fe.n_blocks());
  DoFTools::count_dofs_per_block(mgdof, ndofs);
  
  for (unsigned int l=0;l<tr.n_levels();++l)
    DoFRenumbering::component_wise(mgdof, l);
  std::vector<std::vector<unsigned int> > mg_ndofs(mgdof.get_tria().n_levels());
  MGTools::count_dofs_per_block(mgdof, mg_ndofs);

  deallog << "Global  dofs:";
  for (unsigned int i=0;i<ndofs.size();++i)
    deallog << ' ' << ndofs[i];
  deallog << std::endl;
  for (unsigned int l=0;l<mg_ndofs.size();++l)
    {
      deallog << "Level " << l << " dofs:";
      for (unsigned int i=0;i<mg_ndofs[l].size();++i)
	deallog << ' ' << mg_ndofs[l][i];
      deallog << std::endl;
    }  
  
  MGTransferBlockSelect<double> transfer;
  transfer.build_matrices(dof, mgdof, selected);

				   // First, prolongate the constant
				   // function from the coarsest mesh
				   // to the finer ones. Since this is
				   // the embedding, we obtain the
				   // constant one and the l2-norm is
				   // the number of degrees of freedom.
  Vector<double> u2(mg_ndofs[2][selected]);
  Vector<double> u1(mg_ndofs[1][selected]);
  Vector<double> u0(mg_ndofs[0][selected]);

  u0 = 1;
  transfer.prolongate(1,u1,u0);
  transfer.prolongate(2,u2,u1);
  deallog << "u0\t" << (int) (u0*u0+.5) << std::endl
	  << "u1\t" << (int) (u1*u1+.5) << std::endl
	  << "u2\t" << (int) (u2*u2+.5) << std::endl;
				   // Now restrict the same vectors.
  u1 = 0.;
  u0 = 0.;
  transfer.restrict_and_add(2,u1,u2);
  transfer.restrict_and_add(1,u0,u1);
  deallog << "u1\t" << (int) (u1*u1+.5) << std::endl
	  << "u0\t" << (int) (u0*u0+.5) << std::endl;

				   // Fill a global vector by counting
				   // from one up
  BlockVector<double> u;
  u.reinit (ndofs);
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i+1;

				   // See what part gets copied to mg
  MGLevelObject<Vector<double> > v;
  v.resize(0,tr.n_levels()-1);
  MGTools::reinit_vector_by_blocks(mgdof, v, selected, mg_ndofs);
  
  transfer.copy_to_mg(mgdof, v, u);
  for (unsigned int i=0; i<v[2].size();++i)
    deallog << ' ' << (int) v[2](i);
  deallog << std::endl;

				   // Now do the opposite: fill a
				   // multigrid vector counting the
				   // dofs and see where the numbers go
  u = 0.;
  for (unsigned int i=0;i<v[2].size();++i)
    v[2](i) = i+1;
  transfer.copy_from_mg_add(mgdof, u, v);
  for (unsigned int i=0; i<u.size();++i)
    deallog << ' ' << (int) u(i);
  deallog << std::endl;
  
}


int main()
{
  std::ofstream logfile("transfer_block_select/output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  FE_DGQ<2> q0(0);
  FE_DGQ<2> q1(1);
  FE_RaviartThomasNodal<2> rt0(0);
  FE_RaviartThomasNodal<2> rt1(1);
  
  FESystem<2> fe0(rt1, 1, q1, 1);
  FESystem<2> fe1(rt0, 2, q0, 2);
  
  check_select(fe0, 0);
  check_select(fe0, 1);
  
  check_select(fe1, 0);
  check_select(fe1, 1);
  check_select(fe1, 2);
  check_select(fe1, 3);
}

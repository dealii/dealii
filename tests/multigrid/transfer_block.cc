//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 - 2006 by the deal.II authors
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
#include <fe/fe_system.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_level_object.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void check_block(const FiniteElement<dim>& fe,
		 const vector<bool>& /*selected*/,
		 const vector<bool>& /*mg_selected*/,
		 const vector<double>& factors,
		 std::vector<unsigned int> target_component,
		 std::vector<unsigned int> mg_target_component)
{
  deallog << fe.get_name() << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  
  MGDoFHandler<dim> mgdof(tr);
  DoFHandler<dim>& dof=mgdof;
  mgdof.distribute_dofs(fe);
  DoFRenumbering::component_wise(mgdof, target_component);
  vector<unsigned int> ndofs(fe.n_components());
  DoFTools::count_dofs_per_component(mgdof, ndofs, true, target_component);
  
  for (unsigned int l=0;l<tr.n_levels();++l)
    DoFRenumbering::component_wise(mgdof, l, mg_target_component);
  
  std::vector<std::vector<unsigned int> > mg_ndofs(mgdof.get_tria().n_levels());
  MGTools::count_dofs_per_component(mgdof, mg_ndofs, true, mg_target_component);

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
  
  PrimitiveVectorMemory<Vector<double> > mem;
  MGTransferBlock<double> transfer;
  transfer.build_matrices(dof, mgdof, target_component, mg_target_component);
  if (factors.size()>0)
    transfer.initialize(factors, mem);

  BlockVector<double> u2(mg_ndofs[2]);
  BlockVector<double> u1(mg_ndofs[1]);
  BlockVector<double> u0(mg_ndofs[0]);

  u0 = 1;
  transfer.prolongate(1,u1,u0);
  transfer.prolongate(2,u2,u1);
  deallog << "u0";
  for (unsigned int b=0;b<u0.n_blocks();++b)
    deallog << '\t' << (int) (u0.block(b)*u0.block(b)+.5);
  deallog << std::endl << "u1";
  for (unsigned int b=0;b<u1.n_blocks();++b)
    deallog << '\t' << (int) (u1.block(b)*u1.block(b)+.5);
  deallog << std::endl << "u2";
  for (unsigned int b=0;b<u2.n_blocks();++b)
    deallog << '\t' << (int) (u2.block(b)*u2.block(b)+.5);
  deallog << std::endl;
  
  u1 = 0.;
  u0 = 0.;
  transfer.restrict_and_add(2,u1,u2);
  transfer.restrict_and_add(1,u0,u1);

  deallog << "u1";
  for (unsigned int b=0;b<u1.n_blocks();++b)
    deallog << '\t' << (int) (u1.block(b)*u1.block(b)+.5);
  deallog << std::endl << "u0";
  for (unsigned int b=0;b<u0.n_blocks();++b)
    deallog << '\t' << (int) (u0.block(b)*u0.block(b)+.5);
  deallog << std::endl;
  				   // Check copy to mg and back
  BlockVector<double> u;
  u.reinit (ndofs);
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i;
  
  MGLevelObject<BlockVector<double> > v;
  v.resize(2,2);
  v[2].reinit(mg_ndofs[2]);

  transfer.copy_to_mg(mgdof, v, u);
  for (unsigned int i=0; i<v[2].size();++i)
    deallog << ' ' << (int) v[2](i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("transfer_block/output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  std::vector<unsigned int> v1(4);
  std::vector<unsigned int> v2(4);
  std::vector<unsigned int> v3(4);
  for (unsigned int i=0;i<4;++i)
    {
      v1[i] = i;
      v2[i] = 0;
      v3[i] = i;
    }
  v3[0] = 0;
  v3[1] = 1;
  v3[2] = 1;
  v3[3] = 2;

  std::vector<double> factors;
  
  FESystem<2> fe1(FE_DGQ<2>(1), 4);
  
  vector<bool> s1(4, true);
  
  deallog << "s1 s1 v1 v1" << std::endl;
  check_block(fe1, s1, s1, factors, v1, v1);
//   deallog << "s1 s1 v1 v2" << std::endl;
//   check_block(fe1, s1, s1, factors, v1, v2);
//   deallog << "s1 s1 v1 v3" << std::endl;
//   check_block(fe1, s1, s1, factors, v1, v3);
  
  factors.resize(4, 1.);
  factors[1] = 2.;
  deallog << "s1 s1 v1 v1" << std::endl;
  check_block(fe1, s1, s1, factors, v1, v1);  
//   deallog << "s1 s1 v1 v3" << std::endl;
//   check_block(fe1, s1, s1, factors, v1, v3);  
}

//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 - 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------


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
void check_simple(const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  
  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);

  MGTransferPrebuilt<Vector<double> > transfer;
  transfer.build_matrices(mgdof);

  Vector<double> u2(mgdof.n_dofs(2));
  Vector<double> u1(mgdof.n_dofs(1));
  Vector<double> u0(mgdof.n_dofs(0));

  u0 = 1;
  transfer.prolongate(1,u1,u0);
  transfer.prolongate(2,u2,u1);
  deallog << "u0 " << u0.l2_norm()
	  << "\tu1 " << u1.l2_norm()
	  << "\tu2 " << u2.l2_norm();
  u1 = 0.;
  u0 = 0.;
  transfer.restrict_and_add(2,u1,u2);
  transfer.restrict_and_add(1,u0,u1);
  deallog << "\tu1 " << u1.l2_norm()
	  << "\tu0 " << u0.l2_norm()
	  << std::endl;
}

template <int dim>
void check_select(const FiniteElement<dim>& fe,
		  unsigned int selected,
		  unsigned int mg_selected,
		  std::vector<unsigned int> target_component,
		  std::vector<unsigned int> mg_target_component)
{
  deallog << fe.get_name()
	  << " select " << selected
	  << " (global) and " << mg_selected
	  << " (mg)" << std::endl;
  
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
  std::vector<std::vector<unsigned int> > mg_ndofs(mgdof.get_tria().n_levels());  MGTools::count_dofs_per_component(mgdof, mg_ndofs, true, mg_target_component);

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
  
  
  MGTransferSelect<double> transfer;
  transfer.build_matrices(dof, mgdof, selected, mg_selected,
			  target_component, mg_target_component);

  Vector<double> u2(mg_ndofs[2][mg_selected]);
  Vector<double> u1(mg_ndofs[1][mg_selected]);
  Vector<double> u0(mg_ndofs[0][mg_selected]);

  u0 = 1;
  transfer.prolongate(1,u1,u0);
  transfer.prolongate(2,u2,u1);
  deallog << "u0 " << u0*u0
	  << "\tu1 " << u1*u1
	  << "\tu2 " << u2*u2;
  u1 = 0.;
  u0 = 0.;
  transfer.restrict_and_add(2,u1,u2);
  transfer.restrict_and_add(1,u0,u1);
  deallog << "\tu1 " << u1*u1
	  << "\tu0 " << u0*u0
	  << std::endl;

  				   // Check copy to mg and back
  BlockVector<double> u;
  u.reinit (ndofs);
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i;
  
  MGLevelObject<Vector<double> > v;
  v.resize(2,2);
  v[2].reinit(mg_ndofs[2][mg_selected]);

  transfer.copy_to_mg(mgdof, v, u);
  for (unsigned int i=0; i<v[2].size();++i)
    deallog << ' ' << (int) v[2](i);
  deallog << std::endl;
}

/*
 * The following code must be adapted to the new structures.

template <int dim>
void check_block(const FiniteElement<dim>& fe,
		 const vector<bool>& selected,
		 const vector<bool>& mg_selected,
		 std::vector<unsigned int> target_component,
		 std::vector<unsigned int> mg_target_component)
{
  deallog << fe.get_name() << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  
  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  DoFRenumbering::component_wise(mgdof, target_component);
  for (unsigned int l=0;l<tr.n_levels();++l)
    DoFRenumbering::component_wise(mgdof, l, mg_target_component);
  
  std::vector<std::vector<unsigned int> > mg_ndofs(mgdof.get_tria().n_levels());
  MGTools::count_dofs_per_component(mgdof, mg_ndofs, mg_target_component);
  for (unsigned int l=0;l<mg_ndofs.size();++l)
    {
      deallog << "Level " << l << " dofs:";
      for (unsigned int i=0;i<mg_ndofs[l].size();++i)
	deallog << ' ' << mg_ndofs[l][i];
      deallog << std::endl;
    }
  
  
  MGTransferBlock<double> transfer;
  transfer.build_matrices(mgdof, selected, mg_selected,
			  target_component, mg_target_component);

  BlockVector<double> u2(mg_ndofs[2]);
  BlockVector<double> u1(mg_ndofs[1]);
  BlockVector<double> u0(mg_ndofs[0]);

  u0 = 1;
  transfer.prolongate(1,u1,u0);
  transfer.prolongate(2,u2,u1);
  deallog << "u0 " << u0*u0
	  << "\tu1 " << u1*u1
	  << "\tu2 " << u2*u2;
  u1 = 0.;
  u0 = 0.;
  transfer.restrict_and_add(2,u1,u2);
  transfer.restrict_and_add(1,u0,u1);
  deallog << "\tu1 " << u1*u1
	  << "\tu0 " << u0*u0
	  << std::endl;

  				   // Check copy to mg and back
  std::vector<unsigned int> ndofs (target_component.size());
  DoFTools::count_dofs_per_component(mgdof, ndofs, target_component);
  BlockVector<double> u;
  u.reinit (ndofs);
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i;
  
  MGLevelObject<Vector<double> > v;
  v.resize(2,2);
  v[2].reinit(mg_ndofs[2][mg_selected]);

  transfer.copy_to_mg(mgdof, v, u);
  for (unsigned int i=0; i<v[2].size();++i)
    deallog << ' ' << (int) v[2](i);
  deallog << std::endl;
}
*/

int main()
{
  std::ofstream logfile("transfer.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
//  check_simple (FE_DGP<2>(0));
//  check_simple (FE_DGP<2>(1));
  check_simple (FE_DGQ<2>(1));
//  check_simple (FE_DGQ<2>(2));
//  check_simple (FE_Q<2>(1));
//  check_simple (FE_Q<2>(2));

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

  FESystem<2> fe1(FE_DGQ<2>(1), 4);
  
  check_select (fe1, 0, 0, v1, v1);
  check_select (fe1, 1, 1, v1, v1);
  check_select (fe1, 1, 2, v1, v1);
  check_select (fe1, 2, 1, v1, v1);
  check_select (fe1, 0, 0, v2, v2);
  check_select (fe1, 0, 0, v3, v3);
  check_select (fe1, 1, 1, v3, v3);
  check_select (fe1, 0, 0, v1, v3);
  check_select (fe1, 1, 1, v1, v3);
  check_select (fe1, 2, 1, v1, v3);
  check_select (fe1, 0, 1, v1, v3);
  check_select (fe1, 0, 0, v1, v2);

//  vector<bool> s1(4, true);

//  check_block(fe1, s1, s1, v1, v1);
}

//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 - 2007, 2009 by the deal.II authors
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
#include <multigrid/mg_transfer_component.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_level_object.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

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
  DoFRenumbering::component_wise(static_cast<DoFHandler<dim>&>(mgdof),
                                 target_component);
  vector<unsigned int> ndofs(*std::max_element (target_component.begin(),
						target_component.end()) + 1);
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
  
  
  MGTransferSelect<double> transfer;
  transfer.build_matrices(dof, mgdof, selected, mg_selected,
			  target_component, mg_target_component);

  Vector<double> u2(mg_ndofs[2][mg_selected]);
  Vector<double> u1(mg_ndofs[1][mg_selected]);
  Vector<double> u0(mg_ndofs[0][mg_selected]);

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


int main()
{
  std::ofstream logfile("transfer_select/output");
  deallog << std::setprecision(3);
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
}

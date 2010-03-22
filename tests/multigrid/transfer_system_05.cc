//----------------------------------------------------------------------------
//    transfer.cc,v 1.13 2005/12/30 16:07:03 guido Exp
//    Version:
//
//    Copyright (C) 2000 - 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// like _04, but checks the copy_to_mg and copy_from_mg of MGTransferSelect

#include "../tests.h"
#include <base/logstream.h>
#include <base/mg_level_object.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_transfer_component.h>
#include <multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;


template <int dim>
void check (const FiniteElement<dim>& fe, const unsigned int selected_block)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MGDoFHandler<dim> mg_dof_handler(tr);
  mg_dof_handler.distribute_dofs(fe);

  std::vector<unsigned int> block_component(5,0);
  block_component[2]=1;
  block_component[3]=1;
  block_component[4]=2;

  DoFRenumbering::component_wise (mg_dof_handler, block_component);
  for (unsigned int level=0; level<tr.n_levels(); ++level)
    DoFRenumbering::component_wise (mg_dof_handler, level, block_component);


  MGTransferSelect<double> transfer;

  transfer.build_matrices(mg_dof_handler, mg_dof_handler,
			  selected_block, selected_block,
			  block_component, block_component);

  std::vector<unsigned int> dofs_per_block(3);
  DoFTools::count_dofs_per_block(mg_dof_handler, dofs_per_block, block_component);
  std::vector<std::vector<unsigned int> > mg_dofs_per_block(tr.n_levels(), std::vector<unsigned int>(3));
  MGTools::count_dofs_per_block(mg_dof_handler, mg_dofs_per_block, block_component);

  deallog << "Global  dofs:";
  for (unsigned int i=0;i<dofs_per_block.size();++i)
    deallog << ' ' << dofs_per_block[i];
  deallog << std::endl;
  for (unsigned int l=0;l<mg_dofs_per_block.size();++l)
    {
      deallog << "Level " << l << " dofs:";
      for (unsigned int i=0;i<mg_dofs_per_block[l].size();++i)
	deallog << ' ' << mg_dofs_per_block[l][i];
      deallog << std::endl;
    }

  BlockVector<double> u;
  u.reinit (dofs_per_block);
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i;
  
  MGLevelObject<Vector<double> > v(0,tr.n_levels()-1);
  for(unsigned int l=0; l<tr.n_levels()-1; ++l)
    v[l].reinit(mg_dofs_per_block[l][2]);


  transfer.copy_to_mg(mg_dof_handler, v, u);
  for(unsigned int l=0; l<tr.n_levels(); ++l)
  {
    deallog << "Level " << l << std::endl;
    for (unsigned int i=0; i<v[l].size();++i)
      deallog << ' ' << (int) v[l](i);
    deallog << std::endl;
  }

  for(unsigned int l=0; l<tr.n_levels(); ++l)
    for (unsigned int i=0;i<v[l].size();++i)
      v[l](i) = i;
  
  u=0;
  transfer.copy_from_mg(mg_dof_handler, u, v);
  deallog << "Global" << std::endl;
  for (unsigned int i=0; i<u.block(selected_block).size();++i)
      deallog << ' ' << (int) u.block(selected_block)(i);
  deallog << std::endl;
  for (unsigned int i=0; i<u.size();++i)
      deallog << ' ' << (int) u(i);
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("transfer_system_05/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

//TODO: do in 1d
  check (FESystem<2>(FE_Q<2>(1),5),0);
  check (FESystem<2>(FE_Q<2>(1),5),1);
  check (FESystem<2>(FE_Q<2>(1),5),2);
}

//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>

#include <deal.II/lac/block_list.h>

#include <algorithm>

using namespace dealii;


void
print_patches (const SparsityPattern& bl)
{
  for (unsigned int i=0;i<bl.n_rows();++i)
    {
      deallog << "Block " << std::setw(3) << i;
      std::vector<unsigned int> entries;
      for (SparsityPattern::row_iterator b = bl.row_begin(i);b != bl.row_end(i);++b)
	entries.push_back(*b);
      
      std::sort(entries.begin(), entries.end());
      
      for (unsigned int i=0;i<entries.size();++i)
	deallog << ' ' << std::setw(4) << entries[i];
      deallog << std::endl;
    }
}


template <int dim>
void
test_block_list(const Triangulation<dim>& tr, const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  
  MGDoFHandler<dim> dof;
  dof.initialize(tr, fe); 
  
  const unsigned int level = tr.n_levels()-1;

  {
    deallog.push("ttt");
    SparsityPattern bl;
    DoFTools::make_vertex_patches(bl, dof, level, true, true, true);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ttf");
    SparsityPattern bl;
    DoFTools::make_vertex_patches(bl, dof, level, true, true, false);
    bl.compress();
    print_patches(bl);
    deallog.pop(); 
    deallog << std::endl;
  }
  {
    deallog.push("tff");
    SparsityPattern bl;
    DoFTools::make_vertex_patches(bl, dof, level, true, false, false);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("ftt");
    SparsityPattern bl;
    DoFTools::make_vertex_patches(bl, dof, level, false, true, true);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("fff");
    SparsityPattern bl;
    DoFTools::make_vertex_patches(bl, dof, level, false, false, false);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
}


int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  Triangulation<2> trc, trl;
  GridGenerator::hyper_cube(trc);
  trc.refine_global(2);
  GridGenerator::hyper_L(trl);
  trl.refine_global(2);
  
  FE_DGQ<2> fe1(0);
  FE_DGQ<2> fe2(1);
  FE_RaviartThomas<2> fe3(0);
  FE_RaviartThomas<2> fe4(1);
  FE_Nedelec<2> fe5(0);
  FE_Nedelec<2> fe6(1);
  
  deallog.push("Square");
  test_block_list(trc, fe1);
  test_block_list(trc, fe2);
  test_block_list(trc, fe3);
  test_block_list(trc, fe4);
  test_block_list(trc, fe5);
  test_block_list(trc, fe6);
  deallog.pop();
  deallog.push("L");
  test_block_list(trl, fe1);
  test_block_list(trl, fe2);
  test_block_list(trl, fe3);
  test_block_list(trl, fe4);
  test_block_list(trl, fe5);
  test_block_list(trl, fe6);
  deallog.pop();
}

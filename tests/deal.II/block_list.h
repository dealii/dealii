//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2011 by the deal.II authors
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
#include <deal.II/dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>

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
void test_global_refinement(
  void (*test_block_list)(const Triangulation<dim>& tr, const FiniteElement<dim>& fe))
{
  Triangulation<dim> trc, trl;
  GridGenerator::hyper_cube(trc);
  trc.refine_global(2);
  GridGenerator::hyper_L(trl);
  trl.refine_global(1);
  
  FE_DGQ<dim> fe1(0);
  FE_DGQ<dim> fe2(1);
//  FE_Q<dim> fe3(1);
  FE_RaviartThomas<dim> fe4(0);
  FE_RaviartThomas<dim> fe5(1);
  FE_Nedelec<dim> fe6(0);
  FE_Nedelec<dim> fe7(1);
  
  deallog.push("Square");
  test_block_list(trc, fe1);
  test_block_list(trc, fe2);
//  test_block_list(trc, fe3);
  test_block_list(trc, fe4);
  test_block_list(trc, fe5);
  test_block_list(trc, fe6);
  test_block_list(trc, fe7);
  deallog.pop();
  deallog.push("L");
  test_block_list(trl, fe1);
  test_block_list(trl, fe2);
//  test_block_list(trl, fe3);
  test_block_list(trl, fe4);
  test_block_list(trl, fe5);
  test_block_list(trl, fe6);
  test_block_list(trl, fe7);
  deallog.pop();
}

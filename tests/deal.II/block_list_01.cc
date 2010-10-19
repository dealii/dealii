//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_raviart_thomas.h>

#include <lac/block_list.h>

using namespace dealii;


template <int dim>
void
test_block_list(const Triangulation<dim>& tr, const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  
  DoFHandler<dim> dof;
  dof.initialize(tr, fe);
  
  BlockList bl;
  bl.initialize(tr.n_active_cells(), dof.begin_active(), dof.end());

  for (unsigned int i=0;i<bl.size();++i)
    {
      deallog << "Block " << std::setw(3) << i;
      for (BlockList::const_iterator b = bl.begin(i);b != bl.end(i);++b)
	deallog << ' ' << std::setw(4) << (*b);
      deallog << std::endl;
    }
}


int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  Triangulation<2> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  FE_Q<2> fe1(1);
  FE_Q<2> fe2(2);
  FE_RaviartThomas<2> fe3(0);
  
  test_block_list(tr, fe1);
  test_block_list(tr, fe2);
  test_block_list(tr, fe3);
}

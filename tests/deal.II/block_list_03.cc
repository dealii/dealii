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

#include "block_list.h"

template <int dim>
void
test_block_list(const Triangulation<dim>& tr, const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  
  MGDoFHandler<dim> dof;
  dof.initialize(tr, fe); 
  
  const unsigned int level = tr.n_levels()-1;

  {
    deallog.push("t");
    SparsityPattern bl;
    DoFTools::make_single_patch(bl, dof, level, true);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
  {
    deallog.push("f");
    SparsityPattern bl;
    DoFTools::make_single_patch(bl, dof, level, false);
    bl.compress();
    print_patches(bl);
    deallog.pop();
    deallog << std::endl;
  }
}


int main()
{
  initlog(__FILE__);
  deallog.push("2D");
  test_global_refinement<2>(&test_block_list<2>);
  deallog.pop();
  deallog.push("3D");
  test_global_refinement<3>(&test_block_list<3>);
  deallog.pop();
}

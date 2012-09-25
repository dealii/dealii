//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// tests for the ComponentMask class
//
// here: test conversion from component mask to block mask and back. start with
// blocks where the _12 test starts with components, and use an element that
// does have fewer blocks than components


#include "../tests.h"
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include <fstream>
#include <iomanip>




void test ()
{
				   // the following element has more
				   // components than blocks
  FESystem<2> fe (FE_RaviartThomas<2>(1), 3);

				   // try all possible block
				   // masks, which we encode as bit
				   // strings
  for (unsigned int int_mask=0; int_mask<(1U<<fe.n_blocks()); ++int_mask)
  {
    std::vector<bool> block_mask (fe.n_blocks());
    for (unsigned int c=0; c<fe.n_blocks(); ++c)
      block_mask[c] = (int_mask & (1<<c));

				     // make sure that the round-trip works
    Assert (BlockMask(block_mask) == fe.block_mask(fe.component_mask(BlockMask(block_mask))),
	    ExcInternalError());
  }

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("component_mask_13/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}

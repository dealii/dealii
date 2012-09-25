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
// here: test conversion from component mask to block mask does indeed
// not work if we try to split a block


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
  FESystem<2> fe (FE_RaviartThomas<2>(1), 2);

				   // try all possible component
				   // masks, which we encode as bit
				   // strings
  for (unsigned int int_mask=0; int_mask<(1U<<fe.n_components()); ++int_mask)
  {
    std::vector<bool> component_mask (fe.n_components());
    for (unsigned int c=0; c<fe.n_components(); ++c)
      component_mask[c] = (int_mask & (1<<c));

				     // now try to convert to a block
				     // mask. this should not always
				     // work
    deallog << fe.block_mask (component_mask) << std::endl;
  }

  deallog << "OK" << std::endl;
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("component_mask_14/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}

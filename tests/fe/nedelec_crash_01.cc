//--------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------

// A test extracted from integrators/cochain_01 that crashed with the
// FE_Nedelec at the time

#include "../tests.h"
#include "../lib/test_grids.h"

#include <base/logstream.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>

#include <fe/fe_nedelec.h>


int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

				   // generate a version of the
				   // Nedelec element but force it to
				   // use the old-style constraints
  struct MyFE : FE_Nedelec<3> {
      MyFE () : FE_Nedelec<3>(0) {};
      virtual bool hp_constraints_are_implemented () const
	{
	  return false;
	}
  } fe;

  Triangulation<3> tr(Triangulation<3>::limit_level_difference_at_vertices);
  TestGrids::hypercube(tr, 2, true);

  DoFHandler<3> dof(tr);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();

  constraints.print (deallog.get_file_stream ());
}

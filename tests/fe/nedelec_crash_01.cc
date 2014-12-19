// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// A test extracted from integrators/cochain_01 that crashed with the
// FE_Nedelec at the time

#include "../tests.h"
#include "../test_grids.h"

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>


int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // generate a version of the
  // Nedelec element but force it to
  // use the old-style constraints
  struct MyFE : FE_Nedelec<3>
  {
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

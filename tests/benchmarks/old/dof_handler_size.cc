// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


#include <iomanip>
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;


// Fill dof handlers for different elements and see how large they get.
template <class DOF>
void check_dofs(DOF &dof)
{
  FE_Q<DOF::dimension> q1(1);
  FE_Q<DOF::dimension> q3(3);
  FESystem<DOF::dimension> sys1(q3, 1);
  FESystem<DOF::dimension> sys2(q3, 10);

  dof.distribute_dofs(q1);
  deallog << "Dofs Q1   " << std::setw(12) << dof.n_dofs()
          << " memory " << std::setw(12) << dof.memory_consumption()
          << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
          << std::endl;

  dof.distribute_dofs(q3);
  deallog << "Dofs Q3   " << std::setw(12) << dof.n_dofs()
          << " memory " << std::setw(12) << dof.memory_consumption()
          << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
          << std::endl;

  dof.distribute_dofs(sys1);
  deallog << "Dofs Sys1 " << std::setw(12) << dof.n_dofs()
          << " memory " << std::setw(12) << dof.memory_consumption()
          << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
          << std::endl;

  dof.distribute_dofs(sys2);
  deallog << "Dofs Sys2 " << std::setw(12) << dof.n_dofs()
          << " memory " << std::setw(12) << dof.memory_consumption()
          << " quotient " << (1./dof.n_dofs()*dof.memory_consumption())
          << std::endl;

  dof.clear();
}


template <int dim>
void check (bool local)
{
  deallog << "Dimension " << dim << std::endl;
  Triangulation<dim> tr(Triangulation<dim>::maximum_smoothing);
  GridGenerator::hyper_cube(tr);

  if (local)
    for (unsigned int i=0; i<99/dim; ++i)
      {
        tr.last_active()->set_refine_flag();
        tr.execute_coarsening_and_refinement();
      }
  else
    tr.refine_global(18/dim);

  deallog << "Levels " << tr.n_levels() << "Cells "<< tr.n_cells()
          << std::endl
          << " active " << std::setw(12)  << tr.n_active_cells()
          << " memory " << std::setw(12)  << tr.memory_consumption()
          << " quotient " << (1./tr.n_cells()*tr.memory_consumption())
          << std::endl;

  FE_Q<dim> q1(1);
  FE_Q<dim> q3(3);
  FESystem<dim> sys1(q3, 1);
  FESystem<dim> sys2(q3, 10);
  deallog.push("DoF");
  DoFHandler<dim> dof(tr);
  check_dofs(dof);
  deallog.pop();
  deallog.push("MGDoF");
  MGDoFHandler<dim> mgdof(tr);
  check_dofs(mgdof);
  deallog.pop();
}

int main()
{
  std::ofstream out("dof_handler_size/output");
  deallog.attach(out);
  deallog.push("local");
  check<2>(true);
  check<3>(true);
  deallog.pop();
  deallog.push("global");
  check<2>(false);
  check<3>(false);
  deallog.pop();
}


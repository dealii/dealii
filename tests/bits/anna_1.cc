// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



// check some things about Nedelec elements, basically also that the
// DoFRenumbering::component_wise function also works for
// non_primitive elements, for which it did not work previously since
// there is no component to associate a non-primitive shape function
// with
//
// this program is a modified version of one by Anna Schneebeli,
// University of Basel

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_base.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <iomanip>
#include <fstream>


template <int dim>
class SystemTest
{
public:
  SystemTest ();
  void run ();

private:
  void make_grid_and_dofs ();
  void shape_to_components ();
  void check_numbering ();


  Triangulation<dim>     triangulation;
  FESystem<dim>          fe;
  DoFHandler<dim>        dof_handler;


};

template <int dim>
SystemTest<dim>::SystemTest () :
  fe (FE_Nedelec<dim>(0), 2,
      FE_Q<dim>(1), 1),
  dof_handler (triangulation)
{}


template <int dim>
void SystemTest<dim>::make_grid_and_dofs ()
{

  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (0);
  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: " << triangulation.n_cells()
          << std::endl;

  dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

}

template <int dim>
void SystemTest<dim>::shape_to_components ()
{
  // testing, if the shape function
  // with index i is of type Nedelec:
  // (i.e. the first component of the FESystem)
  // 1 for yes, 0 for no.

  for (unsigned int i = 0; i<fe.dofs_per_cell; i++)
    deallog <<"  shapefunction "<< i << " is Nedelec:  "
            << (fe.is_primitive(i) ? "false" : "true") << std::endl;
}



template <int dim>
void SystemTest<dim>::check_numbering ()
{
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  std::vector<types::global_dof_index>  local_dof_indices(fe.dofs_per_cell);

  for (; cell!=endc; ++cell)
    {
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<fe.dofs_per_cell; i++)
        deallog <<"  DoF "<< local_dof_indices[i] << " belongs to base element "
                << fe.system_to_base_index(i).first.first
                << ", instance " << fe.system_to_base_index(i).first.second
                << std::endl;
      deallog<< std::endl;
    };


  //Now: Componentwise reodering of the dofs

  deallog << "  Now we renumber the DoFs component-wise:" << std::endl;
  deallog << "  ****************************************" << std::endl;
  DoFRenumbering::component_wise (dof_handler);

  cell = dof_handler.begin_active();
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<fe.dofs_per_cell; i++)
        deallog <<"  DoF "<< local_dof_indices[i] << " belongs to base "
                << fe.system_to_base_index(i).first.first
                << ", instance " << fe.system_to_base_index(i).first.second
                << std::endl;
      deallog << std::endl;
    };
}


template <int dim>
void SystemTest<dim>::run ()
{
  make_grid_and_dofs ();
  shape_to_components ();
  check_numbering();
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SystemTest<2>().run();
  SystemTest<3>().run();
  return 0;
}

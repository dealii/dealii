// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check some things about Nedelec elements, here that
// DoFTools::component_select and DoFTools::count_dofs_per_fe_component
// works
//
// this program is a modified version of one by Anna Schneebeli,
// University of Basel

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_base.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
class SystemTest
{
public:
  SystemTest();
  void
  run();

private:
  void
  make_grid_and_dofs();
  void
  check();


  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;
};

template <int dim>
SystemTest<dim>::SystemTest()
  : fe(FE_Nedelec<dim>(0), 2, FE_Q<dim>(1), 1)
  , dof_handler(triangulation)
{}


template <int dim>
void
SystemTest<dim>::make_grid_and_dofs()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(0);
  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: " << triangulation.n_cells() << std::endl;

  dof_handler.distribute_dofs(fe);
  deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;
}


template <int dim>
void
SystemTest<dim>::check()
{
  for (unsigned int c = 0; c < fe.n_components(); ++c)
    {
      deallog << "Checking for component " << c << std::endl;
      std::vector<bool> x(fe.n_components(), false);
      x[c] = true;
      const IndexSet sel =
        DoFTools::extract_dofs(dof_handler, ComponentMask(x));

      for (unsigned int i = 0; i < sel.size(); ++i)
        if (sel.is_element(i))
          deallog << "  DoF " << i << std::endl;
    };

  const std::vector<types::global_dof_index> dofs_per_component =
    DoFTools::count_dofs_per_fe_component(dof_handler);
  deallog << "DoFs per component: ";
  for (unsigned int i = 0; i < fe.n_components(); ++i)
    deallog << dofs_per_component[i] << ' ';
  deallog << std::endl;
}


template <int dim>
void
SystemTest<dim>::run()
{
  deallog << "************* " << dim << "D *************" << std::endl;
  make_grid_and_dofs();
  check();

  // renumber degrees of freedom and
  // try again
  deallog << std::endl << "*** Renumbering ***" << std::endl;
  DoFRenumbering::component_wise(dof_handler);
  check();
}



int
main()
{
  initlog();

  SystemTest<2>().run();
  SystemTest<3>().run();
  return 0;
}

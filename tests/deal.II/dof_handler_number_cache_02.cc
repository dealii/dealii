// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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



// Check that DoFHandler::clear() clears the NumberCache (a bug that is now fixed)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>
#include <cstdlib>


template<int dim>
void test()
{
  Triangulation<dim> triangulation (Triangulation<dim>::limit_level_difference_at_vertices);

  FESystem<dim> fe (FE_Q<dim>(3),2,
                    FE_DGQ<dim>(1),1);

  DoFHandler<dim> dof_handler (triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (2);

  const unsigned int n_refinements[] = { 0, 4, 3, 2 };
  for (unsigned int i=0; i<n_refinements[dim]; ++i)
    {
      // refine one-fifth of cells randomly
      std::vector<bool> flags (triangulation.n_active_cells(), false);
      for (unsigned int k=0; k<flags.size()/5 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      // refine triangulation
      unsigned int index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell, ++index)
        if (flags[index])
          cell->set_refine_flag();
      Assert (index == triangulation.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell, ++index)
        if (!flags[index])
          cell->set_coarsen_flag();

      triangulation.execute_coarsening_and_refinement ();
      dof_handler.distribute_dofs (fe);

      const unsigned int N = dof_handler.n_dofs();
      deallog << N << std::endl;

      IndexSet all (N);
      all.add_range (0, N);

      Assert (dof_handler.n_locally_owned_dofs() == N,
              ExcInternalError());
      Assert (dof_handler.locally_owned_dofs() == all,
              ExcInternalError());
      Assert (dof_handler.n_locally_owned_dofs_per_processor() ==
              std::vector<types::global_dof_index> (1,N),
              ExcInternalError());
      Assert (dof_handler.locally_owned_dofs_per_processor() ==
              std::vector<IndexSet>(1,all),
              ExcInternalError());

      dof_handler.clear();
      deallog << "those should be zero: "
              << dof_handler.n_locally_owned_dofs()
              << " " << dof_handler.n_locally_owned_dofs_per_processor().size()
              << " " << dof_handler.n_dofs()
              << std::endl;
    }
}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check DoFHandler<1>::reserve_space for continuous elements, but
// where we use the same element on all cells


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(1);

  const hp::FECollection<dim> fe_collection(FE_Q<dim>(1));

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe_collection);

  std::vector<types::global_dof_index> local_dof_indices(
    fe_collection[0].dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      cell->get_dof_indices(local_dof_indices);

      deallog << cell << std::endl;
      for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
        deallog << local_dof_indices[i] << ' ';
      deallog << std::endl;
    }
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();

  deallog << "OK" << std::endl;
}

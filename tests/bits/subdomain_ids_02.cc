// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check DoFTools::get_subdomain_association


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << dim << 'D' << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(2);

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (; cell != endc; ++cell)
    {
      unsigned int subdomain = 0;
      for (unsigned int d = 0; d < dim; ++d)
        if (cell->center()[d] > 0)
          subdomain |= (1 << d);
      AssertThrow(subdomain < (1 << dim), ExcInternalError());

      cell->set_subdomain_id(subdomain);
    };

  // distribute some degrees of freedom and
  // output some information on them
  FESystem<dim>   fe(FE_Q<dim>(2), dim, FE_DGQ<dim>(1), 1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  deallog << dof_handler.n_dofs() << std::endl;

  std::vector<types::subdomain_id> subdomain_association(dof_handler.n_dofs());
  DoFTools::get_subdomain_association(dof_handler, subdomain_association);
  for (unsigned int subdomain = 0; subdomain < (1 << dim); ++subdomain)
    {
      // count number on dofs on
      // subdomain. this time it should add
      // up, since each dof is uniquely
      // associated
      deallog << std::count(subdomain_association.begin(),
                            subdomain_association.end(),
                            subdomain)
              << std::endl;
    }

  // make sure that all subdomain_ids are
  // really in the allowed range. if this is
  // the case, then we have also proven that
  // the numbers really add up correctly,
  // since every dof is assigned a valid
  // subdomain id
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    AssertThrow(subdomain_association[i] < (1 << dim), ExcInternalError());
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}

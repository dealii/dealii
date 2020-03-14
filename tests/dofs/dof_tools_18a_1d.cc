// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <string>

#include "../tests.h"

// check
//   DoFTools::
//   make_flux_sparsity_pattern (const DoFHandler<dim> &,
//                           SparsityPattern       &,
//                           AffineConstraints<double>      &);
// this used to fail at one point because we forgot that in 1d
// neighboring cells can be more than one level apart. with
// FE_DGQ(0) elements, we get the following mesh and DoF indices:
//
// *-*-*-*-*-------*
//  1 2 3 4    0
//
// we need to make sure that in row 0 there is a coupling to
// DoF 4. this was previously missing



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  // create sparsity pattern
  SparsityPattern sp(dof_handler.n_dofs(),
                     dof_handler.max_couplings_between_dofs() * 2);
  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       sp,
                                       AffineConstraints<double>());
  sp.compress();

  // write out 20 lines of this
  // pattern (if we write out the
  // whole pattern, the output file
  // would be in the range of 40 MB)
  for (unsigned int l = 0; l < sp.n_rows(); ++l)
    {
      for (unsigned int c = 0; c < sp.row_length(l); ++c)
        deallog << sp.column_number(l, c) << " ";
      deallog << std::endl;
    }

  // write out some other indicators
  deallog << sp.bandwidth() << std::endl
          << sp.max_entries_per_row() << std::endl
          << sp.n_nonzero_elements() << std::endl;

  unsigned int hash = 0;
  for (unsigned int l = 0; l < sp.n_rows(); ++l)
    hash +=
      l * (sp.row_length(l) + (sp.begin(l) - sp.begin()) +
           (sp.row_length(l) > 1 ? ++sp.begin(l) : sp.begin(l))->column());
  deallog << hash << std::endl;
}


void
check_this()
{
  // create a mesh where two cells at levels
  // 1 and 3 are adjacent
  const int          dim = 1;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  for (Triangulation<dim>::active_cell_iterator c = tr.begin_active(2);
       c != tr.end_active(2);
       ++c)
    c->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof_handler(tr);
  dof_handler.distribute_dofs(fe);

  check_this<DoFHandler<dim>>(dof_handler);
}


int
main()
{
  try
    {
      initlog();
      deallog << std::setprecision(2);

      check_this();
      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

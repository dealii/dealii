// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// VectorTools::interpolate_boundary_values produced an exception when
// used with hp::DoFHandler in 1d. Test that this is no longer the case


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <list>
#include <sstream>

#include "../tests.h"


template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution()
  {}
  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    return p[0];
  }
};


template <int dim>
void
test()
{
  Triangulation<dim>    triangulation;
  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));

  DoFHandler<dim> dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);
  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: " << triangulation.n_cells() << std::endl;

  dof_handler.distribute_dofs(fe);
  deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

  ExactSolution<dim>                        exact_solution;
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           exact_solution,
                                           boundary_values);
  if (dim == 1)
    VectorTools::interpolate_boundary_values(dof_handler,
                                             1,
                                             exact_solution,
                                             boundary_values);

  for (std::map<types::global_dof_index, double>::iterator i =
         boundary_values.begin();
       i != boundary_values.end();
       ++i)
    deallog << i->first << ' ' << i->second << std::endl;
}


int
main()
{
  try
    {
      initlog();
      deallog << std::setprecision(2);

      test<1>();
      test<2>();
      test<3>();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };

  return 0;
}

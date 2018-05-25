// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// Observe how the values of the shape functions change as we refine
// the grid. Then, evaluate the values with FEFaceValues, to
// make sure the values scale as in rt_bubbles_10 where we used FEValues.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_rt_bubbles.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 5


std::ofstream logfile("output");

template <int dim>
void
test(const unsigned int degree)
{
  FE_RT_Bubbles<dim> fe_rt_bubbles(degree);

  deallog << "Degree=" << degree << std::endl;

  for (double h = 1; h > 1. / 128; h /= 2)
    {
      deallog << "  h=" << h << std::endl;

      Triangulation<dim> tr;
      GridGenerator::hyper_cube(tr, 0., h);

      DoFHandler<dim> dof(tr);
      dof.distribute_dofs(fe_rt_bubbles);

      QTrapez<dim - 1> quadrature;

      FEFaceValues<dim> fe_values(fe_rt_bubbles, quadrature, update_values);
      fe_values.reinit(dof.begin_active(), 0);
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          deallog << "    Quadrature point " << q << ": ";
          for (unsigned int i = 0; i < fe_rt_bubbles.dofs_per_cell; ++i)
            {
              deallog << '[';
              for (unsigned int c = 0; c < fe_rt_bubbles.n_components(); ++c)
                deallog << fe_values.shape_value_component(i, q, c) << ' ';
              deallog << ']';
            }
          deallog << std::endl;
        }
    }
}



int
main()
{
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  for (unsigned int i = 1; i < 4; ++i)
    {
      test<2>(i);
      test<3>(i);
    }

  return 0;
}

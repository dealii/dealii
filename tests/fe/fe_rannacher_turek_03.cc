// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


// Interfaces being tested
#include <deal.II/fe/fe_rannacher_turek.h>
// Interfaces needed for testing
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

// Same as interpolation test in fe_rannacher_turek_01.cc but for different
// signatures of the interpolate function that were added later

void test_interpolation()
{
  Triangulation<2> tr;
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(2);

  FE_RannacherTurek<2> fe;
  const unsigned int n_dofs = fe.dofs_per_cell;

  DoFHandler<2> dofh;
  dofh.initialize(tr, fe);

  Vector<double> input_vector(dofh.n_dofs());
  for (unsigned int i = 0; i < input_vector.size(); ++i)
    {
      input_vector[i] = double(i);
    }

  Quadrature<2> quadrature(fe.get_generalized_support_points());
  FEValues<2> fev(fe, quadrature, update_values | update_JxW_values);

  typedef DoFHandler<2>::cell_iterator cell_it;
  cell_it cell = dofh.begin_active();
  for (; cell != dofh.end(); ++cell)
    {
      fev.reinit(cell);

      std::vector<double> values(quadrature.size());
      fev.get_function_values(input_vector, values);

      // original dofs
      Vector<double> local_dofs(n_dofs);
      cell->get_dof_values(input_vector, local_dofs);

      const unsigned int n_components = 5;
      for (unsigned int offset = 0; offset < n_components; ++offset)
        {
          // vector<Vector> interpolate function
          std::vector<double> interpolated_local_dofs2(n_dofs);
          std::vector<Vector<double> > offset_values(quadrature.size(),
                                                     Vector<double>(n_components));
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              offset_values[q][offset] = values[q];
            }
          fe.interpolate(interpolated_local_dofs2,
                         offset_values,
                         offset);

          // VectorSlice interpolate function
          std::vector<double> interpolated_local_dofs3(n_dofs);
          std::vector<std::vector<double> > vv_values(
            n_components,
            values);
          VectorSlice<const std::vector<std::vector<double> > > sliced_values(
            vv_values, offset, 1);
          fe.interpolate(interpolated_local_dofs3,
                         sliced_values);

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              deallog << "vector<Vector>@" << offset << ": "
                      << local_dofs[j] - interpolated_local_dofs2[j] << " ";
              deallog << "VectorSlice@" << offset << ": "
                      << local_dofs[j] - interpolated_local_dofs3[j] << " ";
            }
        }
      deallog << std::endl;
    }
}

int main()
{
  initlog();

  test_interpolation();

  return 0;
}

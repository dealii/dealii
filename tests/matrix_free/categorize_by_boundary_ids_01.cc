// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// tests MatrixFreeTools::categorize_by_boundary_ids()

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int n_refinements)
{
  const unsigned int fe_degree = 1;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);

  tria.refine_global(n_refinements);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<1>                                        quadrature(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags                = update_values;
  additional_data.mapping_update_flags_inner_faces    = update_values;
  additional_data.mapping_update_flags_boundary_faces = update_values;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::none;

  AffineConstraints<double> dummy;
  dummy.close();

  const auto process = [&]() {
    bool result = true;

    MatrixFree<dim, double> data;
    data.reinit(dof_handler, dummy, quadrature, additional_data);

    using VectorType = Vector<double>;

    VectorType vec;
    data.initialize_dof_vector(vec);

    data.template cell_loop<VectorType, VectorType>(
      [&](const auto &, auto &, const auto &, const auto cell_range) {
        for (unsigned int cell = cell_range.first; cell < cell_range.second;
             ++cell)
          {
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 ++face)
              {
                bool temp = true;

                for (unsigned int v = 1; v < data.n_components_filled(cell);
                     ++v)
                  temp &= (data.get_faces_by_cells_boundary_id(cell, face)[0] ==
                           data.get_faces_by_cells_boundary_id(cell, face)[v]);

                result &= temp;
              }
          }
      },
      vec,
      vec);

    return result;
  };

  MatrixFreeTools::categorize_by_boundary_ids(tria, additional_data);

  AssertDimension(process(), true); // categorized

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  {
    deallog.push("2D");
    test<2>(3);
    deallog.pop();
  }

  {
    deallog.push("3D");
    test<3>(3);
    deallog.pop();
  }
}

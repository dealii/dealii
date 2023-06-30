/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2022 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// Verify convergence rates for discontinuous simplex elements.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_project.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  FE_SimplexDGP<dim> fe(degree);
  deallog << "FE = " << fe.get_name() << std::endl;
  QGaussSimplex<dim> quadrature(degree + 1);

  double previous_error = 1.0;

  for (unsigned int r = 0; r < 4; ++r)
    {
      Triangulation<dim> tria_hex, tria;
      GridGenerator::hyper_cube(tria_hex);
      tria_hex.refine_global(r);
      GridGenerator::convert_hypercube_to_simplex_mesh(tria_hex, tria);

      ReferenceCell   reference_cell = tria.begin_active()->reference_cell();
      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      Vector<double>                 cell_errors(tria.n_active_cells());
      Vector<double>                 solution(dof_handler.n_dofs());
      Functions::CosineFunction<dim> function;
      AffineConstraints<double>      dummy;
      const auto &                   mapping =
        reference_cell.template get_default_linear_mapping<dim>();
      dummy.close();
      VectorTools::project(
        mapping, dof_handler, dummy, quadrature, function, solution);

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        function,
                                        cell_errors,
                                        Quadrature<dim>(
                                          fe.get_unit_support_points()),
                                        VectorTools::Linfty_norm);

      const double max_error =
        *std::max_element(cell_errors.begin(), cell_errors.end());
      deallog << "max error = " << max_error << std::endl;
      if (max_error != 0.0)
        deallog << "ratio = " << previous_error / max_error << std::endl;
      previous_error = max_error;

#if 0
      if (dim == 2)
        {
          DataOut<dim> data_out;
          data_out.attach_dof_handler(dof_handler);
          data_out.add_data_vector(solution, "u");
          data_out.build_patches(2);

          std::ofstream output("out-" + std::to_string(degree) + "-" +
                               std::to_string(r) + ".vtu");
          data_out.write_vtu(output);
        }
#endif
    }
}

int
main()
{
  initlog();

  test<2>(1);
  test<2>(2);
  test<2>(3);

  test<3>(1);
  test<3>(2);
  test<3>(3);
}

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 by the deal.II authors
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

// Test that FEFieldConvolutionFunction is actually able to interpolate
// between two non-matching grids, both when codimension is zero, and
// when codimension is one

#include <deal.II/base/function_lib.h>
#include <deal.II/base/parsed_convergence_table.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/non_matching/coupling.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/fe_field_convolution_function.templates.h>
#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"

template <int dim>
void
get_immersed_mesh(Triangulation<dim> &tria)
{
  GridGenerator::hyper_ball(tria, Point<dim>(), 2);
}

template <int dim>
void get_immersed_mesh(Triangulation<dim - 1, dim> &tria)
{
  GridGenerator::hyper_sphere(tria, Point<dim>(), 2);
}

using namespace dealii;
template <int dim, int spacedim>
void
test()
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  Triangulation<dim, spacedim>      tria;
  Triangulation<spacedim, spacedim> space_tria;

  get_immersed_mesh(tria);
  GridGenerator::hyper_cube(space_tria, -3, 3);

  ParsedConvergenceTable table({"u"},
                               {{VectorTools::L2_norm, VectorTools::H1_norm}});

  space_tria.refine_global(2);

  GridTools::Cache<dim, spacedim>      cache(tria);
  GridTools::Cache<spacedim, spacedim> space_cache(space_tria);

  FE_Q<dim, spacedim>      fe(1);
  FE_Q<spacedim, spacedim> space_fe(1);

  DoFHandler<dim, spacedim>      dh(tria);
  DoFHandler<spacedim, spacedim> space_dh(space_tria);

  deallog << "FE      : " << fe.get_name() << std::endl
          << "Space FE: " << space_fe.get_name() << std::endl;

  Functions::CosineFunction<spacedim> function;

  Functions::CutOffFunctionC1<spacedim> kernel(
    1., Point<spacedim>(), 1, -1, true);

  auto immersed_dofs = [&]() { return dh.n_dofs(); };
  table.add_extra_column("immersed_dofs", immersed_dofs, false);

  for (unsigned int cycle = 0; cycle < 5 - dim; ++cycle)
    {
      tria.refine_global(1);
      space_tria.refine_global(1);

      dh.distribute_dofs(fe);
      space_dh.distribute_dofs(space_fe);

      Vector<double> immersed_vector(dh.n_dofs());
      Vector<double> space_vector(space_dh.n_dofs());

      VectorTools::interpolate(space_dh, function, space_vector);

      double radius =
        2 * std::max(GridTools::minimal_cell_diameter(tria),
                     GridTools::minimal_cell_diameter(space_tria));
      kernel.set_radius(radius);

      Functions::FEFieldConvolutionFunction<spacedim> space_function(
        space_dh, space_cache, space_vector, kernel, QGauss<spacedim>(2));

      VectorTools::interpolate(dh, space_function, immersed_vector);
      table.error_from_exact(dh, immersed_vector, function);
    }
  table.output_table(deallog.get_file_stream());
}

int
main()
{
  initlog();
  test<1, 2>();
  test<2, 2>();
}
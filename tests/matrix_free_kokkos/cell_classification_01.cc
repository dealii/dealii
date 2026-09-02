// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test the geometry classification of cells (Cartesian, affine, general) in
// Portable::MatrixFree and the consistency of the associated compressed
// mapping-data storage.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/portable_matrix_free.h>

#include "../tests.h"

#include "Kokkos_Core.hpp"


template <int dim, int fe_degree>
void
test(const Triangulation<dim> &tria, const std::string &mesh_name)
{
  const MappingQ1<dim> mapping;
  const FE_Q<dim>      fe(fe_degree);
  DoFHandler<dim>      dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  Portable::MatrixFree<dim, double>                          mf_data;
  typename Portable::MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    update_values | update_gradients | update_JxW_values;
  mf_data.reinit(mapping,
                 dof_handler,
                 constraints,
                 QGauss<1>(fe_degree + 1),
                 additional_data);

  const QGauss<dim>  quadrature(fe_degree + 1);
  const unsigned int n_q_points_per_cell = quadrature.size();
  FEValues<dim>      fe_values(mapping, fe, quadrature, update_JxW_values);

  unsigned int n_cartesian = 0;
  unsigned int n_affine    = 0;
  unsigned int n_general   = 0;

  const auto        &graph    = mf_data.get_colored_graph();
  const unsigned int n_colors = graph.size();
  for (unsigned int color = 0; color < n_colors; ++color)
    {
      const typename Portable::MatrixFree<dim, double>::PrecomputedData
                 gpu_data  = mf_data.get_data(color);
      const auto host_data = Portable::copy_mf_data_to_host<dim, double>(
        gpu_data, additional_data.mapping_update_flags);

      unsigned int expected_offset = 0;
      for (unsigned int cell = 0; cell < host_data.n_cells; ++cell)
        {
          // The data index offsets must match those generated from the process
          // below.
          AssertThrow(host_data.data_index_offsets(cell) == expected_offset,
                      ExcInternalError());

          const auto cell_type = host_data.cell_type(cell);
          if (cell_type == dealii::internal::MatrixFreeFunctions::cartesian)
            {
              ++n_cartesian;
              expected_offset += 1;
            }
          else if (cell_type == dealii::internal::MatrixFreeFunctions::affine)
            {
              ++n_affine;
              expected_offset += 1;
            }
          else if (cell_type == dealii::internal::MatrixFreeFunctions::general)
            {
              ++n_general;
              expected_offset += n_q_points_per_cell;
            }
          else
            DEAL_II_ASSERT_UNREACHABLE();

          // The JxW values reconstructed from the compressed storage must
          // match the ones computed by FEValues on the same cell.
          fe_values.reinit(graph[color][cell]);
          for (unsigned int q = 0; q < n_q_points_per_cell; ++q)
            AssertThrow(std::abs(host_data.JxW_value(cell, q) -
                                 fe_values.JxW(q)) <
                          1e-12 * std::abs(fe_values.JxW(q)) + 1e-30,
                        ExcInternalError());
        }

      // The compressed storage must have exactly the size determined by the
      // cell types.
      AssertThrow(host_data.JxW.extent(0) == expected_offset,
                  ExcInternalError());
      AssertThrow(host_data.inv_jacobian.extent(0) == expected_offset,
                  ExcInternalError());
    }

  // All cells must be classified into one of the three types.
  AssertThrow(n_cartesian + n_affine + n_general == tria.n_active_cells(),
              ExcInternalError());

  // The cells should be classified in the predicted pattern.
  deallog << mesh_name << ": cartesian=" << n_cartesian
          << " affine=" << n_affine << " general=" << n_general << std::endl;
}


template <int dim, int fe_degree>
void
run()
{
  deallog << "dim=" << dim << " degree=" << fe_degree << std::endl;

  // All cells are Cartesian.
  {
    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);
    tria.refine_global(4 - dim);
    test<dim, fe_degree>(tria, "hyper_cube");
  }

  // All cells are affine.
  {
    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);
    tria.refine_global(4 - dim);
    GridTools::transform(
      [](const Point<dim> &p) {
        Point<dim> q = p;
        q[0] += 0.5 * p[dim - 1];
        return q;
      },
      tria);
    test<dim, fe_degree>(tria, "sheared");
  }

  // All cells are general.
  {
    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);
    tria.refine_global(4 - dim);
    GridTools::transform(
      [](const Point<dim> &p) {
        Point<dim> q = p;
        q[0] += 0.1 * p[0] * p[dim - 1];
        return q;
      },
      tria);
    test<dim, fe_degree>(tria, "warped");
  }

  // Ball: Cartesian center cell, curved general cells at the boundary
  {
    Triangulation<dim> tria;
    GridGenerator::hyper_ball(tria);
    test<dim, fe_degree>(tria, "hyper_ball");
  }
}


int
main()
{
  initlog();
  Kokkos::initialize();

  run<2, 1>();
  run<2, 2>();
  run<3, 1>();
  run<3, 2>();

  deallog << "OK" << std::endl;
  Kokkos::finalize();
  return 0;
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Portable::MatrixFree::get_cell_iterator() with coloring enabled

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/portable_matrix_free.templates.h>

#include <Kokkos_Core.hpp>

#include "../tests.h"


// Device functor to compute cell centers on GPU
template <int dim>
struct CellCenterKernel
{
  // Store computed centers: centers_out[global_cell_id * dim + d]
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> centers_out;

  static constexpr unsigned int n_q_points = (1u << dim);

  KOKKOS_FUNCTION void
  operator()(const typename Portable::MatrixFree<dim, double>::Data *data,
             const Portable::DeviceVector<double> &,
             Portable::DeviceVector<double> &) const
  {
    const unsigned int cell_id = data->cell_index;
    const unsigned int global_cell_id =
      data->precomputed_data[0].row_start /
        data->precomputed_data[0].padding_length +
      cell_id;

    Point<dim, double> cell_center;
    for (unsigned int d = 0; d < dim; ++d)
      cell_center[d] = 0.0;

    // Average all quadrature points to get cell center
    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const Point<dim, double> &q_point =
          data->get_quadrature_point(cell_id, q);
        for (unsigned int d = 0; d < dim; ++d)
          cell_center[d] += q_point[d];
      }

    for (unsigned int d = 0; d < dim; ++d)
      cell_center[d] /= n_q_points;

    // Store in output array
    const unsigned int base_idx = global_cell_id * dim;
    for (unsigned int d = 0; d < dim; ++d)
      centers_out[base_idx + d] = cell_center[d];
  }
};

template <int dim>
void
test(const unsigned int level = numbers::invalid_unsigned_int)
{
  Triangulation<dim> tria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  AffineConstraints<double> constraints;
  constraints.close();

  const QGauss<1> quad(2);

  Portable::MatrixFree<dim, double> mf_data;

  typename Portable::MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    update_gradients | update_JxW_values | update_quadrature_points;
  additional_data.mg_level = level;

  additional_data.use_coloring = true;

  mf_data.reinit(
    MappingQ1<dim>(), dof_handler, constraints, quad, additional_data);

  const unsigned int n_colors = level == numbers::invalid_unsigned_int ?
                                  mf_data.get_colored_graph().size() :
                                  mf_data.get_colored_level_graph().size();

  // Calculate total number of cells
  unsigned int n_total_cells = 0;
  for (unsigned int color = 0; color < n_colors; ++color)
    n_total_cells += mf_data.n_cells_per_color(color);

  // Allocate device array to store GPU-computed centers
  // Layout: centers[global_cell_id * dim + d]
  using double_view =
    Kokkos::View<double *, MemorySpace::Default::kokkos_space>;
  double_view centers_gpu("centers_gpu", n_total_cells * dim);

  // Run kernel to compute centers on GPU
  using device_vector_type =
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>;
  device_vector_type src, dst;
  mf_data.initialize_dof_vector(src);
  mf_data.initialize_dof_vector(dst);

  CellCenterKernel<dim> kernel{centers_gpu};
  mf_data.cell_loop(kernel, src, dst);

  // Copy GPU results to host
  auto centers_gpu_host = Kokkos::create_mirror_view(centers_gpu);
  Kokkos::deep_copy(centers_gpu_host, centers_gpu);

  // Compare with host-computed centers using get_cell_iterator()
  unsigned int n_matches      = 0;
  unsigned int n_mismatches   = 0;
  unsigned int global_cell_id = 0;

  for (unsigned int color = 0; color < n_colors; ++color)
    {
      const unsigned int n_cells = mf_data.n_cells_per_color(color);

      for (unsigned int cell_id = 0; cell_id < n_cells; ++cell_id)
        {
          // Get cell iterator and compute center on host
          auto               cell = mf_data.get_cell_iterator(color, cell_id);
          Point<dim, double> center_host = cell->center();

          // Get GPU-computed center
          const unsigned int base_idx = global_cell_id * dim;
          Point<dim, double> center_gpu;
          for (unsigned int d = 0; d < dim; ++d)
            center_gpu[d] = centers_gpu_host[base_idx + d];

          if ((center_host - center_gpu).norm() < 1e-8)
            ++n_matches;
          else
            ++n_mismatches;

          global_cell_id++;
        }
    }

  deallog << "Total cells: " << tria.n_active_cells() << std::endl;
  deallog << "Matches: " << n_matches << ", Mismatches: " << n_mismatches
          << std::endl;
}


int
main()
{
  initlog();

  Kokkos::initialize();
  {
    test<2>(1);
    test<2>(2);
    test<2>();
    test<3>(1);
    test<3>(2);
    test<3>();
  }
  Kokkos::finalize();

  return 0;
}

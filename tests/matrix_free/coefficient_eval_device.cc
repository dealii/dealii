// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
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

// Test that we execute the loop in the same order on the CPU and the GPU

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/cuda_matrix_free.templates.h>

#include "../tests.h"

#include "Kokkos_Core.hpp"


template <int dim, int fe_degree>
class DummyOperator
{
public:
  DummyOperator() = default;

  DEAL_II_HOST_DEVICE void
  operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
    CUDAWrappers::SharedData<dim, double> *                     shared_data,
    const double *                                              src,
    double *                                                    dst) const;

  static const unsigned int n_dofs_1d = fe_degree + 1;
  static const unsigned int n_local_dofs =
    dealii::Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points =
    dealii::Utilities::pow(fe_degree + 1, dim);
};



template <int dim, int fe_degree>
DEAL_II_HOST_DEVICE void
DummyOperator<dim, fe_degree>::operator()(
  const unsigned int                                          cell,
  const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
  CUDAWrappers::SharedData<dim, double> *                     shared_data,
  const double *,
  double *dst) const
{
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(shared_data->team_member, n_q_points),
    [&](const int q_point) {
      const unsigned int pos =
        gpu_data->local_q_point_id(cell, n_q_points, q_point);

      auto point = gpu_data->get_quadrature_point(cell, q_point);
      dst[pos] =
        dim == 2 ? point(0) + point(1) : point(0) + point(1) + point(2);
    });
}



template <int dim, int fe_degree>
class DummyMatrixFree : public Subscriptor
{
public:
  DummyMatrixFree(const CUDAWrappers::MatrixFree<dim, double> &data_in,
                  const unsigned int                           size);
  void
  eval(LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &dst)
    const;

private:
  const CUDAWrappers::MatrixFree<dim, double> &data;
};

template <int dim, int fe_degree>
DummyMatrixFree<dim, fe_degree>::DummyMatrixFree(
  const CUDAWrappers::MatrixFree<dim, double> &data_in,
  const unsigned int                           size)
  : data(data_in)
{}


template <int dim, int fe_degree>
void
DummyMatrixFree<dim, fe_degree>::eval(
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &dst) const
{
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> src(dst);
  DummyOperator<dim, fe_degree> dummy_operator;
  data.cell_loop(dummy_operator, src, dst);
}

template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  // Computation on the device
  MappingQ<dim>                         mapping(fe_degree);
  CUDAWrappers::MatrixFree<dim, double> mf_data;
  typename CUDAWrappers::MatrixFree<dim, double>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         update_quadrature_points;
  const QGauss<1> quad(fe_degree + 1);
  mf_data.reinit(mapping, dof, constraints, quad, additional_data);
  constexpr unsigned int n_q_points_per_cell =
    dealii::Utilities::pow(fe_degree + 1, dim);

  const unsigned int              n_dofs = dof.n_dofs();
  DummyMatrixFree<dim, fe_degree> mf(mf_data,
                                     tria.n_active_cells() *
                                       n_q_points_per_cell);
  const unsigned int size = tria.n_active_cells() * n_q_points_per_cell;
  LinearAlgebra::ReadWriteVector<double>                           coef(size);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> coef_device(
    size);

  mf.eval(coef_device);
  Kokkos::fence();
  coef.import_elements(coef_device, VectorOperation::insert);

  // Computation the host
  auto               graph    = mf_data.get_colored_graph();
  unsigned int const n_colors = graph.size();
  for (unsigned int color = 0; color < n_colors; ++color)
    {
      typename CUDAWrappers::MatrixFree<dim, double>::Data gpu_data =
        mf_data.get_data(color);
      unsigned int const n_cells = gpu_data.n_cells;
      auto gpu_data_host = CUDAWrappers::copy_mf_data_to_host<dim, double>(
        gpu_data, additional_data.mapping_update_flags);
      for (unsigned int cell_id = 0; cell_id < n_cells; ++cell_id)
        {
          for (unsigned int i = 0; i < n_q_points_per_cell; ++i)
            {
              unsigned int const pos =
                gpu_data_host.local_q_point_id(cell_id, n_q_points_per_cell, i);
              auto         p = gpu_data_host.get_quadrature_point(cell_id, i);
              const double p_val = dim == 2 ? p(0) + p(1) : p(0) + p(1) + p(2);
              AssertThrow(std::abs(coef[pos] - p_val) < 1e-12,
                          ExcInternalError());
            }
        }
    }
}

int
main()
{
  initlog();
  Kokkos::initialize();

  test<2, 3>();
  test<3, 3>();

  deallog << "OK" << std::endl;
  Kokkos::finalize();
  return 0;
}

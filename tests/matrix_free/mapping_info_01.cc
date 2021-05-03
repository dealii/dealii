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



// Let MatrixFree initialize the geometry for ridiculously high numbers of
// quadrature points and compute the volume and surface of a hyper ball. In an
// initial implementation we used to consume a lot of memory and take a long
// time, whereas the new state is relatively cheap.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria, Point<dim>(), 1.);

  FE_DGQ<dim>     fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  MappingQGeneric<dim> mapping(degree);

  MatrixFree<dim, double>                          mf_data;
  const QGauss<1>                                  quad(degree);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(mapping,
                 std::vector<const DoFHandler<dim> *>{&dof},
                 std::vector<const AffineConstraints<double> *>{&constraints},
                 std::vector<Quadrature<1>>{
                   {QGauss<1>(std::max(degree / 2, 1U)),
                    QGauss<1>(degree + 1),
                    QGauss<1>(3 * degree / 2)}},
                 data);

  deallog << std::setw(5) << degree;
  for (unsigned int index = 0; index < 3; ++index)
    {
      double                volume = 0;
      FEEvaluation<dim, -1> fe_eval(mf_data, 0, index);
      for (unsigned int cell = 0; cell < mf_data.n_cell_batches(); ++cell)
        {
          fe_eval.reinit(cell);
          VectorizedArray<double> local = 0;
          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            local += fe_eval.JxW(q);
          for (unsigned int v = 0;
               v < mf_data.n_active_entries_per_cell_batch(cell);
               ++v)
            volume += local[v];
        }
      deallog << std::setw(11)
              << std::abs((dim == 2 ? numbers::PI : 4. / 3. * numbers::PI) -
                          volume);
    }
  deallog << "  ";
  for (unsigned int index = 0; index < 3; ++index)
    {
      double                    area = 0;
      FEFaceEvaluation<dim, -1> fe_eval(mf_data, true, 0, index);
      for (unsigned int face = mf_data.n_inner_face_batches();
           face <
           mf_data.n_boundary_face_batches() + mf_data.n_inner_face_batches();
           ++face)
        {
          fe_eval.reinit(face);
          VectorizedArray<double> local = 0;
          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            local += fe_eval.JxW(q);
          for (unsigned int v = 0;
               v < mf_data.n_active_entries_per_face_batch(face);
               ++v)
            area += local[v];
        }
      deallog << std::setw(11) << std::abs(2. * (dim - 1) * numbers::PI - area);
    }
  deallog << std::endl;
}



int
main()
{
  initlog();

  deallog << std::setprecision(4);
  deallog << "2D   volume                             area" << std::endl;
  deallog << "deg  err deg/2  err deg+1  err 3deg/2   "
          << "err deg/2  err deg+1  err 3deg/2" << std::endl;
  for (unsigned int q = 1; q < 20; ++q)
    test<2>(q);
  deallog << std::endl;

  deallog << "3D   volume                             area" << std::endl;
  deallog << "deg  err deg/2  err deg+1  err 3deg/2   "
          << "err deg/2  err deg+1  err 3deg/2" << std::endl;
  for (unsigned int q = 1; q < 17; ++q)
    test<3>(q);
  deallog << std::endl;

  return 0;
}

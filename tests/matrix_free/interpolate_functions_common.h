//------------------  interpolate_functions_common.h  ------------------------
//
// Copyright (C) 2018 by the deal.II authors
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
//------------------  interpolate_functions_common.h  ------------------------


// this is a template for getting the function values and comparing them with
// an analytical function for different meshes

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree>
void
test();

template <int dim>
class CompareFunction;



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in,
                 const Mapping<dim> &           mapping)
    : data(data_in){};

  virtual ~MatrixFreeTest()
  {}

  // make function virtual to allow derived classes to define a different
  // function
  virtual void
  cell(const MatrixFree<dim, Number> &data,
       Vector<Number> &,
       const Vector<Number> &                       src,
       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(data);

    CompareFunction<dim> function;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, true, true);

        for (unsigned int j = 0; j < data.n_components_filled(cell); ++j)
          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            {
              ++cell_times;
              Point<dim> p;
              for (unsigned int d = 0; d < dim; ++d)
                p[d] = fe_eval.quadrature_point(q)[d][j];
              cell_errors[0] +=
                std::abs(fe_eval.get_value(q)[j] - function.value(p, 0));
              for (unsigned int d = 0; d < dim; ++d)
                cell_errors[1] += std::abs(fe_eval.get_gradient(q)[d][j] -
                                           function.gradient(p, 0)[d]);
              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  cell_errors[2] += std::abs(fe_eval.get_hessian(q)[d][e][j] -
                                             function.hessian(p, 0)[d][e]);
            }
      }
  }

  virtual void
  face(const MatrixFree<dim, Number> &data,
       Vector<Number> &,
       const Vector<Number> &                       src,
       const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_evalm(data,
                                                                        true);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_evalp(data,
                                                                        false);

    CompareFunction<dim> function;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_evalm.reinit(face);
        fe_evalm.read_dof_values(src);
        fe_evalm.evaluate(true, true);
        fe_evalp.reinit(face);
        fe_evalp.read_dof_values(src);
        fe_evalp.evaluate(true, true);

        for (unsigned int j = 0; j < VectorizedArray<Number>::size(); ++j)
          {
            // skip empty components in VectorizedArray
            if (data.get_face_info(face).cells_interior[j] ==
                numbers::invalid_unsigned_int)
              break;
            for (unsigned int q = 0; q < fe_evalm.n_q_points; ++q)
              {
                ++facem_times;
                ++facep_times;
                Point<dim> p;
                for (unsigned int d = 0; d < dim; ++d)
                  p[d] = fe_evalm.quadrature_point(q)[d][j];
                facem_errors[0] +=
                  std::abs(fe_evalm.get_value(q)[j] - function.value(p, 0));
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    // std::cout << fe_evalm.get_gradient(q)[d][j] << " ";
                    facem_errors[1] += std::abs(fe_evalm.get_gradient(q)[d][j] -
                                                function.gradient(p, 0)[d]);
                  }
                double normal_derivative = 0;
                for (unsigned int d = 0; d < dim; ++d)
                  normal_derivative += function.gradient(p, 0)[d] *
                                       fe_evalm.get_normal_vector(q)[d][j];
                facem_errors[3] += std::abs(
                  fe_evalm.get_normal_derivative(q)[j] - normal_derivative);
                facep_errors[0] +=
                  std::abs(fe_evalp.get_value(q)[j] - function.value(p, 0));
                for (unsigned int d = 0; d < dim; ++d)
                  facep_errors[1] += std::abs(fe_evalp.get_gradient(q)[d][j] -
                                              function.gradient(p, 0)[d]);
                facep_errors[3] += std::abs(
                  fe_evalp.get_normal_derivative(q)[j] - normal_derivative);
                // hessians not currently implemented in FEFaceEvaluation
                // for (unsigned int d=0; d<dim; ++d)
                //  for (unsigned int e=0; e<dim; ++e)
                //    facem_errors[2] +=
                //    std::abs(fe_evalm.get_hessian(q)[d][e][j]-
                //                                function.hessian(q)[d][e]);
              }
            // std::cout << std::endl;
          }
      }
  }

  virtual void
  boundary(const MatrixFree<dim, Number> &data,
           Vector<Number> &,
           const Vector<Number> &                       src,
           const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_evalm(data,
                                                                        true);

    CompareFunction<dim> function;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_evalm.reinit(face);
        fe_evalm.read_dof_values(src);
        fe_evalm.evaluate(true, true);

        for (unsigned int j = 0; j < VectorizedArray<Number>::size(); ++j)
          {
            // skip empty components in VectorizedArray
            if (data.get_face_info(face).cells_interior[j] ==
                numbers::invalid_unsigned_int)
              break;
            for (unsigned int q = 0; q < fe_evalm.n_q_points; ++q)
              {
                ++boundary_times;
                Point<dim> p;
                for (unsigned int d = 0; d < dim; ++d)
                  p[d] = fe_evalm.quadrature_point(q)[d][j];
                boundary_errors[0] +=
                  std::abs(fe_evalm.get_value(q)[j] - function.value(p, 0));
                for (unsigned int d = 0; d < dim; ++d)
                  boundary_errors[1] +=
                    std::abs(fe_evalm.get_gradient(q)[d][j] -
                             function.gradient(p, 0)[d]);
                double normal_derivative = 0;
                for (unsigned int d = 0; d < dim; ++d)
                  normal_derivative += function.gradient(p, 0)[d] *
                                       fe_evalm.get_normal_vector(q)[d][j];
                boundary_errors[3] += std::abs(
                  fe_evalm.get_normal_derivative(q)[j] - normal_derivative);
                // hessians not currently implemented in FEFaceEvaluation
                // for (unsigned int d=0; d<dim; ++d)
                //  for (unsigned int e=0; e<dim; ++e)
                //    boundary_errors[2] +=
                //    std::abs(fe_evalm.get_hessian(q)[d][e][j]-
                //                                   function.hessian(q)[d][e]);
              }
          }
      }
  }



  void
  test_functions(const Vector<Number> &src) const
  {
    for (unsigned int i = 0; i < 4; ++i)
      {
        cell_errors[i]     = 0;
        facem_errors[i]    = 0;
        facep_errors[i]    = 0;
        boundary_errors[i] = 0;
      }
    cell_times = facem_times = facep_times = boundary_times = 0;

    Vector<Number> dst_dummy;
    data.loop(&MatrixFreeTest::cell,
              &MatrixFreeTest::face,
              &MatrixFreeTest::boundary,
              this,
              dst_dummy,
              src);

    if (std::is_same<Number, float>::value)
      for (unsigned int i = 0; i < 4; ++i)
        {
          if (cell_errors[i] / cell_times < 1e-5)
            cell_errors[i] = 0;
          if (facem_errors[i] / cell_times < 1e-5)
            facem_errors[i] = 0;
          if (facep_errors[i] / cell_times < 1e-5)
            facep_errors[i] = 0;
          if (boundary_errors[i] / cell_times < 1e-5)
            boundary_errors[i] = 0;
        }
    deallog << "Error cell values:      " << cell_errors[0] / cell_times
            << std::endl;
    deallog << "Error cell gradients:   " << cell_errors[1] / cell_times
            << std::endl;
    deallog << "Error cell Hessians:    " << cell_errors[2] / cell_times
            << std::endl;
    deallog << "Error face- values:     " << facem_errors[0] / facem_times
            << std::endl;
    deallog << "Error face- gradients:  " << facem_errors[1] / facem_times
            << std::endl;
    // deallog << "Error face- Hessians:   " << facem_errors[2]/facem_times <<
    // std::endl;
    deallog << "Error face- grad*normal:" << facem_errors[3] / facem_times
            << std::endl;
    deallog << "Error face+ values:     " << facep_errors[0] / facep_times
            << std::endl;
    deallog << "Error face+ gradients:  " << facep_errors[1] / facep_times
            << std::endl;
    // deallog << "Error face+ Hessians:   " << facep_errors[2]/facep_times <<
    // std::endl;
    deallog << "Error face+ grad*normal:" << facem_errors[3] / facep_times
            << std::endl;
    deallog << "Error face b values:    " << boundary_errors[0] / boundary_times
            << std::endl;
    deallog << "Error face b gradients: " << boundary_errors[1] / boundary_times
            << std::endl;
    // deallog << "Error face b Hessians:  " <<
    // boundary_errors[2]/boundary_times << std::endl;
    deallog << "Error face b grad*norm: " << boundary_errors[3] / boundary_times
            << std::endl;
    deallog << std::endl;
  };

protected:
  const MatrixFree<dim, Number> &data;
  mutable double cell_errors[4], facem_errors[4], facep_errors[4],
    boundary_errors[4];
  mutable unsigned long long cell_times, facem_times, facep_times,
    boundary_times;
};



template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim> &          dof,
        const AffineConstraints<double> &constraints)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  // use this for info on problem
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells()
  //          << std::endl;
  // std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  // std::cout << "Number of constraints: " << constraints.n_constraints() <<
  // std::endl;

  Vector<number> interpolated(dof.n_dofs());
  VectorTools::interpolate(dof, CompareFunction<dim>(), interpolated);

  constraints.distribute(interpolated);
  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1> quad(dof.get_fe().degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.mapping_update_flags =
      update_gradients | update_second_derivatives | update_quadrature_points;
    data.mapping_update_flags_boundary_faces =
      update_gradients | update_second_derivatives | update_quadrature_points;
    data.mapping_update_flags_inner_faces =
      update_gradients | update_second_derivatives | update_quadrature_points;
    mf_data.reinit(dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, fe_degree + 1, number> mf(mf_data);
  mf.test_functions(interpolated);
}


int
main()
{
  initlog();

  deallog << std::setprecision(3);
  {
    deallog.push("1d");
    test<1, 0>();
    test<1, 2>();
    deallog.pop();
    deallog.push("2d");
    test<2, 0>();
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}

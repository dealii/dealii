// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this function tests the correctness of the implementation of matrix free
// operations in getting the function values, the function gradients, and the
// function Laplacians on a cartesian mesh (hyper cube). This tests whether
// cartesian meshes are treated correctly. The test case is without any
// constraints

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_raviart_thomas.h>
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


template <int dim>
class CompareFunction : public Function<dim>
{
public:
  CompareFunction()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    double value = (1.2 - 0.1 * component) * p[0] + 0.4;
    for (unsigned int d = 1; d < dim; ++d)
      value -= (2.7 - 0.2 * component) * d * p[d];
    return value;
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component) const
  {
    Tensor<1, dim> grad;
    grad[0] = 1.2 - 0.1 * component;
    for (unsigned int d = 1; d < dim; ++d)
      grad[d] = -(2.7 - 0.2 * component) * d;
    return grad;
  }
};



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
                 const Mapping<dim>            &mapping)
    : data(data_in){};

  virtual ~MatrixFreeTest()
  {}

  // make function virtual to allow derived classes to define a different
  // function
  virtual void
  cell(const MatrixFree<dim, Number> &data,
       Vector<Number> &,
       const Vector<Number>                        &src,
       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_eval(data);

    CompareFunction<dim> function;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        for (unsigned int j = 0; j < data.n_active_entries_per_cell_batch(cell);
             ++j)
          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            {
              ++cell_times;
              Point<dim> p;
              for (unsigned int d = 0; d < dim; ++d)
                p[d] = fe_eval.quadrature_point(q)[d][j];
              for (unsigned int d = 0; d < dim; ++d)
                {
                  cell_errors[0][d] +=
                    std::abs(fe_eval.get_value(q)[d][j] - function.value(p, d));
                  for (unsigned int e = 0; e < dim; ++e)
                    cell_errors[1][d] +=
                      std::abs(fe_eval.get_gradient(q)[d][e][j] -
                               function.gradient(p, d)[e]);
                }
              double divergence = 0;
              for (unsigned int d = 0; d < dim; ++d)
                divergence += function.gradient(p, d)[d];
              cell_errors[2][0] +=
                std::abs(fe_eval.get_divergence(q)[j] - divergence);
            }
      }
  }

  virtual void
  face(const MatrixFree<dim, Number> &data,
       Vector<Number> &,
       const Vector<Number>                        &src,
       const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_evalm(data,
                                                                          true);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_evalp(
      data, false);

    CompareFunction<dim> function;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_evalm.reinit(face);
        fe_evalm.read_dof_values(src);
        fe_evalm.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                          EvaluationFlags::hessians);
        fe_evalp.reinit(face);
        fe_evalp.read_dof_values(src);
        fe_evalp.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                          EvaluationFlags::hessians);

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

                // interior face
                for (unsigned int d = 0; d < dim; ++d)
                  p[d] = fe_evalm.quadrature_point(q)[d][j];
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    facem_errors[0][d] += std::abs(fe_evalm.get_value(q)[d][j] -
                                                   function.value(p, d));
                    for (unsigned int e = 0; e < dim; ++e)
                      {
                        facem_errors[1][d] +=
                          std::abs(fe_evalm.get_gradient(q)[d][e][j] -
                                   function.gradient(p, d)[e]);
                      }
                  }
                double divergence = 0;
                for (unsigned int d = 0; d < dim; ++d)
                  divergence += function.gradient(p, d)[d];
                facem_errors[2][0] +=
                  std::abs(fe_evalm.get_divergence(q)[j] - divergence);

                // exterior face
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    facep_errors[0][d] += std::abs(fe_evalp.get_value(q)[d][j] -
                                                   function.value(p, d));
                    for (unsigned int e = 0; e < dim; ++e)
                      facep_errors[1][d] +=
                        std::abs(fe_evalp.get_gradient(q)[d][e][j] -
                                 function.gradient(p, d)[e]);
                  }
                facep_errors[2][0] +=
                  std::abs(fe_evalp.get_divergence(q)[j] - divergence);
              }
          }
      }
  }

  virtual void
  boundary(const MatrixFree<dim, Number> &data,
           Vector<Number> &,
           const Vector<Number>                        &src,
           const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_evalm(data,
                                                                          true);

    CompareFunction<dim> function;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_evalm.reinit(face);
        fe_evalm.read_dof_values(src);
        fe_evalm.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                          EvaluationFlags::hessians);

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
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    boundary_errors[0][d] += std::abs(
                      fe_evalm.get_value(q)[d][j] - function.value(p, d));
                    for (unsigned int e = 0; e < dim; ++e)
                      boundary_errors[1][d] +=
                        std::abs(fe_evalm.get_gradient(q)[d][e][j] -
                                 function.gradient(p, d)[e]);
                  }
                double divergence = 0;
                for (unsigned int d = 0; d < dim; ++d)
                  divergence += function.gradient(p, d)[d];
                boundary_errors[2][0] +=
                  std::abs(fe_evalm.get_divergence(q)[j] - divergence);
              }
          }
      }
  }



  static void
  print_error(const std::string                     &text,
              const dealii::ndarray<double, 3, dim> &array,
              const unsigned long long               n_times)
  {
    deallog << "Error " << std::left << std::setw(6) << text << " values:     ";
    for (unsigned int d = 0; d < dim; ++d)
      deallog << array[0][d] / n_times << " ";
    deallog << std::endl;
    deallog << "Error " << std::left << std::setw(6) << text << " gradients:  ";
    for (unsigned int d = 0; d < dim; ++d)
      deallog << array[1][d] / n_times << " ";
    deallog << std::endl;
    deallog << "Error " << std::left << std::setw(6) << text << " divergence: ";
    deallog << array[2][0] / n_times << " ";
    deallog << std::endl;
  }

  void
  test_functions(const Vector<Number> &src) const
  {
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int i = 0; i < 3; ++i)
        {
          cell_errors[i][d]     = 0;
          facem_errors[i][d]    = 0;
          facep_errors[i][d]    = 0;
          boundary_errors[i][d] = 0;
        }
    cell_times = facem_times = facep_times = boundary_times = 0;

    Vector<Number> dst_dummy;
    data.loop(&MatrixFreeTest::cell,
              &MatrixFreeTest::face,
              &MatrixFreeTest::boundary,
              this,
              dst_dummy,
              src);

    print_error("cell", cell_errors, cell_times);
    print_error("face-", facem_errors, facem_times);
    print_error("face+", facep_errors, facep_times);
    print_error("face b", boundary_errors, boundary_times);
    deallog << std::endl;
  };

protected:
  const MatrixFree<dim, Number>          &data;
  mutable dealii::ndarray<double, 3, dim> cell_errors, facem_errors,
    facep_errors, boundary_errors;
  mutable unsigned long long cell_times, facem_times, facep_times,
    boundary_times;
};



template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim>           &dof,
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
    data.mapping_update_flags  = update_gradients | update_quadrature_points;
    data.mapping_update_flags_boundary_faces =
      update_gradients | update_quadrature_points;
    data.mapping_update_flags_inner_faces =
      update_gradients | update_quadrature_points;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, fe_degree + 1, number> mf(mf_data);
  mf.test_functions(interpolated);
}



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  // 1 refinement will make sure that the Piola transform is 1
  tria.refine_global(1);

  {
    FE_RaviartThomasNodal<dim> fe(fe_degree - 1);
    DoFHandler<dim>            dof(tria);
    dof.distribute_dofs(fe);

    AffineConstraints<double> constraints;
    constraints.close();
    if (fe_degree > 2)
      do_test<dim, -1, double>(dof, constraints);
    else
      do_test<dim, fe_degree, double>(dof, constraints);
  }
}



int
main()
{
  initlog();

  deallog << std::setprecision(5);
  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    test<3, 3>();
    deallog.pop();
  }
}

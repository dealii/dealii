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

#include <deal.II/fe/fe_nedelec.h>
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
  bool constant = false;
  CompareFunction()
    : Function<dim>(dim)
  {}

  CompareFunction(bool constant)
    : Function<dim>(dim)
  {
    this->constant = constant;
  }

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    if (constant)
      {
        return 23.4 + component;
      }
    double value;
    value = (1.2 - 0.1 * component) * p[0] + 0.4;
    for (unsigned int d = 1; d < dim; ++d)
      value -= (2.7 - 0.2 * component) * d * p[d];
    return value;
  }


  virtual Tensor<1, (dim == 2 ? 1 : dim)>
  curl(const Point<dim> &p) const
  {
    Tensor<1, dim == 2 ? 1 : dim> curl;
    if (dim == 2)
      {
        if (constant)
          {
            curl[0] = 0;
          }
        else
          curl[0] = (1.2 - 0.1) + 2.7;
        return curl;
      }
    if (constant)
      {
        curl[0] = curl[1] = curl[2] = 0;
      }
    else
      {
        curl[0] = -(2.7 - 0.2 * 2) + (2.7 - 0.2) * 2;
        curl[1] = -2.7 * 2 - (1.2 - 0.1 * 2);
        curl[2] = (1.2 - 0.1) + 2.7;
      }
    return curl;
  }
};

// quadratic function for comparison
template <int dim>
class CompareFunctionQuad : public Function<dim>
{
public:
  bool constant = false;
  CompareFunctionQuad()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    double value;
    value = (1.2 - 0.1 * component) * p[0] * p[0] + 0.4;
    for (unsigned int d = 1; d < dim; ++d)
      value -= (2.7 - 0.2 * component) * d * p[d] * p[d];
    return value;
  }


  virtual Tensor<1, (dim == 2 ? 1 : dim)>
  curl(const Point<dim> &p) const
  {
    Tensor<1, dim == 2 ? 1 : dim> curl;
    if (dim == 2)
      {
        curl[0] = (1.2 - 0.1) * 2 * p[0] + 2.7 * 2 * p[1];
        return curl;
      }

    curl[0] = -(2.7 - 0.2 * 2) * 2 * p[1] + (2.7 - 0.2) * 4 * p[2];
    curl[1] = -2.7 * 4 * p[2] - (1.2 - 0.1 * 2) * 2 * p[0];
    curl[2] = (1.2 - 0.1) * 2 * p[0] + 2.7 * 2 * p[1];

    return curl;
  }
};



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 2,
          typename Number   = double>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in, bool Constant = false)
    : data(data_in)
    , constFunction(Constant){};

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in,
                 const Mapping<dim>            &mapping,
                 bool                           Constant = false)
    : data(data_in)
    , constFunction(Constant){};

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

    CompareFunction<dim> function(constFunction);

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
                }

              for (unsigned int d = 0; d < (dim == 2 ? 1 : dim); d++)
                {
                  cell_errors[1][d] +=
                    std::abs(fe_eval.get_curl(q)[d][j] - function.curl(p)[d]);
                }
            }
      }
  }



  static void
  print_error(const std::string                     &text,
              const dealii::ndarray<double, 2, dim> &array,
              const unsigned long long               n_times)
  {
    deallog << "Error " << std::left << std::setw(6) << text << " values:     ";
    for (unsigned int d = 0; d < dim; ++d)
      deallog << array[0][d] / n_times << " ";
    deallog << std::endl;

    deallog << "Error " << std::left << std::setw(6) << text << " curl: ";
    for (unsigned int d = 0; d < dim; ++d)
      deallog << array[1][d] / n_times << " ";
    deallog << std::endl;
  }

  void
  test_functions(const Vector<Number> &src) const
  {
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int i = 0; i < 2; ++i)
        {
          cell_errors[i][d] = 0;
        }
    cell_times = 0;

    Vector<Number> dst_dummy;
    data.cell_loop(&MatrixFreeTest::cell, this, dst_dummy, src);

    print_error("cell", cell_errors, cell_times);
    deallog << std::endl;
  };

protected:
  const MatrixFree<dim, Number>          &data;
  mutable dealii::ndarray<double, 2, dim> cell_errors;
  mutable unsigned long long              cell_times;
  bool                                    constFunction;
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
  Vector<number> interpolatedConst(dof.n_dofs());
  VectorTools::interpolate(dof, CompareFunction<dim>(false), interpolated);
  VectorTools::interpolate(dof, CompareFunction<dim>(true), interpolatedConst);

  constraints.distribute(interpolated);
  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1> quad(dof.get_fe().degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.mapping_update_flags  = update_gradients | update_quadrature_points;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, fe_degree + 2, number> mf(mf_data);
  MatrixFreeTest<dim, fe_degree, fe_degree + 2, number> mfConstFunction(mf_data,
                                                                        true);
  mf.test_functions(interpolated);
  mfConstFunction.test_functions(interpolatedConst);
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
    FE_NedelecNodal<dim> fe(fe_degree);
    DoFHandler<dim>      dof(tria);
    dof.distribute_dofs(fe);

    AffineConstraints<double> constraints;
    constraints.close();
    if (fe_degree > 1)
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
    test<2, 0>();
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    deallog.pop();
    deallog.push("3d");
    test<3, 0>();
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}

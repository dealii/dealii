// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// test FE_Enriched class by comparing with FE_System for multiple enrichment
// functions.

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <iostream>

#include "../tests.h"

const double eps = 1e-10;

// create a mesh without hanging nodes,
// move nodes randomly
// and compare FEEnriched to FESystem with explicit treatment of
// the product rule for
// valus, gradients, hessians on
// elements and faces.
// The comparison is straight forward because local dofs are enumerated
// in the same way for FE_System and FEEnriched.
// Same as fe_enriched_compare_to_fe_system.cc but now multiple enrichment
// functions.
const unsigned int patches = 10;


template <int dim>
class EnrichmentFunction : public Function<dim>
{
public:
  EnrichmentFunction(const Point<dim> &origin)
    : Function<dim>(1)
    , origin(origin)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> d = point - origin;
    return std::exp(-d.norm());
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> d = point - origin;
    Assert(d.norm() > 0, dealii::ExcMessage("gradient is not defined at zero"));
    d *= -std::exp(-d.norm()) / d.norm();
    return d;
  }

  virtual SymmetricTensor<2, dim>
  hessian(const Point<dim> &p, const unsigned int component = 0) const
  {
    Tensor<1, dim> dir = p - origin;
    const double   r   = dir.norm();
    Assert(r > 0.0, ExcMessage("r is not positive"));
    dir /= r;
    SymmetricTensor<2, dim> dir_x_dir;
    for (unsigned int i = 0; i < dim; i++)
      for (unsigned int j = i; j < dim; j++)
        dir_x_dir[i][j] = dir[i] * dir[j];

    return std::exp(-r) *
           ((1.0 + 0.5 / r) * dir_x_dir - unit_symmetric_tensor<dim>() / r);
  }

private:
  Point<dim> origin;
};

/**
 *
 * @param p     point
 * @param func  function
 * @param v_e   value of enriched FE
 * @param g_e   gradient of enriched FE
 * @param h_e   hessian of enriched FE
 * @param v_s0  value of 0th component of FE_System
 * @param g_s0  gradient of 0th component of FE_System
 * @param h_s0  hessian of 0th component of FE_System
 * @param v_s1  value of 1st component of FE_System
 * @param g_s1  gradient of 1st component of FE_System
 * @param h_s1  hessian of 1st component of FE_System
 */
template <int dim>
void
check_consistency(const Point<dim> &    p,
                  const Function<dim> & func1,
                  const Function<dim> & func2,
                  const Function<dim> & func3,
                  const double &        v_e,
                  const Tensor<1, dim> &g_e,
                  const Tensor<2, dim> &h_e,
                  const double &        v_s0,
                  const Tensor<1, dim> &g_s0,
                  const Tensor<2, dim> &h_s0,
                  const double &        v_s1,
                  const Tensor<1, dim> &g_s1,
                  const Tensor<2, dim> &h_s1,
                  const double &        v_s2,
                  const Tensor<1, dim> &g_s2,
                  const Tensor<2, dim> &h_s2,
                  const double &        v_s3,
                  const Tensor<1, dim> &g_s3,
                  const Tensor<2, dim> &h_s3)
{
  const double                  v_f1 = func1.value(p);
  const Tensor<1, dim>          g_f1 = func1.gradient(p);
  const SymmetricTensor<2, dim> h_f1 = func1.hessian(p);

  const double                  v_f2 = func2.value(p);
  const Tensor<1, dim>          g_f2 = func2.gradient(p);
  const SymmetricTensor<2, dim> h_f2 = func2.hessian(p);

  const double                  v_f3 = func3.value(p);
  const Tensor<1, dim>          g_f3 = func3.gradient(p);
  const SymmetricTensor<2, dim> h_f3 = func3.hessian(p);



  // product rule:
  const double         v_s = v_s0 + v_s1 * v_f1 + v_s2 * v_f2 + v_s3 * v_f3;
  const Tensor<1, dim> g_s = g_s0 + g_s1 * v_f1 + v_s1 * g_f1 + g_s2 * v_f2 +
                             v_s2 * g_f2 + g_s3 * v_f3 + v_s3 * g_f3;
  Tensor<2, dim> op1 = outer_product(g_s1, g_f1),
                 op2 = outer_product(g_s2, g_f2),
                 op3 = outer_product(g_s3, g_f3);

  const SymmetricTensor<2, dim> sp1 =
    symmetrize(op1) * 2.0; // symmetrize does [s+s^T]/2
  const SymmetricTensor<2, dim> sp2 =
    symmetrize(op2) * 2.0; // symmetrize does [s+s^T]/2
  const SymmetricTensor<2, dim> sp3 =
    symmetrize(op3) * 2.0; // symmetrize does [s+s^T]/2
  // Hessians should be symmetric, but due to round off errors and the fact
  // that deal.II stores hessians in tensors, we may end up with failing check
  // for symmetry in SymmetricTensor class.
  // For a moment stick with usual Tensor
  const Tensor<2, dim> h_s = h_s0 + h_s1 * v_f1 + sp1 + v_s1 * h_f1 +
                             h_s2 * v_f2 + sp2 + v_s2 * h_f2 + h_s3 * v_f3 +
                             sp3 + v_s3 * h_f3;

  const double v_d = v_s - v_e;
  AssertThrow(std::abs(v_d) < eps, ExcInternalError());
  const Tensor<1, dim> g_d = g_s - g_e;
  AssertThrow(g_d.norm() < eps, ExcInternalError());

  // see note above.
  const Tensor<2, dim> h_d = h_s - h_e;
  AssertThrow(h_d.norm() < eps, ExcInternalError());
}

/**
 * {fe_en1,1}, {fe_en2,2} => 3 functions
 */
template <int dim>
void
test(const FiniteElement<dim> & fe_base,
     const FiniteElement<dim> & fe_en1,
     const FiniteElement<dim> & fe_en2,
     const Quadrature<dim> &    volume_quad,
     const Quadrature<dim - 1> &face_quad,
     const bool                 distort)
{
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler_enriched(triangulation);
  DoFHandler<dim>    dof_handler_system(triangulation);

  Point<dim> p1, p2, p3;
  for (unsigned int d = 0; d < dim; d++)
    {
      p1[d] = 0.0;
      p2[d] = 0.0;
      p3[d] = 0.0;
    }
  p2[0] = 10.0;
  p3[0] = 5.0;
  EnrichmentFunction<dim>                 fun1(p1), fun2(p2), fun3(p3);
  std::vector<const FiniteElement<dim> *> fe_enrichements(3);
  fe_enrichements[0] = &fe_en1;
  fe_enrichements[1] = &fe_en2;
  FE_Nothing<dim> fe_nothing(1, true);
  fe_enrichements[2] = &fe_nothing; // should be ignored
  std::vector<std::vector<std::function<const Function<dim> *(
    const typename Triangulation<dim, dim>::cell_iterator &)>>>
    functions(3);
  functions[0].resize(1);
  functions[0][0] =
    [&](const typename Triangulation<dim, dim>::cell_iterator &) {
      return &fun1;
    };
  functions[1].resize(2);
  functions[1][0] =
    [&](const typename Triangulation<dim, dim>::cell_iterator &) {
      return &fun2;
    };
  functions[1][1] =
    [&](const typename Triangulation<dim, dim>::cell_iterator &) {
      return &fun3;
    };
  functions[2].resize(1);
  functions[2][0] = [&](
                      const typename Triangulation<dim, dim>::cell_iterator &) {
    AssertThrow(false, ExcMessage("Called enrichment function for FE_Nothing"));
    return nullptr;
  };

  FE_Enriched<dim> fe_enriched(&fe_base, fe_enrichements, functions);
  FESystem<dim>    fe_system(fe_base, 1, fe_en1, 1, fe_en2, 2);

  {
    Point<dim>                p1, p2;
    std::vector<unsigned int> repetitions(dim);
    for (unsigned int d = 0; d < dim; d++)
      {
        p1[d]          = -1.0;
        p2[d]          = 2.0;
        repetitions[d] = 3;
      }
    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              repetitions,
                                              p1,
                                              p2);

    if (distort)
      GridTools::distort_random(0.1, triangulation);
  }

  dof_handler_enriched.distribute_dofs(fe_enriched);
  dof_handler_system.distribute_dofs(fe_system);

  FEValues<dim> fe_values_enriched(fe_enriched,
                                   volume_quad,
                                   update_values | update_gradients |
                                     update_hessians);

  FEValues<dim> fe_values_system(fe_system,
                                 volume_quad,
                                 update_values | update_gradients |
                                   update_hessians | update_quadrature_points);

  FEFaceValues<dim> fe_face_values_enriched(
    fe_enriched, face_quad, update_values | update_gradients | update_hessians);

  FEFaceValues<dim> fe_face_values_system(fe_system,
                                          face_quad,
                                          update_values | update_gradients |
                                            update_hessians |
                                            update_quadrature_points);

  AssertThrow(dof_handler_enriched.n_dofs() == dof_handler_system.n_dofs(),
              ExcInternalError());
  Vector<double> solution(dof_handler_enriched.n_dofs());

  const unsigned int dofs_per_cell = fe_enriched.dofs_per_cell;
  Assert(fe_enriched.dofs_per_cell == fe_system.dofs_per_cell,
         ExcInternalError());

  typename DoFHandler<dim>::active_cell_iterator
    cell_enriched = dof_handler_enriched.begin_active(),
    endc_enriched = dof_handler_enriched.end(),
    cell_system   = dof_handler_system.begin_active(),
    endc_system   = dof_handler_system.end();
  for (; cell_enriched != endc_enriched; ++cell_enriched, ++cell_system)
    {
      fe_values_enriched.reinit(cell_enriched);
      fe_values_system.reinit(cell_system);
      const unsigned int                     n_q_points = volume_quad.size();
      const std::vector<dealii::Point<dim>> &q_points =
        fe_values_system.get_quadrature_points();

      // check shape functions
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          check_consistency(
            q_points[q_point],
            fun1,
            fun2,
            fun3,
            fe_values_enriched.shape_value(i, q_point),
            fe_values_enriched.shape_grad(i, q_point),
            fe_values_enriched.shape_hessian(i, q_point),
            fe_values_system.shape_value_component(i, q_point, 0),
            fe_values_system.shape_grad_component(i, q_point, 0),
            fe_values_system.shape_hessian_component(i, q_point, 0),
            fe_values_system.shape_value_component(i, q_point, 1),
            fe_values_system.shape_grad_component(i, q_point, 1),
            fe_values_system.shape_hessian_component(i, q_point, 1),
            fe_values_system.shape_value_component(i, q_point, 2),
            fe_values_system.shape_grad_component(i, q_point, 2),
            fe_values_system.shape_hessian_component(i, q_point, 2),
            fe_values_system.shape_value_component(i, q_point, 3),
            fe_values_system.shape_grad_component(i, q_point, 3),
            fe_values_system.shape_hessian_component(i, q_point, 3));

      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        {
          fe_face_values_enriched.reinit(cell_enriched, face);
          fe_face_values_system.reinit(cell_system, face);
          const unsigned int n_q_points_face = face_quad.size();
          const std::vector<dealii::Point<dim>> &q_points =
            fe_face_values_system.get_quadrature_points();

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int q_point = 0; q_point < n_q_points_face; ++q_point)
              check_consistency(
                q_points[q_point],
                fun1,
                fun2,
                fun3,
                fe_face_values_enriched.shape_value(i, q_point),
                fe_face_values_enriched.shape_grad(i, q_point),
                fe_face_values_enriched.shape_hessian(i, q_point),
                fe_face_values_system.shape_value_component(i, q_point, 0),
                fe_face_values_system.shape_grad_component(i, q_point, 0),
                fe_face_values_system.shape_hessian_component(i, q_point, 0),
                fe_face_values_system.shape_value_component(i, q_point, 1),
                fe_face_values_system.shape_grad_component(i, q_point, 1),
                fe_face_values_system.shape_hessian_component(i, q_point, 1),
                fe_face_values_system.shape_value_component(i, q_point, 2),
                fe_face_values_system.shape_grad_component(i, q_point, 2),
                fe_face_values_system.shape_hessian_component(i, q_point, 2),
                fe_face_values_system.shape_value_component(i, q_point, 3),
                fe_face_values_system.shape_grad_component(i, q_point, 3),
                fe_face_values_system.shape_hessian_component(i, q_point, 3));
        }
    }


  deallog << "Ok" << std::endl;
  dof_handler_enriched.clear();
  dof_handler_system.clear();
}



int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4) << std::fixed;
  deallog.depth_console(0);


  try
    {
      {
        const unsigned int dim = 2;
        test(FE_Q<dim>(3),
             FE_Q<dim>(2),
             FE_Q<dim>(1),
             QGauss<dim>(2),
             QGauss<dim - 1>(2),
             false);
        test(FE_Q<dim>(3),
             FE_Q<dim>(2),
             FE_Q<dim>(1),
             QGauss<dim>(2),
             QGauss<dim - 1>(2),
             true);
      }

      {
        const unsigned int dim = 3;
        test(FE_Q<dim>(3),
             FE_Q<dim>(2),
             FE_Q<dim>(1),
             QGauss<dim>(2),
             QGauss<dim - 1>(2),
             false);
        test(FE_Q<dim>(3),
             FE_Q<dim>(2),
             FE_Q<dim>(1),
             QGauss<dim>(2),
             QGauss<dim - 1>(2),
             true);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// test kinematic tensor definitions

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/transformations.h>

#include "../tests.h"


using namespace dealii::Physics;
using namespace dealii::Physics::Elasticity;

template <int dim>
void
test_kinematic_tensors()
{
  const FESystem<dim> fe(FE_Q<dim>(1), dim);
  const QGauss<dim>   qf(2);
  Triangulation<dim>  tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> soln_t(dof_handler.n_dofs());
  Vector<double> soln_t1(dof_handler.n_dofs());

  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        if (std::abs(cell->vertex(v)[0] - 1.0) < 1e-9)
          soln_t[cell->vertex_dof_index(v, 0)] = 1.0;
    }

  const double   delta_t    = 2.0;
  Vector<double> dot_soln_t = soln_t;
  dot_soln_t -= soln_t1;
  dot_soln_t *= (1.0 / delta_t);

  FEValuesExtractors::Vector  u_fe(0);
  std::vector<Tensor<2, dim>> qp_Grad_u_t;
  std::vector<Tensor<2, dim>> qp_Grad_u_t1;
  std::vector<Tensor<2, dim>> qp_dot_Grad_u_t;
  std::vector<Tensor<2, dim>> qp_dot_grad_u_t;

  FEValues<dim>         fe_values(fe, qf, update_gradients);
  MappingQEulerian<dim> q1_mapping(1, dof_handler, soln_t);
  FEValues<dim>         fe_values_mapped(q1_mapping, fe, qf, update_gradients);

  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      fe_values.reinit(cell);
      fe_values_mapped.reinit(cell);

      const unsigned int n_q_points = fe_values.get_quadrature().size();
      qp_Grad_u_t.resize(n_q_points);
      qp_Grad_u_t1.resize(n_q_points);
      qp_dot_Grad_u_t.resize(n_q_points);
      qp_dot_grad_u_t.resize(n_q_points);

      fe_values[u_fe].get_function_gradients(soln_t, qp_Grad_u_t);
      fe_values[u_fe].get_function_gradients(soln_t1, qp_Grad_u_t1);
      fe_values[u_fe].get_function_gradients(dot_soln_t, qp_dot_Grad_u_t);

      fe_values_mapped[u_fe].get_function_gradients(dot_soln_t,
                                                    qp_dot_grad_u_t);

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          static const double tol = 1e-12;

          // Material gradients
          const Tensor<2, dim> &Grad_u    = qp_Grad_u_t[q_point];
          const Tensor<2, dim> &Grad_u_t1 = qp_Grad_u_t1[q_point];

          // --- Rate independent ---

          // Deformation gradient tensor
          const Tensor<2, dim> F_t1 = Kinematics::F(Grad_u_t1);
          Assert((F_t1 - unit_symmetric_tensor<dim>()).norm() < tol,
                 ExcMessage("Incorrect computation of F_t1"));
          const Tensor<2, dim> F = Kinematics::F(Grad_u);
          Assert((F -
                  (static_cast<Tensor<2, dim>>(unit_symmetric_tensor<dim>()) +
                   Grad_u))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of F"));

          // Volumetric / isochoric split of deformation gradient
          Assert(determinant(F) != 1.0,
                 ExcMessage("No volume change - cannot test vol/iso split"));
          Assert(std::abs(determinant(Kinematics::F_iso(F)) - 1.0) < tol,
                 ExcMessage("F_iso is not volume preserving"));
          Assert(std::abs(determinant(Kinematics::F_vol(F)) - determinant(F)) <
                   tol,
                 ExcMessage("F_vol has no dilatating action"));

          // Right Cauchy-Green tensor
          Assert((static_cast<Tensor<2, dim>>(Kinematics::C(F)) -
                  transpose(F) * F)
                     .norm() < tol,
                 ExcMessage("Incorrect computation of C"));

          // Left Cauchy-Green tensor
          Assert((static_cast<Tensor<2, dim>>(Kinematics::b(F)) -
                  F * transpose(F))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of b"));

          // Small strain tensor
          Assert((static_cast<Tensor<2, dim>>(Kinematics::epsilon(Grad_u)) -
                  0.5 * (Grad_u + transpose(Grad_u)))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of epsilon"));

          // Green-Lagrange strain tensor
          Assert((static_cast<Tensor<2, dim>>(Kinematics::E(F)) -
                  0.5 *
                    (Grad_u + transpose(Grad_u) + transpose(Grad_u) * Grad_u))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of E"));

          // Almansi strain tensor
          // Holzapfel 2.82
          Assert((static_cast<Tensor<2, dim>>(Kinematics::e(F)) -
                  transpose(invert(F)) *
                    static_cast<Tensor<2, dim>>(Kinematics::E(F)) * invert(F))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of e"));

          // --- Rate dependent ---

          // Material rates
          const Tensor<2, dim> &F_dot = qp_dot_Grad_u_t[q_point];

          // Material rate of deformation gradient tensor
          Assert((F_dot - (1.0 / delta_t) * (Grad_u - Grad_u_t1)).norm() < tol,
                 ExcMessage("Incorrect computation of F_dot"));

          // Spatial gradients
          const Tensor<2, dim> &dot_grad_u = qp_dot_grad_u_t[q_point];

          // Spatial velocity gradient
          Assert((static_cast<Tensor<2, dim>>(Kinematics::l(F, F_dot)) -
                  dot_grad_u)
                     .norm() < tol,
                 ExcMessage("Incorrect computation of l"));

          // Rate of deformation tensor
          Assert((static_cast<Tensor<2, dim>>(Kinematics::d(F, F_dot)) -
                  0.5 * (dot_grad_u + transpose(dot_grad_u)))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of d"));

          // Rate of spin tensor
          Assert((static_cast<Tensor<2, dim>>(Kinematics::w(F, F_dot)) -
                  0.5 * (dot_grad_u - transpose(dot_grad_u)))
                     .norm() < tol,
                 ExcMessage("Incorrect computation of w"));
        }
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_kinematic_tensors<2>();
  test_kinematic_tensors<3>();

  deallog << "OK" << std::endl;
}

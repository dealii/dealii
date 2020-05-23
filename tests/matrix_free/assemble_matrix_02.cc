// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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



// test FEEvaluation for assembling the Stokes matrix with Q2-Q1 finite elements

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"



const unsigned int degree_p = 1;

template <int dim>
void
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MappingQ<dim> mapping(degree_p + 2);

  // compute matrix with (\nabla v, \nabla u) + (v, 10 * u)
  {
    QGauss<dim> quadrature_formula(degree_p + 2);

    FEValues<dim> fe_values(mapping,
                            dof.get_fe(),
                            quadrature_formula,
                            update_values | update_gradients |
                              update_JxW_values);

    FEEvaluation<dim, degree_p + 1, degree_p + 2, dim> phi_u(
      mapping,
      dof.get_fe(),
      QGauss<1>(degree_p + 2),
      update_values | update_gradients | update_JxW_values,
      0);
    FEEvaluation<dim, degree_p, degree_p + 2> phi_p(dof.get_fe(), phi_u, dim);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> test_matrix(dofs_per_cell, dofs_per_cell);

    // The matrix part is from step-22
    const FEValuesExtractors::Vector     velocities(0);
    const FEValuesExtractors::Scalar     pressure(dim);
    std::vector<SymmetricTensor<2, dim>> phi_grads_u(dofs_per_cell);
    std::vector<double>                  div_phi_u(dofs_per_cell);
    std::vector<double>                  phi_pres(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                   endc = dof.end();
    for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        test_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                phi_grads_u[k] = fe_values[velocities].symmetric_gradient(k, q);
                div_phi_u[k]   = fe_values[velocities].divergence(k, q);
                phi_pres[k]    = fe_values[pressure].value(k, q);
              }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    cell_matrix(i, j) += (phi_grads_u[i] * phi_grads_u[j] -
                                          div_phi_u[i] * phi_pres[j] -
                                          phi_pres[i] * div_phi_u[j]) *
                                         fe_values.JxW(q);
                  }
              }
          }

        phi_u.reinit(cell);
        phi_p.reinit(cell);
        const unsigned int dofs_per_cell_u = phi_u.dofs_per_cell;
        const unsigned int dofs_per_cell_p = phi_p.dofs_per_cell;
        for (unsigned int i = 0; i < dofs_per_cell_u;
             i += VectorizedArray<double>::size())
          {
            const unsigned int n_items =
              i + VectorizedArray<double>::size() > dofs_per_cell_u ?
                (dofs_per_cell_u - i) :
                VectorizedArray<double>::size();
            for (unsigned int j = 0; j < dofs_per_cell_u; ++j)
              phi_u.begin_dof_values()[j] = VectorizedArray<double>();
            for (unsigned int v = 0; v < n_items; ++v)
              phi_u.begin_dof_values()[i + v][v] = 1.;

            phi_u.evaluate(false, true);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                VectorizedArray<double> div_u = phi_u.get_divergence(q);
                phi_u.submit_symmetric_gradient(phi_u.get_symmetric_gradient(q),
                                                q);
                phi_p.submit_value(-div_u, q);
              }
            phi_u.integrate(false, true);
            phi_p.integrate(true, false);

            for (unsigned int v = 0; v < n_items; ++v)
              {
                for (unsigned int j = 0; j < dofs_per_cell_u; ++j)
                  test_matrix(phi_u.get_internal_dof_numbering()[j],
                              phi_u.get_internal_dof_numbering()[i + v]) =
                    phi_u.begin_dof_values()[j][v];
                for (unsigned int j = 0; j < dofs_per_cell_p; ++j)
                  test_matrix(phi_p.get_internal_dof_numbering()[j],
                              phi_u.get_internal_dof_numbering()[i + v]) =
                    phi_p.begin_dof_values()[j][v];
              }
          }

        for (unsigned int i = 0; i < dofs_per_cell_p;
             i += VectorizedArray<double>::size())
          {
            const unsigned int n_items =
              i + VectorizedArray<double>::size() > dofs_per_cell_p ?
                (dofs_per_cell_p - i) :
                VectorizedArray<double>::size();
            for (unsigned int j = 0; j < dofs_per_cell_p; ++j)
              phi_p.begin_dof_values()[j] = VectorizedArray<double>();
            for (unsigned int v = 0; v < n_items; ++v)
              phi_p.begin_dof_values()[i + v][v] = 1.;

            phi_p.evaluate(true, false);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                phi_u.submit_divergence(-phi_p.get_value(q), q);
              }
            phi_u.integrate(false, true);

            for (unsigned int v = 0; v < n_items; ++v)
              for (unsigned int j = 0; j < dofs_per_cell_u; ++j)
                test_matrix(phi_u.get_internal_dof_numbering()[j],
                            phi_p.get_internal_dof_numbering()[i + v]) =
                  phi_u.begin_dof_values()[j][v];
          }
        test_matrix.add(-1., cell_matrix);
        deallog << test_matrix.frobenius_norm() << " ";
      }
    deallog << std::endl;
  }
}



template <int dim>
void
test()
{
  const SphericalManifold<dim> manifold;
  Triangulation<dim>           tria;
  GridGenerator::hyper_ball(tria);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (; cell != endc; ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  if (dim == 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FESystem<dim>   fe(FE_Q<dim>(degree_p + 1), dim, FE_Q<dim>(degree_p), 1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  do_test<dim>(dof);
}



int
main()
{
  initlog();

  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2>();
    deallog.pop();
    deallog.push("3d");
    test<3>();
    deallog.pop();
  }
}

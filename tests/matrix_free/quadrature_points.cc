// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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



// this function tests the correctness of cached quadrature points

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim, int fe_degree>
void
test()
{
  typedef double               number;
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
  tria.refine_global(5 - dim);

  MappingQ<dim>   mapping(4);
  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  deallog << "Testing " << fe.get_name() << std::endl;
  // std::cout << "Number of cells: " << tria.n_active_cells() << std::endl;
  // std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.mapping_update_flags  = update_quadrature_points;
    mf_data.reinit(mapping, dof, constraints, quad, data);
  }

  double                       error_points = 0, abs_points = 0;
  const unsigned int           n_cells = mf_data.n_macro_cells();
  FEEvaluation<dim, fe_degree> fe_eval(mf_data);
  FEValues<dim>                fe_values(mapping,
                          fe,
                          mf_data.get_quadrature(),
                          update_quadrature_points);

  typedef VectorizedArray<double> vector_t;
  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      fe_eval.reinit(cell);
      for (unsigned int j = 0; j < mf_data.n_components_filled(cell); ++j)
        {
          fe_values.reinit(mf_data.get_cell_iterator(cell, j));
          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            {
              abs_points += fe_values.quadrature_point(q).norm();
              for (unsigned int d = 0; d < dim; ++d)
                error_points += std::fabs(fe_values.quadrature_point(q)[d] -
                                          fe_eval.quadrature_point(q)[d][j]);
            }
        }
    }

  deallog << "Norm of difference: " << error_points / abs_points << std::endl
          << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 3>();
    deallog.pop();
  }
}

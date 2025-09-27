// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test FETools::get_projection_matrix for simplices.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/householder.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    if (component == 0)
      return p[0];
    else
      return 0.0;
  }
};

template <int dim, int spacedim = dim>
void
test()
{
  FE_SimplexP<dim, spacedim> fe_coarse(1);
  FE_SimplexP<dim, spacedim> fe_fine(2);
  MappingFE<dim>             mapping(FE_SimplexP<dim>(1));

  FullMatrix<double> matrix(fe_fine.n_dofs_per_cell(),
                            fe_coarse.n_dofs_per_cell());
  FETools::get_projection_matrix(fe_coarse, fe_fine, matrix);

  const unsigned int n_refinements = 2;

  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(
    tria, Utilities::pow(2, n_refinements));

  DoFHandler<dim, spacedim> dof_handler_coarse(tria);
  dof_handler_coarse.distribute_dofs(fe_coarse);

  DoFHandler<dim, spacedim> dof_handler_fine(tria);
  dof_handler_fine.distribute_dofs(fe_fine);

  Vector<double> vec_coarse(dof_handler_coarse.n_dofs());
  Vector<double> vec_fine(dof_handler_fine.n_dofs());

  Vector<double> temp_coarse(fe_coarse.n_dofs_per_cell());
  Vector<double> temp_fine(fe_fine.n_dofs_per_cell());

  VectorTools::interpolate(mapping,
                           dof_handler_coarse,
                           RightHandSideFunction<dim>(1),
                           vec_coarse);

  for (const auto &cell_coarse : dof_handler_coarse.active_cell_iterators())
    {
      cell_coarse->get_dof_values(vec_coarse, temp_coarse);

      matrix.vmult(temp_fine, temp_coarse);

      DoFCellAccessor<dim, spacedim, false>(&tria,
                                            cell_coarse->level(),
                                            cell_coarse->index(),
                                            &dof_handler_fine)
        .set_dof_values(temp_fine, vec_fine);
    }

  vec_fine.print(deallog.get_file_stream());

#if false
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler_coarse);
    data_out.add_data_vector(vec_coarse, "solution");

    data_out.build_patches(mapping, 2);

    std::ofstream output("test_coarse.vtk");
    data_out.write_vtk(output);
  }

  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler_fine);
    data_out.add_data_vector(vec_fine, "solution");

    data_out.build_patches(mapping, 2);

    std::ofstream output("test_fine.vtk");
    data_out.write_vtk(output);
  }
#endif
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
}

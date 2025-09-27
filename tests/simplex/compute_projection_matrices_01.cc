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


// Test FE_SimplexPoly::get_prolongation_matrix()
// (and indirectly FETools::compute_embedding_matrices() for simplices).


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

int
main()
{
  initlog();

  const int dim      = 2;
  const int spacedim = 2;

  FE_SimplexP<dim, spacedim> fe(2);
  MappingFE<dim>             mapping(FE_SimplexP<dim>(1));

  const unsigned int n_refinements = 2;

  Triangulation<dim, spacedim> tria_coarse, tria_fine;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria_coarse, 1);
  tria_coarse.refine_global(n_refinements);
  GridGenerator::subdivided_hyper_cube_with_simplices(tria_fine, 1);
  tria_fine.refine_global(n_refinements + 1);

  DoFHandler<dim, spacedim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(fe);

  DoFHandler<dim, spacedim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(fe);

  Vector<double> vec_coarse(dof_handler_coarse.n_dofs());
  Vector<double> vec_fine(dof_handler_fine.n_dofs());

  Vector<double> temp_coarse(fe.n_dofs_per_cell());
  Vector<double> temp_fine(fe.n_dofs_per_cell());

  // interpolate function onto coarse grid.
  VectorTools::interpolate(mapping,
                           dof_handler_coarse,
                           RightHandSideFunction<dim>(1),
                           vec_coarse);

  // project the result onto fine grid (cell by cell)
  for (const auto &cell_coarse : dof_handler_coarse.active_cell_iterators())
    {
      DoFCellAccessor<dim, spacedim, false> cell_coarse_on_fine_tria(
        &tria_fine,
        cell_coarse->level(),
        cell_coarse->index(),
        &dof_handler_fine);

      cell_coarse->get_dof_values(vec_coarse, temp_coarse);

      for (unsigned int c = 0; c < cell_coarse_on_fine_tria.n_children(); ++c)
        {
          fe.get_prolongation_matrix(c).vmult(temp_fine, temp_coarse);
          cell_coarse_on_fine_tria.child(c)->set_dof_values(temp_fine,
                                                            vec_fine);
        }
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

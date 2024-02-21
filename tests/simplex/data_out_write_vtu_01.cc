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



// Test DataOut::write_vtu() for simplex meshes.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

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
    return p[component % dim] * p[component % dim];
  }
};

template <int dim, int spacedim = dim>
void
test(const FiniteElement<dim, spacedim> &fe, const unsigned int n_components)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, dim == 2 ? 4 : 2);

  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());

  MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  VectorTools::interpolate(mapping,
                           dof_handler,
                           RightHandSideFunction<dim>(n_components),
                           solution);

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_compression;

  for (unsigned int n_subdivisions = 1; n_subdivisions <= 2; ++n_subdivisions)
    {
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");
      data_out.set_flags(vtk_flags);


      data_out.build_patches(mapping, n_subdivisions);

      data_out.write_vtu(deallog.get_file_stream());
    }
}

int
main()
{
  initlog();

  {
    const unsigned int dim = 2;
    test<dim>(FE_SimplexP<dim>(2) /*=degree*/, 1);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/), dim), dim);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/),
                            dim,
                            FE_SimplexP<dim>(1 /*=degree*/),
                            1),
              dim + 1);
  }
  {
    const unsigned int dim = 3;
    test<dim>(FE_SimplexP<dim>(2) /*=degree*/, 1);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/), dim), dim);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/),
                            dim,
                            FE_SimplexP<dim>(1 /*=degree*/),
                            1),
              dim + 1);
  }
}

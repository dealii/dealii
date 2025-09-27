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



// Test DataOut::write_vtk() for simplex meshes.

#include <deal.II/base/quadrature_lib.h>

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
test(const FiniteElement<dim, spacedim> &fe,
     const unsigned int                  n_components,
     const bool                          do_high_order)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, dim == 2 ? 4 : 2);

  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());

  MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  AffineConstraints<double> dummy;
  dummy.close();

  VectorTools::project(mapping,
                       dof_handler,
                       dummy,
                       QGaussSimplex<dim>(fe.tensor_degree() + 1),
                       RightHandSideFunction<dim>(n_components),
                       solution);

  static unsigned int counter = 0;

  for (unsigned int n_subdivisions = 1;
       n_subdivisions <= (do_high_order ? 3 : 2);
       ++n_subdivisions)
    {
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = do_high_order;

      DataOut<dim> data_out;
      data_out.set_flags(flags);

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");


      data_out.build_patches(mapping, n_subdivisions);

#if false
      std::ofstream output("test." + std::to_string(dim) + "." +
                           std::to_string(counter++) + ".vtk");
      data_out.write_vtk(output);
#else
      data_out.write_vtk(deallog.get_file_stream());
#endif
    }
}

int
main()
{
  initlog();

  for (unsigned int i = 0; i < 2; ++i)
    {
      const bool do_high_order = (i == 1);

      if (do_high_order)
        {
          const unsigned int dim = 2;
          test<dim>(FE_SimplexP<dim>(2) /*=degree*/, 1, do_high_order);
          test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/), dim),
                    dim,
                    do_high_order);
          test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/),
                                  dim,
                                  FE_SimplexP<dim>(1 /*=degree*/),
                                  1),
                    dim + 1,
                    do_high_order);
        }

      if (do_high_order ==
          false /*TODO: higher-order output not working for 3D*/)
        {
          const unsigned int dim = 3;
          test<dim>(FE_SimplexP<dim>(2) /*=degree*/, 1, do_high_order);
          test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/), dim),
                    dim,
                    do_high_order);
          test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/),
                                  dim,
                                  FE_SimplexP<dim>(1 /*=degree*/),
                                  1),
                    dim + 1,
                    do_high_order);
        }
    }
}

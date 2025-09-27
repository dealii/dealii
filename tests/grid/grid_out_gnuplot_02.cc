
// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Verify that the GNUPLOT output works with mappings that do not preserve
// vertex locations.

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int spacedim>
class Shift : public Function<spacedim>
{
public:
  Shift()
    : Function<spacedim>(spacedim)
  {}

  virtual double
  value(const Point<spacedim> &p,
        const unsigned int     component = 0) const override
  {
    return 2.0 * p[component];
  }
};

// overloads to get multiple grids for multiple dim and spacedim combinations
void
make_grid(Triangulation<2, 2> &triangulation)
{
  GridGenerator::hyper_shell(triangulation, Point<2>(), 2.0, 6.0, 12);
}

void
make_grid(Triangulation<2, 3> &triangulation)
{
  GridGenerator::hyper_sphere(triangulation, Point<3>(), 6.0);
  triangulation.refine_global(1); // need more cells
}

void
make_grid(Triangulation<3, 3> &triangulation)
{
  GridGenerator::hyper_shell(triangulation, Point<3>(), 2.0, 6.0, 12);
  triangulation.refine_global(0);
}

template <int dim, int spacedim = dim>
void
gnuplot_output(const GridOutFlags::Gnuplot &flags)
{
  deallog << "Triangulation<" << dim << ", " << spacedim << '>' << std::endl;

  Triangulation<dim, spacedim> triangulation;
  make_grid(triangulation);

  auto cell = triangulation.begin_active();
  cell->set_refine_flag(); // 0
  ++cell;
  ++cell;
  cell->set_refine_flag(); // 2
  ++cell;
  ++cell;
  cell->set_refine_flag(); // 4
  triangulation.execute_coarsening_and_refinement();

  FESystem<dim, spacedim>   displacement_fe(FE_Q<dim, spacedim>(1), spacedim);
  DoFHandler<dim, spacedim> displacement_dof_handler(triangulation);
  displacement_dof_handler.distribute_dofs(displacement_fe);

  Vector<double> displacements(displacement_dof_handler.n_dofs());
  VectorTools::interpolate(displacement_dof_handler,
                           Shift<spacedim>(),
                           displacements);
  MappingQ1Eulerian<dim, Vector<double>, spacedim> mapping(
    displacement_dof_handler, displacements);

  GridOut grid_out;
  grid_out.set_flags(flags);
  grid_out.write_gnuplot(triangulation, deallog.get_file_stream(), &mapping);
  deallog << std::endl << std::endl;
}

int
main()
{
  initlog();

  {
    GridOutFlags::Gnuplot flags;
    deallog << "curve lines and output extra boundary lines" << std::endl;
    flags.n_extra_curved_line_points      = 2;
    flags.curved_inner_cells              = true;
    flags.write_additional_boundary_lines = true;
    gnuplot_output<2>(flags);
    gnuplot_output<2, 3>(flags);
    // write one fewer line in 3D
    flags.n_extra_curved_line_points = 1;
    gnuplot_output<3>(flags);
  }

  deallog << "OK" << std::endl;
}

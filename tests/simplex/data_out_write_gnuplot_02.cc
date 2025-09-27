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



// Test DataOut::write_gnuplot() for a pyramid mesh

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
test(const FiniteElement<dim, spacedim> &fe)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::reference_cell(tria, ReferenceCells::Pyramid);

  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());

  MappingFE<dim> mapping(FE_PyramidP<dim>(1));

  AffineConstraints<double> dummy;
  dummy.close();

  VectorTools::project(mapping,
                       dof_handler,
                       dummy,
                       QGaussPyramid<dim>(fe.tensor_degree() + 1),
                       RightHandSideFunction<dim>(1),
                       solution);


  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");


  data_out.build_patches(mapping, 1);

#if false
  static unsigned int counter = 0;
      std::ofstream output("test." + std::to_string(dim) + "." +
                           std::to_string(counter++) + ".gnuplot");
      data_out.write_gnuplot(output);
#else
  data_out.write_gnuplot(deallog.get_file_stream());
#endif
}

int
main()
{
  initlog();

  const int dim = 3;
  test<dim>(FE_PyramidP<dim>(1));
}

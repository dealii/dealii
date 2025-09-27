// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This header contains a non-trivial testfunction F and a small test
// procedure shared among all generalized_interpolation_* tests.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim>
class F : public Function<dim>
{
public:
  F(const unsigned int n_comp, const unsigned int q)
    : Function<dim>(n_comp)
    , q(q)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double v = 0;
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int i = 0; i <= q; ++i)
        v += (d + 1) * (i + 1) * std::pow(p[d], 1. * i);
    return v;
  }

private:
  const unsigned int q;
};

template <int dim, typename T>
void
test(const FiniteElement<dim> &fe,
     const T                  &f,
     const unsigned int        order_mapping,
     bool                      distort_mesh,
     bool                      print_function_values = false)
{
  deallog << "dim " << dim << " " << fe.get_name() << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -0.3, 0.7);
  triangulation.refine_global(dim == 2 ? 2 : 1);
  if (distort_mesh)
    GridTools::distort_random(0.03, triangulation);

  MappingQ<dim> mapping(order_mapping);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, f, interpolant);

  // Print function value in origin:
  if (print_function_values)
    {
      Functions::FEFieldFunction<dim> f2(dof_handler, interpolant, mapping);
      deallog << "Function value at (0.0,0.0): ";
      for (unsigned int i = 0; i < fe.n_components(); ++i)
        deallog << f2.value(Point<dim>(), i) << " ";
      deallog << std::endl;
    }

  // Check that VectorTools::interpolate is in fact a
  // projection, i.e. applying the interpolation twice results in the same
  // vector:

  Functions::FEFieldFunction<dim> f2(dof_handler, interpolant, mapping);

  Vector<double> interpolant2(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, f2, interpolant2);

  interpolant2 -= interpolant;
  deallog << "Check projection property: " << interpolant2.linfty_norm()
          << std::endl;
}

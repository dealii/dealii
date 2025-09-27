// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check for matching values on cartesian mesh with UpdateFlags
// update_values and update_values|update_jacobian

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class F : public Function<dim>
{
public:
  F(const unsigned int q)
    : q(q)
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

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const unsigned int q_test = 2;

  MappingCartesian<dim> mapping;
  FE_Q<dim>             elem(q_test);
  // choose a point that is not right in the middle of the cell so that the
  // Jacobian contains many nonzero entries
  Point<dim> quad_p;
  for (int d = 0; d < dim; ++d)
    quad_p[d] = 0.42 + 0.11 * d;
  Quadrature<dim> quad(quad_p);

  Point<dim - 1> f_quad_p;
  for (int d = 0; d < dim - 1; ++d)
    f_quad_p[d] = 0.42 + 0.11 * d;
  Quadrature<dim - 1> f_quad(f_quad_p);


  {
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(elem);

    Vector<double> interpolant(dof_handler.n_dofs());
    VectorTools::interpolate(dof_handler, F<dim>(q_test), interpolant);

    const UpdateFlags    no_m(update_values);
    const UpdateFlags    w_m(update_values | update_jacobians);
    std::vector<double>  values(quad.size());
    std::vector<double>  values_m(quad.size());
    FEValues<dim>        fe_val(mapping, elem, quad, no_m);
    FEValues<dim>        fe_val_m(mapping, elem, quad, w_m);
    FEFaceValues<dim>    fe_f_val(mapping, elem, f_quad, no_m);
    FEFaceValues<dim>    fe_f_val_m(mapping, elem, f_quad, w_m);
    FESubfaceValues<dim> fe_subf_val(mapping, elem, f_quad, no_m);
    FESubfaceValues<dim> fe_subf_val_m(mapping, elem, f_quad, w_m);

    deallog << dim << " Checking no mapping FEValues behavior: " << std::endl;
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);
        fe_val_m.reinit(cell);

        fe_val.get_function_values(interpolant, values);
        fe_val_m.get_function_values(interpolant, values_m);
        Assert(values[0] == values_m[0], ExcInternalError());

        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          {
            fe_f_val.reinit(cell, f);
            fe_f_val_m.reinit(cell, f);

            fe_f_val.get_function_values(interpolant, values);
            fe_f_val_m.get_function_values(interpolant, values_m);
            Assert(values[0] == values_m[0], ExcInternalError());

            // Also check the Jacobian with FESubfaceValues
            if (cell->at_boundary(f) == false &&
                cell->neighbor(f)->level() < cell->level())
              {
                fe_subf_val.reinit(
                  cell->neighbor(f),
                  cell->neighbor_face_no(f),
                  cell->neighbor_of_coarser_neighbor(f).second);
                fe_subf_val_m.reinit(
                  cell->neighbor(f),
                  cell->neighbor_face_no(f),
                  cell->neighbor_of_coarser_neighbor(f).second);

                fe_subf_val.get_function_values(interpolant, values);
                fe_subf_val_m.get_function_values(interpolant, values_m);
                Assert(values[0] == values_m[0], ExcInternalError());
              }
          }
      }
    deallog << "OK" << std::endl;
  }
}


int
main()
{
  initlog();
  deallog << std::setprecision(8) << std::fixed;

  test<2>();
  test<3>();

  return 0;
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Like _03 but with a graph coloring algorithm for the iterators. The
// graph coloring here is simple since all writes conflict (we just
// add to a scalar)

#include <deal.II/base/parallel.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"


template <int dim>
double
value(const Point<dim> &p)
{
  double val = 0;
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int i = 0; i <= 1; ++i)
      val += std::pow(p[d], 1. * i);
  return val;
}

namespace
{
  template <int dim>
  struct Scratch
  {
    Scratch(const FiniteElement<dim> &fe, const Quadrature<dim> &quadrature)
      : fe_collection(fe)
      , quadrature_collection(quadrature)
      , x_fe_values(fe_collection,
                    quadrature_collection,
                    update_quadrature_points)
      , rhs_values(quadrature_collection.size())
    {}

    Scratch(const Scratch &data)
      : fe_collection(data.fe_collection)
      , quadrature_collection(data.quadrature_collection)
      , x_fe_values(fe_collection,
                    quadrature_collection,
                    update_quadrature_points)
      , rhs_values(data.rhs_values)
    {}

    const FiniteElement<dim> &fe_collection;
    const Quadrature<dim>    &quadrature_collection;

    FEValues<dim> x_fe_values;

    std::vector<double> rhs_values;
  };

  struct CopyData
  {
    std::vector<double> cell_rhs;
  };
} // namespace

void
zero_subrange(const unsigned int   begin,
              const unsigned int   end,
              std::vector<double> &dst)
{
  for (unsigned int i = begin; i < end; ++i)
    dst[i] = 0;
}


void
zero_element(std::vector<double> &dst, const unsigned int i)
{
  dst[i] = 0;
}


template <int dim>
void
mass_assembler(const typename Triangulation<dim>::active_cell_iterator &cell,
               Scratch<dim>                                            &data,
               CopyData &copy_data)
{
  data.x_fe_values.reinit(cell);

  const Point<dim> q = data.x_fe_values.quadrature_point(0);

  // this appears to be the key: the following two ways both overwrite some
  // of the memory in which we store the quadrature point location.
  parallel::apply_to_subranges(0U,
                               copy_data.cell_rhs.size(),
                               std::bind(&zero_subrange,
                                         std::placeholders::_1,
                                         std::placeholders::_2,
                                         std::ref(copy_data.cell_rhs)),
                               1);

  AssertThrow(q == data.x_fe_values.quadrature_point(0), ExcInternalError());

  copy_data.cell_rhs[0] = value(data.x_fe_values.quadrature_point(0));
}


void
copy_local_to_global(const CopyData &data, double *sum)
{
  *sum += data.cell_rhs[0];
}


// the function that defines conclicts. we always write into the same
// field, so we need to always return the same index
template <int dim>
std::vector<types::global_dof_index>
conflictor(const typename Triangulation<dim>::active_cell_iterator &)
{
  return std::vector<types::global_dof_index>(1, types::global_dof_index());
}


void
do_project()
{
  static const int dim = 3;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  FE_Nothing<dim> fe;
  QMidpoint<dim>  q;

  for (unsigned int i = 0; i < 12; ++i)
    {
      std::vector<double> tmp;

      double       sum = 0;
      Scratch<dim> assembler_data(fe, q);
      CopyData     copy_data;
      copy_data.cell_rhs.resize(8);
      WorkStream::run(GraphColoring::make_graph_coloring(
                        triangulation.begin_active(),
                        triangulation.end(),
                        std::function<std::vector<types::global_dof_index>(
                          const Triangulation<dim>::active_cell_iterator &)>(
                          &conflictor<dim>)),
                      &mass_assembler<dim>,
                      std::bind(&copy_local_to_global,
                                std::placeholders::_1,
                                &sum),
                      assembler_data,
                      copy_data,
                      8,
                      1);

      Assert(std::fabs(sum - 288.) < 1e-12, ExcInternalError());
      deallog << sum << std::endl;
    }
}


int
main()
{
  initlog();

  do_project();
}

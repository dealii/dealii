// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check InverseQuadraticApproximation used for the initial guess in
// MappingQ::transform_points_real_to_unit_cell

#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_internal.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim, int spacedim>
void
print_result(const unsigned int                  mapping_degree,
             const Triangulation<dim, spacedim> &tria,
             const Point<dim>                    p)
{
  deallog << "Testing " << dim << "D with point " << p << std::endl;

  FE_Q<dim>     dummy(mapping_degree);
  MappingQ<dim> mapping(mapping_degree);

  FEValues<dim> fe_values(mapping,
                          dummy,
                          Quadrature<dim>(dummy.get_unit_support_points()),
                          update_quadrature_points);

  std::vector<Point<1>> mapping_points(
    QGaussLobatto<1>(mapping_degree + 1).get_points());
  std::vector<unsigned int> renumber =
    FETools::lexicographic_to_hierarchic_numbering<dim>(mapping_degree);
  std::vector<Point<dim>> mapping_unit_support_points =
    internal::MappingQImplementation::unit_support_points<dim>(mapping_points,
                                                               renumber);

  for (const auto &cell : tria.active_cell_iterators())
    {
      fe_values.reinit(cell);
      try
        {
          const Point<dim> p_unit =
            mapping.transform_real_to_unit_cell(cell, p);
          deallog << "Testing on cell " << cell->id() << " with center "
                  << cell->center(true) << std::endl;
          deallog << "Exact inverse:                   " << p_unit << std::endl;
          deallog << "Affine approximation:            "
                  << cell->real_to_unit_cell_affine_approximation(p)
                  << std::endl;
          internal::MappingQImplementation::
            InverseQuadraticApproximation<dim, spacedim>
              approx(fe_values.get_quadrature_points(),
                     mapping_unit_support_points);
          deallog << "Inverse quadratic approximation: " << approx.compute(p)
                  << std::endl
                  << std::endl;
        }
      catch (typename Mapping<dim, spacedim>::ExcTransformationFailed &)
        {}
    }
  deallog << std::endl;
}



template <int dim, int spacedim>
void
test(const unsigned mapping_degree, const unsigned int n_ref)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.8, 1., dim == 2 ? 3 : 6);
  tria.refine_global(n_ref);
  {
    const double phi   = 0.56 * numbers::PI;
    const double theta = 0.49 * numbers::PI;
    for (unsigned int i = 0; i <= 4; ++i)
      {
        const double r = 0.8 + 0.2 * i / 4;
        Point<dim>   p;
        p[0] = std::cos(phi) * std::sin(theta) * r;
        p[1] = std::sin(phi) * std::sin(theta) * r;
        if (dim == 3)
          p[2] = std::cos(theta) * r;
        print_result(mapping_degree, tria, p);
      }
  }
  if (dim == 3)
    {
      const double phi   = 0.3 * numbers::PI;
      const double theta = 0.1;
      for (unsigned int i = 0; i <= 6; ++i)
        {
          const double r = 0.8 + 0.2 * i / 6;
          Point<dim>   p;
          p[0] = std::cos(phi) * std::sin(theta) * r;
          p[1] = std::sin(phi) * std::sin(theta) * r;
          if (dim == 3)
            p[2] = std::cos(theta) * r;
          print_result(mapping_degree, tria, p);
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<2, 2>(6, 1);
  test<3, 3>(4, 0);
}

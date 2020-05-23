// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// find_all_active_cells_around_point for a deformed shell mesh

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim, int spacedim>
void
print_result(const Mapping<dim, spacedim> &      mapping,
             const Triangulation<dim, spacedim> &tria,
             const Point<dim>                    p,
             const double                        tolerance)
{
  deallog << "Testing " << dim << "D with point " << p << " tolerance "
          << tolerance << std::endl;
  auto c_p =
    GridTools::find_all_active_cells_around_point(mapping, tria, p, tolerance);
  for (auto i : c_p)
    deallog << "Cell: " << i.first->id() << " unit point " << i.second
            << std::endl;
  deallog << std::endl;
}



template <int dim>
void
do_test(const Mapping<dim> &mapping, const Triangulation<dim> &tria)
{
  {
    const double phi   = 0.;
    const double theta = 0.5 * numbers::PI;
    for (unsigned int i = 0; i <= 4; ++i)
      {
        const double r = 0.8 + 0.2 * i / 4;
        Point<dim>   p;
        p[0] = std::cos(phi) * std::sin(theta) * r;
        p[1] = std::sin(phi) * std::sin(theta) * r;
        if (dim == 3)
          p[2] = std::cos(theta) * r;
        print_result(mapping, tria, p, 1e-8);
      }
  }
  {
    const double phi   = 0.25 * numbers::PI;
    const double theta = 0.5 * numbers::PI;
    for (unsigned int i = 0; i <= 4; ++i)
      {
        const double r = 0.8 + 0.2 * i / 4;
        Point<dim>   p;
        p[0] = std::cos(phi) * std::sin(theta) * r;
        p[1] = std::sin(phi) * std::sin(theta) * r;
        if (dim == 3)
          p[2] = std::cos(theta) * r;
        print_result(mapping, tria, p, 1e-8);
        print_result(mapping, tria, p, 1e-12);
      }
  }
  {
    const double phi   = 0.3 * numbers::PI;
    const double theta = 0.5 * numbers::PI;
    for (unsigned int i = 0; i <= 4; ++i)
      {
        const double r = 0.8 + 0.2 * i / 4;
        Point<dim>   p;
        p[0] = std::cos(phi) * std::sin(theta) * r;
        p[1] = std::sin(phi) * std::sin(theta) * r;
        if (dim == 3)
          p[2] = std::cos(theta) * r;
        print_result(mapping, tria, p, 1e-8);
      }
  }
  {
    const double phi   = 0.5 * numbers::PI;
    const double theta = 0.5 * numbers::PI;
    for (unsigned int i = 0; i <= 4; ++i)
      {
        const double r = 0.8 + 0.2 * i / 4;
        Point<dim>   p;
        p[0] = std::cos(phi) * std::sin(theta) * r;
        p[1] = std::sin(phi) * std::sin(theta) * r;
        if (dim == 3)
          p[2] = std::cos(theta) * r;
        print_result(mapping, tria, p, 1e-8);
      }
  }
  {
    const double phi   = 0.75 * numbers::PI;
    const double theta = 0.5 * numbers::PI;
    for (unsigned int i = 0; i <= 4; ++i)
      {
        const double r = 0.8 + 0.2 * i / 4;
        Point<dim>   p;
        p[0] = std::cos(phi) * std::sin(theta) * r;
        p[1] = std::sin(phi) * std::sin(theta) * r;
        if (dim == 3)
          p[2] = std::cos(theta) * r;
        print_result(mapping, tria, p, 1e-8);
      }
  }
  if (dim == 3)
    {
      const double phi   = 0.75 * numbers::PI;
      const double theta = 0.;
      for (unsigned int i = 0; i <= 4; ++i)
        {
          const double r = 0.8 + 0.2 * i / 4;
          Point<dim>   p;
          p[0] = std::cos(phi) * std::sin(theta) * r;
          p[1] = std::sin(phi) * std::sin(theta) * r;
          if (dim == 3)
            p[2] = std::cos(theta) * r;
          print_result(mapping, tria, p, 1e-8);
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
          print_result(mapping, tria, p, 1e-8);
        }
    }
}



template <int dim, int spacedim>
void
test(unsigned int n_ref)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.8, 1., dim == 2 ? 3 : 6);
  tria.refine_global(n_ref);
  const unsigned int             fe_degree = dim == 2 ? 8 : 4;
  MappingQGeneric<dim, spacedim> mapping_q(fe_degree);

  FE_Q<dim>               fe_q(fe_degree);
  FESystem<dim, dim>      fe_system(fe_q, dim);
  dealii::DoFHandler<dim> dofh;
  dofh.initialize(tria, fe_system);
  dofh.distribute_dofs(fe_system);

  const ComponentMask mask(dim, true);
  Vector<double>      nodes;
  nodes.reinit(dofh.n_dofs());

  VectorTools::get_position_vector(dofh, nodes, mask);
  MappingFEField<dim, dim, Vector<double>, DoFHandler<dim>> mapping(dofh,
                                                                    nodes,
                                                                    mask);

  deallog << "Test with MappingQGeneric in " << dim << "D on "
          << tria.n_active_cells() << " cells:" << std::endl;
  do_test(mapping_q, tria);
  deallog << std::endl;
  deallog << "Test with MappingFEField in " << dim << "D on "
          << tria.n_active_cells() << " cells:" << std::endl;
  do_test(mapping, tria);
  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<2, 2>(1);
  test<2, 2>(2);
  test<3, 3>(1);
  test<3, 3>(2);
}

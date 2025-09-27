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


// Test hp::FEFaceValues for different face quadrature rules.


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class Fu : public Function<dim>
{
public:
  Fu(const unsigned int component = numbers::invalid_unsigned_int)
    : Function<dim>(component == numbers::invalid_unsigned_int ? dim : 1)
    , component(component)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int i = 0) const
  {
    if (component == numbers::invalid_unsigned_int)
      return p[i];
    else
      return p[component];
  }

  const unsigned int component;
};

template <int dim>
void
test()
{}

template <>
void
test<2>()
{
  const unsigned int dim = 2;

  // test FEFaceValues for FE_System(FE_Q)
  {
    const std::vector<hp::QCollection<dim - 1>> quad_ref = {
      hp::QCollection<dim - 1>(QGauss<dim - 1>(1),
                               QGauss<dim - 1>(2),
                               QGauss<dim - 1>(3),
                               QGauss<dim - 1>(4)),
      hp::QCollection<dim - 1>(QGauss<dim - 1>(5),
                               QGauss<dim - 1>(6),
                               QGauss<dim - 1>(7),
                               QGauss<dim - 1>(8))};

    hp::MappingCollection<dim> mapping(MappingFE<dim>(FE_Q<dim>(1)));
    hp::FECollection<dim>      fe(FESystem<dim>(FE_Q<dim>{3}, dim));

    const UpdateFlags flags =
      update_values | update_quadrature_points | update_JxW_values;


    hp::FEFaceValues<dim> fe_face_values(mapping, fe, quad_ref, flags);

    Triangulation<dim> tria;
    GridGenerator::hyper_cube(tria);

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    Vector<double> vector_0(dof_handler.n_dofs());

    VectorTools::interpolate(mapping[0 /*TODO*/],
                             dof_handler,
                             Fu<dim>(),
                             vector_0);

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const auto face_no : cell->face_indices())
        for (unsigned int q_index = 0; q_index < quad_ref.size(); ++q_index)
          {
            fe_face_values.reinit(cell, face_no, q_index);

            std::vector<Vector<double>> values_0(
              fe_face_values.get_present_fe_values().n_quadrature_points,
              Vector<double>(dim));

            fe_face_values.get_present_fe_values().get_function_values(
              vector_0, values_0);

            deallog << "q_index=" << q_index << " face_no=" << face_no << ':'
                    << std::endl;

            for (unsigned int q = 0; q < values_0.size(); ++q)
              deallog << values_0[q][0] << ' ' << values_0[q][1] << ' '
                      << std::endl;

            deallog << std::endl;
          }
  }
}

int
main()
{
  initlog();

  test<2>();
}

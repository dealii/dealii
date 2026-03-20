// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test MappingQ::fill_fe_values when the associated FEValues object is kept
// along for multiple triangulation objects

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
do_test(const unsigned int mapping_degree)
{
  FE_Nothing<dim> fe;
  MappingQ<dim>   mapping(mapping_degree);
  deallog << "dim = " << dim << ", mapping degree " << mapping_degree
          << std::endl;

  FEValues<dim>     fe_values(mapping, fe, QGauss<dim>(1), update_jacobians);
  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   QGauss<dim - 1>(1),
                                   update_jacobians);

  for (unsigned int i = 0; i < 3; ++i)
    {
      const double edge_length = Utilities::pow(0.5, i);
      deallog << "Testing triangulation create of edge length " << edge_length
              << std::endl;
      Triangulation<dim> tria;
      GridGenerator::hyper_cube(tria, 0, edge_length);

      fe_values.reinit(tria.begin_active());
      deallog << "Jacobian cell: " << fe_values.jacobian(0) << std::endl;

      fe_face_values.reinit(tria.begin_active(), 0);
      deallog << "Jacobian face: " << fe_face_values.jacobian(0) << std::endl;
    }

  for (unsigned int i = 0; i < 3; ++i)
    {
      const double edge_length = Utilities::pow(0.5, i);
      deallog << "Testing triangulation via pointer of edge length "
              << edge_length << std::endl;
      std::unique_ptr<Triangulation<dim>> tria =
        std::make_unique<Triangulation<dim>>();
      GridGenerator::hyper_cube(*tria, 0, edge_length);

      fe_values.reinit(tria->begin_active());
      deallog << "Jacobian cell: " << fe_values.jacobian(0) << std::endl;

      fe_face_values.reinit(tria->begin_active(), 0);
      deallog << "Jacobian face: " << fe_face_values.jacobian(0) << std::endl;
    }

  deallog << "Testing two triangulations with edge length 0.5 and 0.25"
          << std::endl;
  Triangulation<dim> tria1;
  GridGenerator::hyper_cube(tria1, 0, 0.5);

  fe_values.reinit(tria1.begin_active());
  deallog << "Jacobian cell: " << fe_values.jacobian(0) << std::endl;

  fe_face_values.reinit(tria1.begin_active(), 0);
  deallog << "Jacobian face: " << fe_face_values.jacobian(0) << std::endl;
  Triangulation<dim> tria2;
  GridGenerator::hyper_cube(tria2, 0, 0.25);

  fe_values.reinit(tria2.begin_active());
  deallog << "Jacobian cell: " << fe_values.jacobian(0) << std::endl;

  fe_face_values.reinit(tria2.begin_active(), 0);
  deallog << "Jacobian face: " << fe_face_values.jacobian(0) << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();

  do_test<1>(1);
  do_test<1>(2);
  do_test<2>(1);
  do_test<2>(2);
  do_test<3>(1);
  do_test<3>(2);
}

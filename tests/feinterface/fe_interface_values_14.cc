// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Set up a two cells with different elements, such that neither element
// dominates. Create an FEInterfaceValues object and call reinit on it.
//
// In this test case, there is only one quadrature and mapping, so
// that is the only reasonable choice.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
make_2_cells(Triangulation<dim> &tria);

template <>
void
make_2_cells<2>(Triangulation<2> &tria)
{
  const unsigned int        dim         = 2;
  std::vector<unsigned int> repetitions = {2, 1};
  Point<dim>                p1;
  Point<dim>                p2(2.0, 1.0);

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
}


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  make_2_cells(tria);

  const unsigned      max_fe_degree = 2;
  const FESystem<dim> fe_left(FE_Q<dim>(1), FE_Q<dim>(2));
  const FESystem<dim> fe_right(FE_Q<dim>(2), FE_Q<dim>(1));

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe_left);
  fe_collection.push_back(fe_right);

  DoFHandler<dim> dofh(tria);

  const unsigned int face_index = 1;

  auto cell     = dofh.begin();
  auto neighbor = cell->neighbor(face_index);

  cell->set_active_fe_index(0);
  neighbor->set_active_fe_index(1);

  dofh.distribute_dofs(fe_collection);

  const QGauss<dim - 1>          quadrature(max_fe_degree + 1);
  const hp::QCollection<dim - 1> q_collection(quadrature);

  const MappingQ<dim>              mapping(1);
  const hp::MappingCollection<dim> mapping_collection(mapping);


  deallog << "Adjacent finite elements: " << cell->get_fe().get_name()
          << " and " << neighbor->get_fe().get_name() << std::endl;

  const UpdateFlags      update_flags = update_values;
  FEInterfaceValues<dim> fiv(mapping_collection,
                             fe_collection,
                             q_collection,
                             update_flags);
  fiv.reinit(cell,
             face_index,
             /* subface_no */ numbers::invalid_unsigned_int,
             neighbor,
             cell->neighbor_of_neighbor(face_index),
             /* subface_no_neighbor */ numbers::invalid_unsigned_int,
             /* q_index */ numbers::invalid_unsigned_int,
             /* mapping_index */ numbers::invalid_unsigned_int,
             /* fe_index */ numbers::invalid_unsigned_int);

  deallog << "Chosen quadrature for the common face has "
          << fiv.get_fe_face_values(0).get_quadrature().size()
          << " quadrature points. These points are located at:" << std::endl;
  for (const auto p : fiv.get_fe_face_values(0).get_quadrature().get_points())
    deallog << "  " << p << std::endl;

  // Also check that we see the same quadrature formula from both
  // sides:
  Assert(fiv.get_fe_face_values(0).get_quadrature().get_points() ==
           fiv.get_fe_face_values(1).get_quadrature().get_points(),
         ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
}

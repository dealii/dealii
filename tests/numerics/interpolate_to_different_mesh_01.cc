// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// test for interpolate_to_different_mesh with dealii::DoFHandler

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace std;

template <int spacedim>
class F1 : public Function<spacedim>
{
public:
  F1()
    : Function<spacedim>()
  {}
  virtual double
  value(const Point<spacedim> &p, const unsigned int component = 0) const
  {
    return p(0) * p(0);
  }
};

template <unsigned int spacedim>
void
check(const unsigned int refinement_1, const unsigned int refinement_2)
{
  MappingQGeneric<spacedim> mapping(1);

  Triangulation<spacedim> tria_1, tria_2;
  GridGenerator::hyper_cube(tria_1);
  GridGenerator::hyper_cube(tria_2);

  tria_1.refine_global(refinement_1);
  tria_2.refine_global(refinement_2);

  DoFHandler<spacedim> dof_handler_1(tria_1);
  DoFHandler<spacedim> dof_handler_2(tria_2);

  FE_Q<spacedim> fe(1);
  dof_handler_1.distribute_dofs(fe);
  dof_handler_2.distribute_dofs(fe);

  Vector<double> u_1(dof_handler_1.n_dofs()), u_2(dof_handler_2.n_dofs());
  F1<spacedim>   f_1;
  VectorTools::interpolate(mapping, dof_handler_1, f_1, u_1);
  VectorTools::interpolate_to_different_mesh(dof_handler_1,
                                             u_1,
                                             dof_handler_2,
                                             u_2);

  Point<spacedim>                 support_point, unit_support_point;
  vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());
  deallog << "dof values of dof_handler_1:" << endl;
  for (const auto &cell : dof_handler_1.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int shapefun = 0;
           shapefun < cell->get_fe().n_dofs_per_cell();
           ++shapefun)
        {
          unit_support_point = cell->get_fe().unit_support_point(shapefun);
          support_point =
            mapping.transform_unit_to_real_cell(cell, unit_support_point);
          deallog << " " << support_point << ":\n "
                  << u_1[local_dof_indices[shapefun]] << endl
                  << endl;
        }
    }
  deallog << endl << "dof values of dof_handler_2:" << endl;
  for (const auto &cell : dof_handler_2.active_cell_iterators())

    {
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int shapefun = 0;
           shapefun < cell->get_fe().n_dofs_per_cell();
           ++shapefun)
        {
          unit_support_point = cell->get_fe().unit_support_point(shapefun);
          support_point =
            mapping.transform_unit_to_real_cell(cell, unit_support_point);
          deallog << " " << support_point << ":\n "
                  << u_2[local_dof_indices[shapefun]] << endl
                  << endl;
        }
    }
}

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog
    << "### 2D-Case, first cell refined once, second cell unrefined###\n\n";
  check<2>(1, 0);
  deallog << endl;

  deallog
    << "### 2D-Case, first cell unrefined, second cell refined once###\n\n";
  check<2>(0, 1);
  deallog << endl;

  deallog
    << "### 2D-Case, first cell refined once, second cell refined once###\n\n";
  check<2>(1, 1);
  deallog << endl;

  deallog
    << "### 3D-Case, first cell refined once, second cell unrefined###\n\n";
  check<3>(1, 0);
  deallog << endl;

  deallog
    << "### 3D-Case, first cell unrefined, second cell refined once###\n\n";
  check<3>(0, 1);
  deallog << endl;

  deallog
    << "### 3D-Case, first cell refined once, second cell refined once###\n\n";
  check<3>(1, 1);
  deallog << endl;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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


// Test basic properties of FEInterfaceValues for multilevel

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim>
void
inspect_fiv(FEInterfaceValues<dim> &fiv)
{
  deallog << "at_boundary(): " << fiv.at_boundary() << "\n"
          << "n_current_interface_dofs(): " << fiv.n_current_interface_dofs()
          << "\n";

  std::vector<types::global_dof_index> indices =
    fiv.get_interface_dof_indices();
  Assert(indices.size() == fiv.n_current_interface_dofs(), ExcInternalError());

  deallog << "interface_dof_indices: ";
  for (auto i : indices)
    deallog << i << ' ';
  deallog << "\n";


  unsigned int idx = 0;
  for (auto v : indices)
    {
      deallog << "  index " << idx << " global_dof_index:" << v << ":\n";

      const auto pair = fiv.interface_dof_to_dof_indices(idx);
      deallog << "    dof indices: " << static_cast<int>(pair[0]) << " | "
              << static_cast<int>(pair[1]) << "\n";

      ++idx;
    }

  deallog << "\n";
  deallog << std::endl;
}

void
make_2_cells(Triangulation<2> &tria)
{
  const unsigned int        dim         = 2;
  std::vector<unsigned int> repetitions = {2, 1};
  Point<dim>                p1;
  Point<dim>                p2(2.0, 1.0);

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
}

void
make_2_cells(Triangulation<3> &tria)
{
  const unsigned int        dim         = 3;
  std::vector<unsigned int> repetitions = {2, 1, 1};
  Point<dim>                p1;
  Point<dim>                p2(2.0, 1.0, 1.0);

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
}


template <int dim>
void
test()
{
  Triangulation<dim> tria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  make_2_cells(tria);
  tria.refine_global();


  DoFHandler<dim> dofh(tria);
  FE_DGQ<dim>     fe(1);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();

  deallog << "n_dofs = " << dofh.n_dofs() << " n_dofs(0)= " << dofh.n_dofs(0)
          << " n_dofs(1)= " << dofh.n_dofs(1) << std::endl;

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;

  FEInterfaceValues<dim> fiv(mapping,
                             fe,
                             QGauss<dim - 1>(fe.degree + 1),
                             update_flags);


  auto cell = dofh.begin_mg(0);

  deallog << "** interface between cell 0 and 1 on level 0: **\n";

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    if (!cell->at_boundary(f))
      {
        fiv.reinit(cell,
                   f,
                   numbers::invalid_unsigned_int,
                   cell->neighbor(f),
                   cell->neighbor_of_neighbor(f),
                   numbers::invalid_unsigned_int);

        Assert(fiv.get_fe_face_values(0).get_cell() == cell,
               ExcInternalError());
        Assert(fiv.get_fe_face_values(1).get_cell() == cell->neighbor(f),
               ExcInternalError());
        Assert(fiv.n_current_interface_dofs() == 2 * fe.n_dofs_per_cell(),
               ExcInternalError());
        Assert(!fiv.at_boundary(), ExcInternalError());

        auto mycell = cell;
        for (unsigned int c = 0; c < 2; ++c)
          {
            std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
            mycell->get_mg_dof_indices(indices);
            deallog << "cell " << c << ": ";
            for (auto i : indices)
              deallog << i << ' ';
            deallog << "\n";
            ++mycell;
          }

        inspect_fiv(fiv);
      }


  deallog << "** interface on level 1: **\n";
  cell = dofh.begin_mg(1);

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    if (!cell->at_boundary(f))
      {
        fiv.reinit(cell,
                   f,
                   numbers::invalid_unsigned_int,
                   cell->neighbor(f),
                   cell->neighbor_of_neighbor(f),
                   numbers::invalid_unsigned_int);

        Assert(fiv.get_fe_face_values(0).get_cell() == cell,
               ExcInternalError());
        Assert(fiv.get_fe_face_values(1).get_cell() == cell->neighbor(f),
               ExcInternalError());
        Assert(fiv.n_current_interface_dofs() == 2 * fe.n_dofs_per_cell(),
               ExcInternalError());
        Assert(!fiv.at_boundary(), ExcInternalError());

        auto mycell = cell;
        for (unsigned int c = 0; c < 2; ++c)
          {
            std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
            mycell->get_mg_dof_indices(indices);
            deallog << "cell " << c << ": ";
            for (auto i : indices)
              deallog << i << ' ';
            deallog << "\n";
            mycell = cell->neighbor(f);
          }

        inspect_fiv(fiv);
        break; // only look at the first internal face
      }

  deallog << "** boundary interface on cell 1 **\n";
  cell = dofh.begin_mg(0);
  {
    ++cell;
    fiv.reinit(cell, 1);
    Assert(fiv.get_fe_face_values(0).get_cell() == cell, ExcInternalError());
    Assert(fiv.n_current_interface_dofs() == fe.n_dofs_per_cell(),
           ExcInternalError());
    Assert(fiv.at_boundary(), ExcInternalError());
    inspect_fiv(fiv);
  }
}



int
main()
{
  initlog();
  test<2>();
}

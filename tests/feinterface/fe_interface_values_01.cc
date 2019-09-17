// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


// Test basic properties of FEInterfaceValues, global refinement, no-hp

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
    deallog << i << " ";
  deallog << "\n";


  unsigned int idx = 0;
  for (auto v : indices)
    {
      deallog << "  index " << idx << " global_dof_index:" << v << ":\n";

      const auto pair = fiv.interface_dof_to_dof_indices(idx);
      deallog << "    dof indices: " << pair[0] << " | " << pair[1] << "\n";

      ++idx;
    }

  deallog << "\n";
  deallog << std::endl;
}


template <int dim>
void
make_2_cells(Triangulation<dim> &tria);

template <>
void make_2_cells<2>(Triangulation<2> &tria)
{
  const unsigned int        dim         = 2;
  std::vector<unsigned int> repetitions = {2, 1};
  Point<dim>                p1;
  Point<dim>                p2(2.0, 1.0);

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
}

template <>
void make_2_cells<3>(Triangulation<3> &tria)
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
  Triangulation<dim> tria;
  make_2_cells(tria);

  DoFHandler<dim> dofh(tria);
  FE_DGQ<dim>     fe(1);
  dofh.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;

  FEInterfaceValues<dim> fiv(mapping,
                             fe,
                             QGauss<dim - 1>(fe.degree + 1),
                             update_flags);


  auto cell = dofh.begin();

  deallog << "** interface between cell 0 and 1 **\n";

  for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
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
            mycell->get_dof_indices(indices);
            deallog << "cell " << c << ": ";
            for (auto i : indices)
              deallog << i << " ";
            deallog << "\n";
            ++mycell;
          }

        inspect_fiv(fiv);
      }

  deallog << "** boundary interface on cell 1 **\n";

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

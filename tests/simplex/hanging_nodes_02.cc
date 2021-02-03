// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// verify hanging node constraints on locally p-refined simplex mesh
//
// dofs will be enumerated as follows
//  scenario 1:    scenario 2:
//   6-------4      2---4---3
//   |\      |      |\      |
//   |  \    |      |  \    |
//   3   2   |      |   5   6
//   |    \  |      |    \  |
//   |      \|      |      \|
//   0---1---5      0-------1


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/grid_generator.h>

#include "../tests.h"


// ----- diagnostics -----

template <int dim>
void
print_dof_indices_on_faces(const DoFHandler<dim> &dofh)
{
  std::vector<types::global_dof_index> dof_indices;

  for (const auto &cell : dofh.active_cell_iterators())
    for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        const auto &face = cell->face(f);

        Assert(!face->has_children(), ExcInternalError());

        const unsigned int fe_index = cell->active_fe_index();
        const auto &       fe       = cell->get_fe();

        dof_indices.resize(fe.n_dofs_per_face(f));
        face->get_dof_indices(dof_indices, fe_index);

        deallog << "cell:" << cell->active_cell_index() << " face:" << f
                << " dofs:";
        for (const auto &i : dof_indices)
          deallog << i << " ";
        deallog << std::endl;
      }
}


template <int dim>
void
print_dof_points(const DoFHandler<dim> &dofh)
{
  hp::MappingCollection<dim> mapping;
  for (unsigned int i = 0; i < dofh.get_fe_collection().size(); ++i)
    mapping.push_back(MappingFE<dim>(dofh.get_fe(i)));

  std::vector<Point<dim>> points(dofh.n_dofs());
  DoFTools::map_dofs_to_support_points(mapping, dofh, points);

  for (unsigned int i = 0; i < dofh.n_dofs(); ++i)
    deallog << "dof:" << i << " point:" << points[i] << std::endl;
}


// ----- test -----

template <int dim>
void
test(const hp::FECollection<dim> &fes)
{
  // setup grid
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

#if false
  GridOut grid_out;
  grid_out.write_vtk(tria, deallog.get_file_stream());
#endif

  DoFHandler<dim> dofh(tria);
  dofh.begin_active()->set_active_fe_index(1);

  dofh.distribute_dofs(fes);
  deallog << "ndofs: " << dofh.n_dofs() << std::endl;

#if false
  print_dof_points(dofh);
  print_dof_indices_on_faces(dofh);
#endif

  // hanging node constraints
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>(hp::FECollection<2>(Simplex::FE_P<2>(1), Simplex::FE_P<2>(2)));
  test<2>(hp::FECollection<2>(Simplex::FE_P<2>(2), Simplex::FE_P<2>(1)));
  deallog.pop();
}

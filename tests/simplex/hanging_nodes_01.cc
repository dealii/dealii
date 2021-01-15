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



// verify hanging node constraints on locally h-refined simplex mesh
//
// dofs will be enumerated as follows
//  scenario 1:
//   1-------0
//   |\      |
//   |  \    |
//   5---6   |
//   |\  |\  |
//   |  \|  \|
//   3---4---2


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

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

        if (face->has_children())
          {
            for (unsigned int sf = 0; sf < face->n_children(); ++sf)
              {
                const auto &subface = face->child(sf);
                Assert(subface->n_active_fe_indices() == 1, ExcInternalError());
                const unsigned int subface_fe_index =
                  subface->nth_active_fe_index(0);
                const auto &subface_fe = dofh.get_fe(subface_fe_index);

                dof_indices.resize(subface_fe.n_dofs_per_face(f));
                subface->get_dof_indices(dof_indices, subface_fe_index);

                deallog << "cell:" << cell->active_cell_index() << " face:" << f
                        << " subface:" << sf << " dofs:";
                for (const auto &i : dof_indices)
                  deallog << i << " ";
                deallog << std::endl;
              }
          }
        else
          {
            Assert(face->n_active_fe_indices() == 1, ExcInternalError());
            const unsigned int face_fe_index = face->nth_active_fe_index(0);
            const auto &       face_fe       = dofh.get_fe(face_fe_index);

            dof_indices.resize(face_fe.n_dofs_per_face(f));
            face->get_dof_indices(dof_indices, face_fe_index);

            deallog << "cell:" << cell->active_cell_index() << " face:" << f
                    << " dofs:";
            for (const auto &i : dof_indices)
              deallog << i << " ";
            deallog << std::endl;
          }
      }
}


template <int dim>
void
print_dof_points(const DoFHandler<dim> &dofh)
{
  std::vector<Point<dim>> points(dofh.n_dofs());
  DoFTools::map_dofs_to_support_points(MappingFE<dim>(dofh.get_fe()),
                                       dofh,
                                       points);

  for (unsigned int i = 0; i < dofh.n_dofs(); ++i)
    deallog << "dof:" << i << " point:" << points[i] << std::endl;
}



// ----- test -----

template <int dim>
void
test()
{
  // setup grid
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

#if false
  GridOut grid_out;
  grid_out.write_vtk(tria, deallog.get_file_stream());
#endif

  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(Simplex::FE_P<dim>(1));
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
  test<2>();
  deallog.pop();
}

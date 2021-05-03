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



// Check MatrixFree::boundary_id() for some large boundary id numbers
// (exceeding a very old case of 8 bit integers)


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  unsigned int count = 0;
  for (unsigned int face = 0; face < 2 * dim; ++face)
    {
      tria.begin()->face(face)->set_boundary_id(count);
      count += 100000;
    }

  MappingQGeneric<dim> mapping(1);
  FE_Q<dim>            fe(1);
  DoFHandler<dim>      dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  MatrixFree<dim>                          mf_data;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(mapping, dof, constraints, QGauss<1>(2), data);
  for (unsigned int i = 0; i < mf_data.n_boundary_face_batches(); ++i)
    deallog << "Face "
            << static_cast<unsigned int>(
                 mf_data.get_face_info(i).interior_face_no)
            << " boundary id " << mf_data.get_boundary_id(i) << std::endl;
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}

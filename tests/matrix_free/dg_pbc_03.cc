// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test matrix-free implementation of face identification in 1d, 2d, 3d with
// adaptivity and periodic boundary conditions, using a similar test tas he
// dg_pbc_02 test, but for the active cells

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"



template <int dim>
void
test(const bool periodicity, const bool adaptive)
{
  deallog.push(std::to_string(dim) + "d");
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  if (periodicity)
    {
      deallog.push("periodic");
      tria.begin()->face(0)->set_all_boundary_ids(10);
      tria.begin()->face(1)->set_all_boundary_ids(11);
      if (dim == 3)
        {
          tria.begin()->face(4)->set_all_boundary_ids(12);
          tria.begin()->face(5)->set_all_boundary_ids(13);
        }

      std::vector<
        GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
        periodic_faces;
      GridTools::collect_periodic_faces(tria, 10, 11, 0, periodic_faces);
      if (dim == 3)
        GridTools::collect_periodic_faces(tria, 12, 13, 2, periodic_faces);

      tria.add_periodicity(periodic_faces);
    }
  else
    deallog.push("boundary");

  tria.refine_global(8 - 2 * dim);

  if (adaptive)
    {
      tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      if (dim == 1 && !periodicity)
        {
          // create a cell that differs by two levels
          tria.last()->set_refine_flag();
          tria.execute_coarsening_and_refinement();
        }
      deallog.push("adaptive");
    }
  else
    deallog.push("uniform ");

  FE_DGQ<dim>     fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  using MF = MatrixFree<dim, double, VectorizedArray<double, 1>>;
  const QGauss<1>             quad(1);
  typename MF::AdditionalData data;
  data.tasks_parallel_scheme = MF::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  deallog << "Active cells: " << tria.n_active_cells() << std::endl;
  data.mg_level = numbers::invalid_unsigned_int;

  MF mf_data;
  mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  std::vector<unsigned int> n_inner_faces(2 * dim),
    n_inner_other_faces(2 * dim), n_boundary_faces(2 * dim);
  for (unsigned int f = 0; f < mf_data.n_inner_face_batches(); ++f)
    {
      if (mf_data.get_face_info(f).cells_interior[0] !=
          numbers::invalid_unsigned_int)
        {
          n_inner_faces[mf_data.get_face_info(f).interior_face_no]++;
          n_inner_other_faces[mf_data.get_face_info(f).exterior_face_no]++;
        }
    }
  for (unsigned int f = mf_data.n_inner_face_batches();
       f < mf_data.n_inner_face_batches() + mf_data.n_boundary_face_batches();
       ++f)
    {
      if (mf_data.get_face_info(f).cells_interior[0] !=
          numbers::invalid_unsigned_int)
        {
          n_boundary_faces[mf_data.get_face_info(f).interior_face_no]++;
        }
    }

  deallog << "Cell batches: " << mf_data.n_cell_batches() << std::endl;
  deallog << "Interior faces: ";
  for (unsigned int f = 0; f < n_inner_faces.size(); ++f)
    deallog << n_inner_faces[f] << ' ';
  deallog << std::endl << "Exterior faces: ";
  for (unsigned int f = 0; f < n_inner_other_faces.size(); ++f)
    deallog << n_inner_other_faces[f] << ' ';
  deallog << std::endl << "Boundary faces: ";
  for (unsigned int f = 0; f < n_boundary_faces.size(); ++f)
    deallog << n_boundary_faces[f] << ' ';
  deallog << std::endl;
  deallog.pop();
  deallog.pop();
  deallog.pop();
  deallog << std::endl;
}



int
main()
{
  initlog();
  test<1>(false, false);
  test<1>(true, false);
  test<1>(false, true);
  test<1>(true, true);
  test<2>(false, false);
  test<2>(true, false);
  test<2>(false, true);
  test<2>(true, true);
  test<3>(false, false);
  test<3>(true, false);
  test<3>(false, true);
  test<3>(true, true);
  return 0;
}

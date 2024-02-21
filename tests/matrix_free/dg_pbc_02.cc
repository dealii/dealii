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



// tests that the matrix-free implementation works correctly with periodic
// boundary conditions by counting the number of faces in the different
// categories

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  std::vector<unsigned int> refinements(dim, 1);
  refinements[0] = 2;
  Point<dim> p2;
  p2[0] = 2;
  for (unsigned int d = 1; d < dim; ++d)
    p2[d] = 1;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            refinements,
                                            Point<dim>(),
                                            p2);

  tria.begin()->face(0)->set_all_boundary_ids(10);
  tria.last()->face(1)->set_all_boundary_ids(11);
  if (dim == 3)
    {
      tria.begin()->face(4)->set_all_boundary_ids(12);
      tria.begin()->face(5)->set_all_boundary_ids(13);
      tria.last()->face(4)->set_all_boundary_ids(12);
      tria.last()->face(5)->set_all_boundary_ids(13);
    }

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(tria, 10, 11, 0, periodic_faces);
  if (dim == 3)
    GridTools::collect_periodic_faces(tria, 12, 13, 2, periodic_faces);

  tria.add_periodicity(periodic_faces);

  tria.refine_global(8 - 2 * dim);

  FE_DGQ<dim>     fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();
  AffineConstraints<double> constraints;
  constraints.close();

  const QGauss<1>                          quad(1);
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    {
      MatrixFree<dim> mf_data;
      data.mg_level = level;
      mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
      std::vector<unsigned int> n_inner_faces(2 * dim),
        n_inner_other_faces(2 * dim), n_boundary_faces(2 * dim);
      for (unsigned int f = 0; f < mf_data.n_inner_face_batches(); ++f)
        {
          for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
            if (mf_data.get_face_info(f).cells_interior[v] !=
                numbers::invalid_unsigned_int)
              {
                n_inner_faces[mf_data.get_face_info(f).interior_face_no]++;
                n_inner_other_faces[mf_data.get_face_info(f)
                                      .exterior_face_no]++;
              }
        }
      for (unsigned int f = mf_data.n_inner_face_batches();
           f <
           mf_data.n_inner_face_batches() + mf_data.n_boundary_face_batches();
           ++f)
        {
          for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
            if (mf_data.get_face_info(f).cells_interior[v] !=
                numbers::invalid_unsigned_int)
              {
                n_boundary_faces[mf_data.get_face_info(f).interior_face_no]++;
              }
        }

      deallog << "Level: " << level << std::endl;
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
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                            argv,
                                            testing_max_num_threads());
  MPILogInitAll                    log;
  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2>();
    deallog.pop();
    deallog.push("3d");
    test<3>();
    deallog.pop();
  }
}

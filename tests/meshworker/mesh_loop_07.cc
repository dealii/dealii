// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Check that mesh loop returns the correct values for integral over
// volumes and boundaries

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include "../tests.h"


template <int dim, int spacedim = dim>
void
run()
{
  deallog << "Dim: " << dim << std::endl;

  const FE_Q<dim, spacedim> fe(1);
  const QGauss<dim>         cell_quadrature(fe.degree + 1);
  const QGauss<dim - 1>     face_quadrature(fe.degree + 1);

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::subdivided_hyper_cube(triangulation, 4, 0.0, 1.0);

  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
  using CopyData         = MeshWorker::CopyData<1, 1, 1>;
  using CellIteratorType = decltype(dof_handler.begin_active());

  // Volume integral
  {
    ScratchData scratch(dof_handler.get_fe(),
                        cell_quadrature,
                        update_JxW_values);
    CopyData    copy(1);

    auto cell_worker = [](const CellIteratorType &cell,
                          ScratchData &           scratch_data,
                          CopyData &              copy_data) {
      const auto &fe_values = scratch_data.reinit(cell);

      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
           ++q_point)
        {
          copy_data.vectors[0][0] += 1.0 * fe_values.JxW(q_point);
        }
    };

    double volume = 0.0;
    auto   copier = [&volume](const CopyData &copy_data) {
      volume += copy_data.vectors[0][0];
    };

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          scratch,
                          copy,
                          MeshWorker::assemble_own_cells);
    deallog << "Volume: " << volume << std::endl;

    double reference_volume = 0.0;
    for (auto &cell : dof_handler.active_cell_iterators())
      reference_volume += cell->measure();

    Assert(std::abs(volume - reference_volume) < 1e-6,
           ExcMessage("Volumes do not match. Reference value: " +
                      Utilities::to_string(reference_volume) +
                      "; Calculated value: " + Utilities::to_string(volume)));
  }

  // Boundary integral
  {
    ScratchData scratch(dof_handler.get_fe(),
                        cell_quadrature,
                        update_default,
                        face_quadrature,
                        update_JxW_values);
    CopyData    copy(1);

    std::function<void(const CellIteratorType &, ScratchData &, CopyData &)>
      empty_cell_worker;

    auto boundary_worker = [](const CellIteratorType &cell,
                              const unsigned int      face,
                              ScratchData &           scratch_data,
                              CopyData &              copy_data) {
      const auto &fe_face_values = scratch_data.reinit(cell, face);

      for (unsigned int q_point = 0;
           q_point < fe_face_values.n_quadrature_points;
           ++q_point)
        {
          copy_data.vectors[0][0] += 1.0 * fe_face_values.JxW(q_point);
        }
    };

    double area   = 0.0;
    auto   copier = [&area](const CopyData &copy_data) {
      // ***
      // The problem is observed here. For some faces, the copy_data that
      // is passed in is reset, even though there is always work being
      // done in the boundary_worker.
      // ***
      area += copy_data.vectors[0][0];
    };

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          empty_cell_worker,
                          copier,
                          scratch,
                          copy,
                          MeshWorker::assemble_boundary_faces,
                          boundary_worker);
    deallog << "Area: " << area << std::endl;

    double reference_area = 0.0;
    for (auto &cell : dof_handler.active_cell_iterators())
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        if (cell->face(face)->at_boundary())
          {
            reference_area += cell->face(face)->measure();
          }

    Assert(std::abs(area - reference_area) < 1e-6,
           ExcMessage("Areas do not match. Reference value: " +
                      Utilities::to_string(reference_area) +
                      "; Calculated value: " + Utilities::to_string(area)));
  }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  run<2>();
  run<3>();

  deallog << "OK" << std::endl;
}

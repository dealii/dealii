// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GridTools::MarchingCubeAlgorithm.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



using VectorType = LinearAlgebra::distributed::Vector<double>;


namespace dealii
{
  namespace GridGenerator
  {
    template <int dim, typename VectorType>
    void
    create_triangulation_with_marching_cube_algorithm(
      Triangulation<dim - 1, dim> &tria,
      const Mapping<dim>          &mapping,
      const DoFHandler<dim>       &background_dof_handler,
      const VectorType            &ls_vector,
      const double                 iso_level,
      const unsigned int           n_subdivisions,
      const double                 tolerance = 1e-10)
    {
      std::vector<Point<dim>>        vertices;
      std::vector<CellData<dim - 1>> cells;
      SubCellData                    subcelldata;

      const GridTools::MarchingCubeAlgorithm<dim, VectorType> mc(
        mapping, background_dof_handler.get_fe(), n_subdivisions, tolerance);

      mc.process(background_dof_handler, ls_vector, iso_level, vertices, cells);

      std::vector<unsigned int> considered_vertices;

      // note: the following operation does not work for simplex meshes yet
      // GridTools::delete_duplicated_vertices (vertices, cells, subcelldata,
      // considered_vertices);

      tria.create_triangulation(vertices, cells, subcelldata);
    }
  } // namespace GridGenerator
} // namespace dealii

template <int dim>
class LSFunction : public Function<dim>
{
public:
  LSFunction()
    : Function<dim>(1)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    const double v = p.distance(Point<dim>()) - 0.5;

    if (v > +0.1)
      return +0.1;

    if (v < -0.1)
      return -0.1;

    return v;
  }
};



template <typename MeshType>
MPI_Comm
get_mpi_comm(const MeshType &mesh)
{
  const auto *tria_parallel = dynamic_cast<
    const parallel::TriangulationBase<MeshType::dimension,
                                      MeshType::space_dimension> *>(
    &(mesh.get_triangulation()));

  return tria_parallel != nullptr ? tria_parallel->get_mpi_communicator() :
                                    MPI_COMM_SELF;
}



template <int dim, int spacedim>
std::shared_ptr<const Utilities::MPI::Partitioner>
create_partitioner(const DoFHandler<dim, spacedim> &dof_handler)
{
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  return std::make_shared<const Utilities::MPI::Partitioner>(
    dof_handler.locally_owned_dofs(),
    locally_relevant_dofs,
    get_mpi_comm(dof_handler));
}



template <int dim>
void
test(const unsigned int n_subdivisions, const double iso_level)
{
  const unsigned int spacedim             = dim + 1;
  const unsigned int background_fe_degree = 3;

  parallel::shared::Triangulation<spacedim> background_tria(
    MPI_COMM_WORLD, Triangulation<spacedim>::none, true);

  if (spacedim == 2)
    GridGenerator::subdivided_hyper_rectangle(background_tria,
                                              {20, 20},
                                              {-1.0, -1.0},
                                              {+1.0, +1.0});
  else
    GridGenerator::subdivided_hyper_rectangle(background_tria,
                                              {20, 20, 20},
                                              {-1.0, -1.0, -1.0},
                                              {+1.0, +1.0, +1.0});

  FE_Q<spacedim>       background_fe(background_fe_degree);
  DoFHandler<spacedim> background_dof_handler(background_tria);
  background_dof_handler.distribute_dofs(background_fe);

  MappingQ1<spacedim> background_mapping;

  VectorType ls_vector(create_partitioner(background_dof_handler));

  VectorTools::interpolate(background_mapping,
                           background_dof_handler,
                           LSFunction<spacedim>(),
                           ls_vector);

  ls_vector.update_ghost_values();

  parallel::shared::Triangulation<dim, spacedim> tria(
    MPI_COMM_WORLD, Triangulation<dim, spacedim>::none, true);

  GridGenerator::create_triangulation_with_marching_cube_algorithm(
    tria,
    background_mapping,
    background_dof_handler,
    ls_vector,
    iso_level,
    n_subdivisions);

  // write computed vectors to Paraview
  if (true /*write surface mesh*/)
    {
      GridOut grid_out;
      grid_out.write_vtk(tria, deallog.get_file_stream());
    }

  if (false /*write background mesh*/)
    {
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;

      DataOut<spacedim> data_out;
      data_out.set_flags(flags);
      data_out.add_data_vector(background_dof_handler, ls_vector, "curvature");

      data_out.build_patches(background_mapping, background_fe_degree + 1);
      data_out.write_vtu_with_pvtu_record(
        "./", "data_background." + std::to_string(spacedim), 0, MPI_COMM_WORLD);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  for (unsigned int i = 1; i <= 3; ++i)
    test<1>(i, -0.1 + i * 0.05);
  for (unsigned int i = 1; i <= 3; ++i)
    test<2>(i, -0.1 + i * 0.05);
}

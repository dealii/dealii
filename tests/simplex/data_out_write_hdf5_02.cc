// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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



// Test DataOut with HDF5 for simplex meshes in parallel.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p[component % dim] * p[component % dim];
  }
};

template <int dim, int spacedim = dim>
void
test(const FiniteElement<dim, spacedim> &fe, const unsigned int n_components)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  parallel::fullydistributed::Triangulation<dim> tria(comm);

  {
    Triangulation<dim, spacedim> tria_serial;
    GridGenerator::subdivided_hyper_cube_with_simplices(tria_serial,
                                                        dim == 2 ? 4 : 2);

    GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                       tria_serial);

    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(tria_serial, comm);

    tria.create_triangulation(construction_data);
  }

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  IndexSet owned_dofs = dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  LinearAlgebra::distributed::Vector<double> solution;

  solution.reinit(owned_dofs, locally_relevant_dofs, comm);

  MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  VectorTools::interpolate(mapping,
                           dof_handler,
                           RightHandSideFunction<dim>(n_components),
                           solution);
  solution.update_ghost_values();

  static unsigned int counter = 0;

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches(mapping);

  const std::string output_basename("test." + std::to_string(dim) + "." +
                                    std::to_string(counter++));

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(true, true));
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter, output_basename + ".h5", comm);

  std::vector<XDMFEntry> xdmf_entries({data_out.create_xdmf_entry(
    data_filter, output_basename + ".h5", 0, comm)});

  data_out.write_xdmf_file(xdmf_entries, output_basename + ".xdmf", comm);

  data_out.clear();

  // Sadly HDF5 is binary and we can not use hd5dump because it might not be
  // in the path. At least we can make sure that both the .xdmf and the .h5
  // files are created.
  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      std::ifstream h5((output_basename + ".h5").c_str());
      AssertThrow(h5.good(), ExcIO());

      std::ifstream xdmf((output_basename + ".xdmf").c_str());
      AssertThrow(h5.good(), ExcIO());

      deallog << "Files " << output_basename + ".h5"
              << " and " << output_basename + ".xdmf"
              << " created succesfully!" << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  {
    const unsigned int dim = 2;
    test<dim>(FE_SimplexP<dim>(2), 1);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2), dim), dim);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2), dim, FE_SimplexP<dim>(1), 1),
              dim + 1);
  }
  {
    const unsigned int dim = 3;
    test<dim>(FE_SimplexP<dim>(2), 1);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2), dim), dim);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2), dim, FE_SimplexP<dim>(1), 1),
              dim + 1);
  }
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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



// Test DataOut with HDF5 for simplex meshes.

#include <deal.II/base/mpi.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
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
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, dim == 2 ? 4 : 2);

  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());

  MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  VectorTools::interpolate(mapping,
                           dof_handler,
                           RightHandSideFunction<dim>(n_components),
                           solution);

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
  data_out.write_hdf5_parallel(data_filter,
                               output_basename + ".h5",
                               MPI_COMM_SELF);

  std::vector<XDMFEntry> xdmf_entries({data_out.create_xdmf_entry(
    data_filter, output_basename + ".h5", 0, MPI_COMM_SELF)});

  data_out.write_xdmf_file(xdmf_entries,
                           output_basename + ".xdmf",
                           MPI_COMM_SELF);

  data_out.clear();

  // Sadly hdf5 is binary and we can not use hd5dump because it might
  // not be in the path. At least we can make sure that both the xdmf and
  // the h5 file are created.
  std::ifstream h5((output_basename + ".h5").c_str());
  AssertThrow(h5.good(), ExcIO());

  std::ifstream xdmf((output_basename + ".xdmf").c_str());
  AssertThrow(h5.good(), ExcIO());

  deallog << "Files " << output_basename + ".h5"
          << " and " << output_basename + ".xdmf"
          << " created succesfully!" << std::endl;
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

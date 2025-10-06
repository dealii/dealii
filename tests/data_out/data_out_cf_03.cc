// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests DataOut with HDF5

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

#include "./data_out_cf_common.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  if constexpr (dim == 2)
    {
      std::vector points{
        Point<dim>{0, 0},
        Point<dim>{0, 1},
        Point<dim>{1, 0},
        Point<dim>{1, 1},
      };
      std::vector cells{2, CellData<dim>(3)};
      cells[0].vertices.assign({0, 1, 2});
      cells[1].vertices.assign({1, 3, 2});
      tria.create_triangulation(points, cells, {});
    }
  else
    {
      AssertThrow(false, StandardExceptions::ExcNotImplemented());
    }
  tria.refine_global(1);

  FE_SimplexP<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  Vector<double> v1(dof1.n_dofs());
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = i;

  DataOut<dim> data_out;
  data_out.add_data_vector(dof1, v1, "linear");
  data_out.build_patches();

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(/*filter_vertices*/ true));
  data_out.write_filtered_data(data_filter);
  data_out.write_cf_parallel(data_filter, "out.nc", MPI_COMM_SELF);

  deallog << "ok" << std::endl;
  dump_nc_file(deallog, "out.nc");
}

int
main(int argc, char *argv[])
{
  initlog();

  try
    {
      test<2>();

      return 0;
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

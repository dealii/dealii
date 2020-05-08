// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// test VectorTools::interpolate_to_different_mesh in parallel


namespace LA
{
  using namespace dealii::LinearAlgebraTrilinos;
}

template <int dim>
class SomeFunction : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int) const
  {
    return 1 + p(0) * p(0);
  }
};

template <int dim>
void
setup(DoFHandler<dim> &dh,
      FE_Q<dim> &      fe,
      LA::MPI::Vector &vec,
      LA::MPI::Vector &lr_vec)
{
  dh.distribute_dofs(fe);
  vec.reinit(dh.locally_owned_dofs(), MPI_COMM_WORLD);
  IndexSet locally_relevant;
  DoFTools::extract_locally_relevant_dofs(dh, locally_relevant);
  lr_vec.reinit(locally_relevant, MPI_COMM_WORLD);
}

template <int dim>
void
output(DoFHandler<dim> & dh,
       LA::MPI::Vector & v,
       unsigned int      loop,
       const std::string filename_)
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  DataOut<dim> data_out;
  data_out.add_data_vector(dh, v, "1");
  data_out.build_patches(1);
  std::ostringstream filename;
  filename << filename_ << Utilities::int_to_string(loop, 2) << "."
           << Utilities::int_to_string(myid, 2) << ".vtu";

  std::ofstream output(filename.str().c_str());
  data_out.write_vtu(output);
  if (myid)
    {
      std::vector<std::string> filenames;
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
           ++i)
        filenames.push_back(filename_ + Utilities::int_to_string(loop, 2) +
                            "." + Utilities::int_to_string(i, 2) + ".vtu");
      const std::string pvtu_master_filename =
        (filename_ + Utilities::int_to_string(loop, 2) + ".pvtu");
      std::ofstream pvtu_master(pvtu_master_filename.c_str());
      data_out.write_pvtu_record(pvtu_master, filenames);
    }
}

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim, dim>::none,
    parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  parallel::distributed::Triangulation<dim> tr2(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim, dim>::none,
    parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

  GridGenerator::hyper_cube(tr2);
  tr2.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dh(tr);
  DoFHandler<dim> dh2(tr2);

  SomeFunction<dim> func;


  for (unsigned int loop = 0; loop < 5; ++loop)
    {
      // randomly refine:
      std::vector<bool> r_flags(tr.n_active_cells() * dim, false);
      std::vector<bool> c_flags(tr.n_active_cells(), false);

      for (unsigned int i = 0; i < c_flags.size(); ++i)
        {
          int roll = Testing::rand() % 4;
          if (roll >= 2)
            {
              for (unsigned int j = 0; j < dim; ++j)
                r_flags[i * dim + j] = true;
            }
          else if (roll == 1)
            c_flags[i] = true;
        }

      tr.load_coarsen_flags(c_flags);
      tr.load_refine_flags(r_flags);

      tr.execute_coarsening_and_refinement();
      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / " << tr.n_global_active_cells() << std::endl;


      LA::MPI::Vector vec1;
      LA::MPI::Vector lr_vec1;
      setup(dh, fe, vec1, lr_vec1);

      LA::MPI::Vector vec2;
      LA::MPI::Vector lr_vec2;
      setup(dh2, fe, vec2, lr_vec2);

      // interpolate func on old mesh:
      VectorTools::interpolate(dh2, func, vec2);
      lr_vec2 = vec2;

      // interpolate from vec2 to vec1
      VectorTools::interpolate_to_different_mesh(dh2, lr_vec2, dh, vec1);
      lr_vec1 = vec1;

      {
        Vector<double> local_errors(tr.n_active_cells());
        VectorTools::integrate_difference(dh,
                                          lr_vec1,
                                          func,
                                          local_errors,
                                          QGauss<dim>(3),
                                          VectorTools::L2_norm);
        double       total_local_error = local_errors.l2_norm();
        const double total_global_error =
          std::sqrt(Utilities::MPI::sum(total_local_error * total_local_error,
                                        MPI_COMM_WORLD));
        if (myid == 0)
          deallog << "err: " << total_global_error << std::endl;
      }

      // output(dh, lr_vec1, loop, "solutionA-");
      // output(dh2, lr_vec2, loop, "solutionB-");

      // also update tr2 to be the same as tr1
      tr2.load_coarsen_flags(c_flags);
      tr2.load_refine_flags(r_flags);
      tr2.execute_coarsening_and_refinement();

      // repartition both in the same way
      tr.repartition();
      tr2.repartition();

      // checks:
      const unsigned int checksum  = tr.get_checksum();
      const unsigned int checksum2 = tr2.get_checksum();
      if (myid == 0)
        deallog << "Checksum: " << checksum << " " << checksum2 << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  test<2>();
}

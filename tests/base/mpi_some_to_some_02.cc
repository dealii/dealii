// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for the internal::send_and_receive function in the case of
// collective communication

#include <deal.II/base/mpi.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <tuple>
#include <vector>

#include "../tests.h"

using namespace Patterns::Tools;

inline unsigned int
random_index(const unsigned int &max_index)
{
  return Testing::rand() * (max_index - 1) / RAND_MAX;
}

template <typename T>
std::string
to_string(const T &object)
{
  return Convert<T>::to_string(object);
}

template <int dim>
void
test(const unsigned int max_particles,
     const double       shared_fraction,
     const unsigned int max_cell_levels,
     const unsigned int max_cell_index)
{
  auto n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // Test some fake structure, similar to particles
  using particle = typename std::pair<unsigned int, Point<dim>>;


  unsigned int          n_local_particles = random_index(max_particles);
  std::vector<particle> particles(n_local_particles);

  auto all_n_local_particles =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, n_local_particles);

  AssertDimension(all_n_local_particles.size(), n_procs);

  // Compute the start of my global indices
  unsigned int my_global_start = 0;
  for (unsigned int i = 0; i < my_proc; ++i)
    my_global_start += all_n_local_particles[i];

  for (unsigned int i = 0; i < n_local_particles; ++i)
    particles[i] = std::make_pair(my_global_start + i, random_point<dim>());

  std::set<unsigned int>                        shared_indices;
  std::set<unsigned int>                        shared_procs_set;
  std::map<unsigned int, std::vector<particle>> shared_particles;

  for (unsigned int i = 0; i < n_procs; ++i)
    {
      auto rank = random_index(n_procs);
      if (rank != my_proc)
        shared_procs_set.insert(rank);
    }

  if (shared_procs_set.size())
    {
      std::vector<unsigned int> shared_procs(shared_procs_set.begin(),
                                             shared_procs_set.end());

      deallog << "Proc " << my_proc << '/' << n_procs << ": sharing with  "
              << Patterns::Tools::to_string(shared_procs_set) << std::endl;
      for (unsigned int i = 0; i < n_local_particles * shared_fraction; ++i)
        shared_particles[shared_procs[random_index(shared_procs.size())]]
          .push_back(*(particles.begin() + i));
    }

  auto received_particles =
    Utilities::MPI::some_to_some(MPI_COMM_WORLD, shared_particles);

  auto original_particles =
    Utilities::MPI::some_to_some(MPI_COMM_WORLD, received_particles);

  // now check that shared_particles and original_particles are the same

  if (original_particles == shared_particles)
    deallog << "OK" << std::endl;
  else
    deallog << "Not OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;


  test<2>(40, .3, 5, 10);
}

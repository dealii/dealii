// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


// Purpose of this test:
// Demonstrate that mesh smoothing option limit_level_difference_at_vertices
// is always effective across different initial cells (different trees in p4est
// context) not matter it is set or not.
//
// Design of this test:
// The initial mesh is set to 4 cells in unite square. Then the center most cell
// in each level is refined consecutively for 4 times. Finally, total active
// cell number in each refine level and cell center coordinates of the final
// active cells are compared for cases with and without flag mesh smoothing
// limit_level_difference_at_vertices.
//
// Expected result:
// In principle, the result of the two refinement runs should be different.
// However, under the current implementation, their results are the same.
// So, if some day this test failed, it might be a good news.

// Now the output of mesh and refine indicators in vtu format is disabled to
// prevent unnecessary outputs during regression test. It you want to see the
// mesh, uncomment the following macro definition.
// #define WRITE_VTU

template <int dim>
class Location;

template <int dim>
class TriaTest
{
private:
  using TypeTria = typename parallel::distributed::Triangulation<dim>;

public:
  TriaTest(const typename dealii::Triangulation<dim>::MeshSmoothing
             smoothing_option = dealii::Triangulation<dim>::none);
  void
  run(std::vector<unsigned int> &n_cell,
      std::set<Location<dim>>   &position_list);

private:
  void
  write_vtu(const unsigned int counter) const;

  const MPI_Comm     mpi_communicator;
  TypeTria           triangulation;
  const unsigned int myid;
  const bool         I_am_host;
  const std::string  case_name;
};

// Equip class Point<dim> with operator < in order to sort the cell center
// coordinates.
template <int dim>
class Location : public Point<dim>
{
public:
  Location(const Point<dim> &in)
    : Point<dim>(in)
  {}

  inline bool
  operator<(const Location<dim> &op) const
  {
    for (unsigned int d = 0; d < dim; ++d)
      {
        if ((*this)[d] == op[d])
          {
            continue;
          }
        return ((*this)[d] < op[d]);
      }
    return (false);
  }
};

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, /* int max_num_threads */ 1);
  const bool I_am_host =
    (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
  // Although most part of this test is designed to run in parallel, there is
  // still one place that doesn't work perfectly in parallel. Now the sorting of
  // final cell centers is just a local operation and won't produce any
  // reasonable result in parallel.
  AssertThrow(I_am_host,
              ExcMessage("Current code works properly only with one process."));

  const unsigned int dim = 2;
  if (I_am_host)
    {
      initlog();
      deallog.push("2d");
    }

  std::vector<unsigned int> n_cell_smooth;
  std::vector<unsigned int> n_cell_no_smooth;

  std::set<Location<dim>> final_cell_center_locations_smooth;
  std::set<Location<dim>> final_cell_center_locations_no_smooth;
  if (I_am_host)
    {
      deallog << "Flag limit_level_difference_at_vertices set:" << std::endl;
    }
  {
    TriaTest<dim> tria_test(
      dealii::Triangulation<dim>::limit_level_difference_at_vertices);
    tria_test.run(n_cell_smooth, final_cell_center_locations_smooth);
  }
  if (I_am_host)
    {
      deallog << "Flag limit_level_difference_at_vertices unset:" << std::endl;
    }
  {
    TriaTest<dim> tria_test(dealii::Triangulation<dim>::none);
    tria_test.run(n_cell_no_smooth, final_cell_center_locations_no_smooth);
  }
  if (I_am_host)
    {
      deallog << "Compare result:" << std::endl;
    }
  {
    bool n_cells_are_same = (n_cell_smooth.size() == n_cell_no_smooth.size());

    for (unsigned int i = 0; n_cells_are_same && (i < n_cell_smooth.size());
         ++i)
      {
        n_cells_are_same =
          n_cells_are_same && (n_cell_smooth[i] == n_cell_no_smooth[i]);
      }
    if (I_am_host)
      {
        deallog << "n_cells_are_same = " << n_cells_are_same << std::endl;
      }
  }
  {
    bool cell_center_locations_are_same =
      (final_cell_center_locations_smooth.size() ==
       final_cell_center_locations_no_smooth.size());

    std::set<Location<dim>>::const_iterator it1 =
      final_cell_center_locations_smooth.begin();
    std::set<Location<dim>>::const_iterator it2 =
      final_cell_center_locations_no_smooth.begin();

    for (; cell_center_locations_are_same &&
           (it1 != final_cell_center_locations_smooth.end());
         ++it1, ++it2)
      {
        cell_center_locations_are_same =
          cell_center_locations_are_same && (*it1 == *it2);
      }
    if (I_am_host)
      {
        deallog << "cell_center_locations_are_same = "
                << cell_center_locations_are_same << std::endl;
      }
  }

  deallog.pop();

  return (0);
}

template <int dim>
TriaTest<dim>::TriaTest(
  const typename dealii::Triangulation<dim>::MeshSmoothing smoothing_option)
  : mpi_communicator(MPI_COMM_WORLD)
  , triangulation(mpi_communicator, smoothing_option)
  , myid(Utilities::MPI::this_mpi_process(mpi_communicator))
  , I_am_host(myid == 0)
  , case_name(smoothing_option ? "smooth-" : "no_smooth-")
{
  std::vector<unsigned int> repetitions;
  Point<dim>                p1;
  Point<dim>                p2;

  for (unsigned int d = 0; d < dim; ++d)
    {
      repetitions.push_back(2);
      p1[d] = 0.0;
      p2[d] = 1.0;
    }
  GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1, p2);

  // GridGenerator::hyper_cube(triangulation);
  // triangulation.refine_global(1);
}

template <int dim>
void
TriaTest<dim>::run(std::vector<unsigned int> &n_cell,
                   std::set<Location<dim>>   &position_list)
{
  n_cell.clear();
  position_list.clear();

  unsigned int counter = 0;
  if (I_am_host)
    {
      deallog << "n_loop  n_cell" << std::endl;
      deallog << counter << "       " << triangulation.n_global_active_cells()
              << std::endl;
    }
  ++counter;

  for (; counter < 5; ++counter)
    {
      {
        Point<dim> p;
        for (unsigned int d = 0; d < dim; ++d)
          {
            p[d] = 0.5 - std::pow(0.5, 1.0 + counter);
          }
        typename TypeTria::active_cell_iterator cell =
          triangulation.begin_active();
        for (; cell != triangulation.end(); ++cell)
          if (cell->is_locally_owned() && ((cell->center()).distance(p) < 1e-4))
            {
              cell->set_refine_flag();
            }
      }

      triangulation.prepare_coarsening_and_refinement();
      write_vtu(counter);
      triangulation.execute_coarsening_and_refinement();

      if (I_am_host)
        {
          deallog << counter << "       "
                  << triangulation.n_global_active_cells() << std::endl;
        }
      n_cell.push_back(triangulation.n_global_active_cells());
    }

  // Output cell center coordinates in a sorted manner to prevent change in
  // cell ordering causing failure of this test.
  // This is just a local operation, you have to re-implement this part if you
  // want to run this test in distributed parallel.
  {
    deallog << " position of cell centers:" << std::endl;

    for (typename TypeTria::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      if (cell->is_locally_owned())
        {
          Location<dim> loc(cell->center());
          position_list.insert(loc);
        }

    for (typename std::set<Location<dim>>::const_iterator it =
           position_list.begin();
         it != position_list.end();
         ++it)
      {
        deallog << *it << std::endl;
      }
  }

  return;
}

template <int dim>
void
TriaTest<dim>::write_vtu(const unsigned int counter) const
{
#ifdef WRITE_VTU
  // save refine flag
  Vector<float> refine_mark;
  {
    refine_mark.reinit(triangulation.n_active_cells());

    typename TypeTria::active_cell_iterator cell = triangulation.begin_active();
    const typename TypeTria::active_cell_iterator endc = triangulation.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          refine_mark[cell->active_cell_index()] = cell->refine_flag_set();
        }
  }

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);
  {
    const std::string data_name("refine_flag");
    data_out.add_data_vector(refine_mark,
                             data_name,
                             DataOut<dim>::type_cell_data);
  }

  data_out.build_patches();

  const std::string output_tag =
    case_name + Utilities::int_to_string(counter, 4);
  const std::string slot_itag = ".slot-" + Utilities::int_to_string(myid, 4);

  std::ofstream output(output_tag + slot_itag + ".vtu");
  data_out.write_vtu(output);

  if (I_am_host)
    {
      std::vector<std::string> filenames;
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
        {
          filenames.push_back(output_tag + ".slot-" +
                              Utilities::int_to_string(i, 4) + ".vtu");
        }
      std::ofstream pvtu_output(output_tag + ".pvtu");
      data_out.write_pvtu_record(pvtu_output, filenames);
    }
#else
  (void)counter;
#endif
  return;
}

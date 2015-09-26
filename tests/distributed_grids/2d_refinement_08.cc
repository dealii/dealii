// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/lac/vector.h>
#include <fstream>

using namespace dealii;

template<int dim>
class TriaTest
{
private:
  typedef typename parallel::distributed::Triangulation<dim> TypeTria;
public:
  TriaTest(const typename dealii::Triangulation<dim>::MeshSmoothing smoothing_option = dealii::Triangulation<dim>::none);
  void run();
private:
  void write_vtu(const unsigned int counter) const;

  const MPI_Comm mpi_communicator;
  TypeTria triangulation;
  const unsigned int myid;
  const bool I_am_host;
};

template<int dim>
class Location: public Point<dim>
{
public:
  Location(const Point<dim> &in)
    :
    Point<dim>(in)
  {}

  inline
  bool operator<(const Location<dim> &op) const
  {
    for (unsigned int d=0; d<dim; ++d)
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


int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, /* int max_num_threads */ 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  TriaTest<2> tria_test;
  tria_test.run();
  deallog.pop();

  return (0);
}

template<int dim>
TriaTest<dim>::TriaTest(const typename dealii::Triangulation<dim>::MeshSmoothing smoothing_option)
  :
  mpi_communicator(MPI_COMM_WORLD),
  triangulation(mpi_communicator, smoothing_option),
  myid (Utilities::MPI::this_mpi_process (mpi_communicator)),
  I_am_host(myid == 0)
{
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
}

template<int dim>
void TriaTest<dim>::run()
{

  unsigned int counter = 0;
  if (I_am_host)
    {
      deallog << "n_loop  n_cell" << std::endl;
      deallog << counter << "       "
              << triangulation.n_global_active_cells() << std::endl;
    }
  ++counter;

  for (; counter < 5; ++counter)
    {
      {
        Point<dim> p;
        for (unsigned int d=0; d<dim; ++d)
          {
            p[d] = 0.5 - std::pow(0.5, 1.0 + counter);
          }
        typename TypeTria::active_cell_iterator
        cell = triangulation.begin_active();
        for (; cell!= triangulation.end(); ++cell)
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
    }

  // Output cell center coordinates in a sorted manner to prevent change in
  // cell ordering causing failure of this test.
  {
    deallog << " position of cell centers:" << std::endl;
    std::set<Location<dim> > position_list;

    for (typename TypeTria::active_cell_iterator cell = triangulation.begin_active();
         cell!= triangulation.end();
         ++cell)
      if (cell->is_locally_owned())
        {
          Location<dim> loc(cell->center());
          position_list.insert(loc);
        }

    for (typename std::set<Location<dim> >::const_iterator it = position_list.begin();
         it != position_list.end();
         ++it)
      {
        deallog << *it << std::endl;
      }
  }

  return;
}

template<int dim>
void TriaTest<dim>::write_vtu(const unsigned int counter) const
{
  // save refine flag
  Vector<float> refine_mark;
  {
    refine_mark.reinit (triangulation.n_active_cells());

    typename TypeTria::active_cell_iterator
    cell = triangulation.begin_active();
    const typename TypeTria::active_cell_iterator
    endc = triangulation.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          refine_mark[cell->active_cell_index()] =
            cell->refine_flag_set();
        }
  }

  DataOut<dim> data_out;
  data_out.attach_triangulation (triangulation);
  {
    const std::string data_name ("refine_flag");
    data_out.add_data_vector (refine_mark,
                              data_name,
                              DataOut<dim>::type_cell_data);
  }

  data_out.build_patches();

  const std::string output_tag = "grid-" +
                                 Utilities::int_to_string (counter, 4);
  const std::string slot_itag = ".slot-" + Utilities::int_to_string (myid, 4);

  std::ofstream output ((output_tag + slot_itag + ".vtu").c_str());
  data_out.write_vtu (output);

  if (I_am_host)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0;
           i<Utilities::MPI::n_mpi_processes (mpi_communicator);
           ++i)
        {
          filenames.push_back (output_tag +
                               ".slot-" +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");
        }
      std::ofstream master_output ((output_tag + ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
  return;
}


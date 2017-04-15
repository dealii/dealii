// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
//

// test DoFTools::extract_dofs_with_support_contained_within() for calculation
// of the RHS function in parallel

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <set>
#include <cstdio>


std::string output_name(const unsigned int subdomain)
{
  return "output_" +
         Utilities::int_to_string(subdomain) +
         ".vtu";
}

template<int dim>
bool
pred_d(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return (cell->center()(0) < 0.49);
}

template<int dim>
bool
pred_r(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return (cell->center()(0) < 0.49 &&
          cell->center()(1) < 0.49) ||
         (cell->center()(0) > 0.49 &&
          cell->center()(1) > 0.49);
}



template <int dim>
void test ()
{

  // Setup system
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  GridGenerator::hyper_rectangle (triangulation,
                                  Point<dim>(0,0),
                                  Point<dim>(1,1));

  triangulation.refine_global(1);

  // Extra refinement to generate hanging nodes
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned() && pred_r<dim>(cell))
      cell->set_refine_flag ();

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement ();

  DoFHandler<dim> dh (triangulation);

  FE_Q<dim> fe(2);
  dh.distribute_dofs (fe);

  const IndexSet &locally_owned_set = dh.locally_owned_dofs();
  IndexSet locally_relevant_set;
  DoFTools::extract_locally_relevant_dofs (dh,
                                           locally_relevant_set);

  ConstraintMatrix cm;
  cm.reinit (locally_relevant_set);
  DoFTools::make_hanging_node_constraints(dh,cm);
  cm.close ();

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  // get support on the predicate
  IndexSet support = DoFTools::extract_dofs_with_support_contained_within(dh, pred_d<dim>, cm);
  IndexSet local_support = support & locally_owned_set;

  // rhs vectors:
  LinearAlgebra::distributed::Vector<double> sparse_rhs;
  LinearAlgebra::distributed::Vector<double> rhs;
  sparse_rhs.reinit(locally_owned_set, locally_relevant_set, MPI_COMM_WORLD);
  rhs.reinit(       locally_owned_set, locally_relevant_set, MPI_COMM_WORLD);

  sparse_rhs.zero_out_ghosts();
  rhs.zero_out_ghosts();

  // assemble RHS which has a local support:
  const std::function< double(const Point<dim> &)> rhs_func
  = [=] (const Point<dim> &p) -> double
  {
    return p[0] > 0.5 ? 0. : 0.5 - p[0];
  };

  Vector<double> local_rhs(fe.dofs_per_cell);
  QGauss<dim> quadrature(3);
  FEValues<dim> fe_values(fe, quadrature,
                          update_values | update_JxW_values | update_quadrature_points);
  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dh.begin_active();
       cell != dh.end(); ++cell)
    if (cell->is_locally_owned() && pred_d<dim>(cell))
      {
        fe_values.reinit(cell);
        cell->get_dof_indices (local_dof_indices);

        const std::vector<Point<dim>> &q_points = fe_values.get_quadrature_points();

        local_rhs = 0.;
        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            local_rhs[i] += fe_values.shape_value (i, q) *
                            rhs_func(q_points[q]) *
                            fe_values.JxW (q);

        cm.distribute_local_to_global (local_rhs,
                                       local_dof_indices,
                                       rhs);

        // copy-paste of CM distribute_local_to_global and
        // add is_element() checks:
        auto local_vector_begin   = local_rhs.begin();
        const auto local_vector_end = local_rhs.end();
        auto local_indices_begin = local_dof_indices.begin();
        const std::vector<std::pair<types::global_dof_index,double> > *line_ptr;
        for ( ; local_vector_begin != local_vector_end;
              ++local_vector_begin, ++local_indices_begin)
          {
            line_ptr = cm.get_constraint_entries(*local_indices_begin);
            if (line_ptr==NULL) // unconstrained
              {
                if (support.is_element(*local_indices_begin))
                  sparse_rhs(*local_indices_begin) += *local_vector_begin;
              }
            else
              {
                const unsigned int line_size = line_ptr->size();
                for (unsigned int j=0; j<line_size; ++j)
                  if (support.is_element((*line_ptr)[j].first))
                    sparse_rhs((*line_ptr)[j].first)
                    += *local_vector_begin * (*line_ptr)[j].second;
              }
          }
        /*
        cm.distribute_local_to_global (local_rhs,
                                       local_dof_indices,
                                       sparse_rhs);
        */
      }

  rhs.compress(VectorOperation::add);
  sparse_rhs.compress(VectorOperation::add);

  rhs.update_ghost_values();
  sparse_rhs.update_ghost_values();

  rhs.add(-1., sparse_rhs);

  // print grid and DoFs for visual inspection
  if (false)
    {
      const unsigned int this_mpi_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      const unsigned int n_mpi_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

      for (unsigned int i = 0; i < n_mpi_processes; ++i)
        {
          MPI_Barrier(MPI_COMM_WORLD);
          if (i == this_mpi_process)
            {
              std::cout << "-------------------- " << this_mpi_process << std::endl;
              cm.print(std::cout);

              std::cout << "local support:" << std::endl;
              local_support.print(std::cout);
              std::cout << "support:" << std::endl;
              support.print(std::cout);
              std::cout << "locally owned:" << std::endl;
              locally_owned_set.print(std::cout);
              std::cout << "locally relevant set:" << std::endl;
              locally_relevant_set.print(std::cout);
            }
        }

      {
        std::map<types::global_dof_index, Point<dim> > support_points;
        MappingQ1<dim> mapping;
        DoFTools::map_dofs_to_support_points(mapping, dh, support_points);

        const std::string filename =
          "grid" + dealii::Utilities::int_to_string(n_mpi_processes) + dealii::Utilities::int_to_string(this_mpi_process);
        std::ofstream f((filename+".gp").c_str());

        f << "set terminal png size 400,410 enhanced font \"Helvetica,8\"" << std::endl
          << "set output \"" << filename << ".png\"" << std::endl
          << "set size square" << std::endl
          << "set view equal xy" << std::endl
          << "unset xtics" << std::endl
          << "unset ytics" << std::endl
          << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle" << std::endl;
        GridOut().write_gnuplot (triangulation, f);
        f << "e" << std::endl;

        DoFTools::write_gnuplot_dof_support_point_info(f,
                                                       support_points);
        f << "e" << std::endl;

      }

      {
        DataOut<dim> data_out;
        data_out.attach_dof_handler (dh);

        std::vector<LinearAlgebra::distributed::Vector<double> > shape_functions(dh.n_dofs());
        for (unsigned int i = 0; i < dh.n_dofs(); ++i)
          {
            LinearAlgebra::distributed::Vector<double> sl(locally_owned_set, MPI_COMM_WORLD);
            sl = 0.;
            if (locally_owned_set.is_element(i))
              sl[i] = 1.0;
            cm.distribute(sl);

            LinearAlgebra::distributed::Vector<double> &s = shape_functions[i];
            s.reinit(locally_owned_set,
                     locally_relevant_set,
                     MPI_COMM_WORLD);
            s = 0.;
            s = sl;

            data_out.add_data_vector (s,
                                      std::string("N_") +
                                      dealii::Utilities::int_to_string(i));
          }

        Vector<float> subdomain (triangulation.n_active_cells());
        for (unsigned int i=0; i<subdomain.size(); ++i)
          subdomain(i) = triangulation.locally_owned_subdomain();
        data_out.add_data_vector (subdomain, "subdomain");
        data_out.build_patches ();

        const std::string filename = output_name(triangulation.locally_owned_subdomain());

        std::ofstream output (filename.c_str());
        data_out.write_vtu (output);

        if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          {
            std::vector<std::string> filenames;
            for (unsigned int i=0; i<dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
              filenames.push_back (output_name(i));

            const std::string master_name = "output.pvtu";
            std::ofstream pvtu_master (master_name.c_str());
            data_out.write_pvtu_record (pvtu_master, filenames);
          }
      }
    }

  for (unsigned int i = 0; i < locally_owned_set.n_elements(); ++i)
    {
      const unsigned int ind = locally_owned_set.nth_index_in_set(i);
      const double v = rhs[ind];
      AssertThrow(std::abs(v) < 1e-12,
                  ExcMessage(
                    "Element " + std::to_string(ind) +
                    " has an error " + std::to_string(v) +
                    " on the process " + std::to_string(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
                  ));
    }

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Ok" << std::endl;


  dh.clear();
}

int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile,/*do not print job id*/false);
  deallog.depth_console(0);
  deallog.precision(4);

  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  test<2>();

  return 0;
}

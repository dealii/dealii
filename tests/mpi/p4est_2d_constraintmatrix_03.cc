// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// check ConstraintMatrix.distribute() for a distributed mesh
// with Trilinos; manual check of the graphical output...
// Mesh: shell with random refinement

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <sstream>


template<int dim>
class FilteredDataOut : public DataOut<dim>
{
public:
  FilteredDataOut (const unsigned int subdomain_id)
    :
    subdomain_id (subdomain_id)
  {}

  virtual typename DataOut<dim>::cell_iterator
  first_cell ()
  {
    typename DataOut<dim>::active_cell_iterator
    cell = this->triangulation->begin_active();
    while ((cell != this->triangulation->end()) &&
           (cell->subdomain_id() != subdomain_id))
      ++cell;

    return cell;
  }

  virtual typename DataOut<dim>::cell_iterator
  next_cell (const typename DataOut<dim>::cell_iterator &old_cell)
  {
    if (old_cell != this->triangulation->end())
      {
        const IteratorFilters::SubdomainEqualTo
        predicate(subdomain_id);

        return
          ++(FilteredIterator
             <typename DataOut<dim>::active_cell_iterator>
             (predicate,old_cell));
      }
    else
      return old_cell;
  }

private:
  const unsigned int subdomain_id;
};




template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  const double R0      = 6371000.-2890000.;
  const double R1      = 6371000.-  35000.;


  GridGenerator::hyper_shell (tr,
                              Point<dim>(),
                              R0,
                              R1,
                              12,
                              true);
  static HyperShellBoundary<dim> boundary;
  tr.set_boundary (0, boundary);
  tr.set_boundary (1, boundary);

  tr.refine_global (1);
  for (unsigned int step=0; step<20; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tr.begin_active(),
      endc = tr.end();

      for (; cell!=endc; ++cell)
        if (Testing::rand()%42==1)
          cell->set_refine_flag ();

      tr.execute_coarsening_and_refinement ();
    }

  DoFHandler<dim> dofh(tr);

  static FESystem<dim> fe (FE_Q<dim>(1+1), dim,
                           FE_Q<dim>(1), 1);

  dofh.distribute_dofs (fe);

  IndexSet owned_set = dofh.locally_owned_dofs();

  IndexSet dof_set;
  DoFTools::extract_locally_active_dofs (dofh, dof_set);

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);
  x=2.0;

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);

  ConstraintMatrix cm(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, cm);
  std::vector<bool> velocity_mask (dim+1, true);

  velocity_mask[dim] = false;

  VectorTools::interpolate_boundary_values (dofh,
                                            0,
                                            ZeroFunction<dim>(dim+1),
                                            cm,
                                            velocity_mask);

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (1);


  VectorTools::compute_no_normal_flux_constraints (dofh, 0,
                                                   no_normal_flux_boundaries,
                                                   cm);

  cm.close ();

  cm.distribute(x);
  x_rel = x;

  std::vector<std::string> joint_solution_names (dim, "vel");
  joint_solution_names.push_back ("p");

  FilteredDataOut<dim> data_out (tr.locally_owned_subdomain());
  data_out.attach_dof_handler (dofh);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  data_component_interpretation
  (dim+1, DataComponentInterpretation::component_is_scalar);
  for (unsigned int i=0; i<dim; ++i)
    data_component_interpretation[i]
      = DataComponentInterpretation::component_is_part_of_vector;

  data_out.add_data_vector(x_rel, joint_solution_names,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);
  data_out.build_patches (4);
  const std::string filename = ("solution." +
                                Utilities::int_to_string
                                (tr.locally_owned_subdomain(), 4) +
                                ".d2");
  std::ofstream output (filename.c_str());
  data_out.write_deal_II_intermediate (output);

  TrilinosWrappers::Vector x_dub;
  x_dub.reinit(dof_set.size());
  x_dub = x_rel;

  if (myid==0)
    {
      std::ofstream file((std::string("dat.") + Utilities::int_to_string(myid)).c_str());
      file << "**** proc " << myid << std::endl;
      x_dub.print(file);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0)
    {
      cat_file((std::string("dat.") + Utilities::int_to_string(0)).c_str());
    }

  tr.set_boundary (0);
  tr.set_boundary (1);
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

}

//---------------------------------------------------------------------------
//    $Id: no_flux_constraints.cc 24175 2011-08-24 12:14:21Z kronbichler $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check that the ConstraintMatrix with hanging nodes and no-normal-flux
// constraints on an adaptively refined hyper_cube are the same independet
// of the number of CPUs

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>
#include <sstream>

template<int dim>
void test()
{
  Assert (dim == 3, ExcNotImplemented());
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  GridGenerator::hyper_cube (triangulation, -1.0, 1.0);
  triangulation.refine_global (3);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (!cell->is_ghost() && !cell->is_artificial())
      if (cell->center().norm() < 0.3)
        {
          cell->set_refine_flag();
        }

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

  if (myid == 0)
    deallog <<  "#cells = " << triangulation.n_global_active_cells()
	    << std::endl;

				// create FE_System and fill in no-normal flux
				// conditions on boundary 1 (outer)
  static const FESystem<dim> fe(FE_Q<dim> (1), dim);
  DoFHandler<dim> dofh(triangulation);
  dofh.distribute_dofs (fe);
  DoFRenumbering::hierarchical(dofh);
  
  if (myid == 0)
    deallog <<  "#dofs = " << dofh.locally_owned_dofs().size()
	    << std::endl;

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  ConstraintMatrix constraints;
  constraints.reinit(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, constraints);
  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  const unsigned int degree = 1;
  VectorTools::compute_no_normal_flux_constraints (dofh, 0,
						   no_normal_flux_boundaries,
						   constraints,
						   MappingQ<dim>(degree));
  constraints.close();

  if (myid==0)
	system("rm -rf no_flux_constraints_03/cm_?.dot");
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  { //write the constraintmatrix to a file on each cpu
	char fname[] = "no_flux_constraints_03/cm_0.dot";
	fname[26]+=myid;
	std::ofstream file(fname);
	constraints.print(file);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  sleep(1);
  if (myid==0)
  {
	//sort and merge the constraint matrices on proc 0, generate a checksum
	//and output that into the deallog
	system("cat no_flux_constraints_03/cm_?.dot|sort -n|uniq >no_flux_constraints_03/cm");
	system("md5sum no_flux_constraints_03/cm >no_flux_constraints_03/cm.check");
	{
	  std::ifstream file("no_flux_constraints_03/cm.check");
	  std::string str;
	  while (!file.eof())
	  {
		std::getline(file, str);
		deallog << str << std::endl;
	  }
	}
  }

				// print the number of constraints. since
				// processors might write info in different
				// orders, copy all numbers to root processor
  std::vector<unsigned int> n_constraints_glob (numprocs);
  unsigned int n_constraints = constraints.n_constraints();
  MPI_Gather (&n_constraints, 1, MPI_UNSIGNED,
	      &n_constraints_glob[0], 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if (myid == 0)
    for (unsigned int i=0; i<numprocs; ++i)
      deallog <<  "#constraints on " << i << ": " << n_constraints_glob[i]
	      << std::endl;

				// dummy assembly: put 1 in all components of
				// the vector
  TrilinosWrappers::MPI::Vector vector;
  vector.reinit (dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  {
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    Vector<double> local_vector (dofs_per_cell);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      local_vector(i) = 1.;
    typename DoFHandler<dim>::active_cell_iterator
      cell = dofh.begin_active(),
      endc = dofh.end();
    for (; cell!=endc; ++cell)
      if (cell->subdomain_id() == triangulation.locally_owned_subdomain())
	{
	  cell->get_dof_indices (local_dof_indices);
	  constraints.distribute_local_to_global (local_vector,
						  local_dof_indices,
						  vector);
	}
    vector.compress (Add);
  }

				// now check that no entries were generated
				// for constrained entries on the locally
				// owned range.
  const std::pair<unsigned int,unsigned int> range = vector.local_range();
  for (unsigned int i=range.first; i<range.second; ++i)
    if (constraints.is_constrained(i))
      Assert (vector(i)==0, ExcInternalError());

  if (myid==0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  {
    Utilities::System::MPI_InitFinalize mpi_init (argc, argv);
    unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);


    deallog.push(Utilities::int_to_string(myid));

    if (myid == 0)
      {
	std::ofstream logfile(output_file_for_mpi("no_flux_constraints_03").c_str());
	deallog.attach(logfile);
	deallog.depth_console(0);
	deallog.threshold_double(1.e-10);

	deallog.push("3d");
	test<3>();
	deallog.pop();
      }
    else
      test<3>();
  }

  return 0;
}

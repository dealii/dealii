//---------------------------------------------------------------------------
//    $Id$
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


// check ConstraintMatrix.distribute() for a distributed mesh
// with Trilinos; manual check of the graphical output...
// Mesh: shell with random refinement

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vectors.h>
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

    virtual typename DoFHandler<dim>::cell_iterator
    first_cell ()
      {
	typename DoFHandler<dim>::active_cell_iterator
	  cell = this->dofs->begin_active();
	while ((cell != this->dofs->end()) &&
	       (cell->subdomain_id() != subdomain_id))
	  ++cell;

	return cell;
      }

    virtual typename DoFHandler<dim>::cell_iterator
    next_cell (const typename DoFHandler<dim>::cell_iterator &old_cell)
      {
	if (old_cell != this->dofs->end())
	  {
	    const IteratorFilters::SubdomainEqualTo
	      predicate(subdomain_id);

	    return
	      ++(FilteredIterator
		 <typename DoFHandler<dim>::active_cell_iterator>
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
  for (unsigned int step=0; step<20;++step)
        {
          typename Triangulation<dim>::active_cell_iterator
            cell = tr.begin_active(),
            endc = tr.end();

            for (; cell!=endc; ++cell)
                if (std::rand()%42==1)
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
  x.compress();

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);
  x_rel = 0;
  x_rel.compress();

  ConstraintMatrix cm(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, cm);
  std::vector<bool> velocity_mask (dim+1, true);

  velocity_mask[dim] = false;

  VectorTools::interpolate_boundary_values (dofh,
					      0,
					      ZeroFunction<dim>(dim+1),
					      cm,
					      velocity_mask);

    std::set<unsigned char> no_normal_flux_boundaries;
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
    const std::string filename = ("p4est_2d_constraintmatrix_03/solution." +
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
	std::ofstream file((std::string("p4est_2d_constraintmatrix_03/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(myid)).c_str());
	file << "**** proc " << myid << std::endl;
	x_dub.print(file);
      }

    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0)
      {
	cat_file((std::string("p4est_2d_constraintmatrix_03/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(0)).c_str());
      }

  tr.set_boundary (0);
  tr.set_boundary (1);
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_constraintmatrix_03").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}

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


// ConstraintMatrix.distribute() produces a different result when using a
// Trilinos::Vector with ghost elements (e.g. owned vs. active), which is a
// bug. Now distribute() throws an Exception when called with a Vector with
// ghost elements. check that this happens.

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
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  
  GridGenerator::hyper_cube(tr);
 
  tr.refine_global (1);
  for (unsigned int step=0; step<5;++step)
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
  
    cm.close ();

    TrilinosWrappers::MPI::Vector x_test;
    x_test.reinit(x_rel);
    
    x_test=x;

    bool throwing=false;
    try
      {	
	cm.distribute(x_test);
      }
    catch (ExcMessage e)
      {
	if (myid==0)
	  deallog << "Exception: " << e.what() << std::endl;
	throwing=true;
      }
    Assert(throwing, ExcInternalError());

    cm.distribute(x);
    x_rel = x;

				     //l2_norm() not possible for ghosted vectors...
    //double a=0;//x_test.l2_norm();
    //double b=0;//x_rel.l2_norm();
    
/*    if (myid==0)
      deallog << a << " vs " << b << std::endl;
*/
/*    Assert (x_test.l2_norm() == x_rel.l2_norm(),
      ExcInternalError());
*/
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("constraint_matrix_trilinos_bug").c_str());
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

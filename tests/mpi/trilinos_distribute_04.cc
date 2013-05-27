//---------------------------------------------------------------------------
//    $Id: trilinos_distribute_02.cc 29593 2013-05-25 12:29:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// document bug in the new ConstraintMatrix::distribute():
/*

 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <sstream>



template <int dim>
void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  if (myid==0)
    {
      typename Triangulation<dim>::active_cell_iterator it = tr.begin_active();
      ++it;++it;++it;
      for (unsigned int i=0;i<5;++i)
	{
	  it->set_refine_flag();
	  ++it;	  
	}
    }
  tr.execute_coarsening_and_refinement();

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  IndexSet locally_relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh,
					   locally_relevant_set);
  
  ConstraintMatrix cm;
  cm.reinit (locally_relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, cm);

  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  VectorTools::interpolate_boundary_values (dofh,
					    dirichlet_boundary,
					    cm);
  
  cm.close();
  
  TrilinosWrappers::MPI::Vector vec (dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  for (unsigned int i=vec.local_range().first; i<vec.local_range().second; ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  deallog <<"--" << std::endl;
  solution.locally_owned_elements().print(deallog);
  deallog <<"--" << std::endl;
  locally_relevant_set.print(deallog);
  deallog <<"--" << std::endl;
  constraints.print(deallog);
  
  cm.distribute (vec);

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_distribute_04").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test<2>();
    }
  else
    test<2>();

}

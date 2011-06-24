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
// with Trilinos

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
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
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <sstream>

template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global (2);
  for (unsigned int step=0; step<15; ++step)
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

  static const FE_Q<dim> fe(1);
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

				   //x.print(std::cout);
  cm.distribute(x);
  x_rel = x;

				   //x.print(std::cout);
  
//  x_rel.print(std::cout);

  TrilinosWrappers::Vector x_dub;
  x_dub.reinit(dof_set.size());
  
  x_dub = x_rel;
  
    {
    std::ofstream file((std::string("p4est_3d_constraintmatrix_02/dat.") + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
    x_dub.print(file);
  }
    MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0)
    {
      for (unsigned int i=0;i<numproc;++i)
        {
          cat_file((std::string("p4est_3d_constraintmatrix_02/dat.") + Utilities::int_to_string(i)).c_str());
        }

    }

  

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
      std::ofstream logfile(output_file_for_mpi("p4est_3d_constraintmatrix_02").c_str());
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

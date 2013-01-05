//---------------------------------------------------------------------------
//    $Id: p4est_data_out_01.cc 25200 2012-03-02 21:11:21Z heister $
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


// distribute dofs and mgdofs for a parallel Triangulation

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
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim> dofh(tr);

    {
      std::ofstream file((JobIdentifier::base_name(__FILE__)+"/ncpu_" + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
          for (unsigned int lvl=0;lvl<tr.n_levels();++lvl)
            {
              file << "level " << lvl << ": ";
              typename DoFHandler<dim>::cell_iterator
              cell = dofh.begin(lvl),
              endc = dofh.end(lvl);

              for (;cell!=endc;++cell)
                {
                if (cell->level_subdomain_id()!=4294967294)
                  file << cell->level_subdomain_id();
                else
                  file << "-";
                }
              file << std::endl;
            }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    
  if (myid==0)
    {
      for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);++i)
        {
          cat_file((JobIdentifier::base_name(__FILE__)+"/ncpu_" + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(i)).c_str());
        }
    }      

  static const FE_DGP<dim> fe(0);
  Assert(dofh.has_active_dofs()==false, ExcInternalError());
  Assert(dofh.has_level_dofs()==false, ExcInternalError());

  dofh.distribute_dofs (fe);

  Assert(dofh.has_active_dofs()==true, ExcInternalError());
  Assert(dofh.has_level_dofs()==false, ExcInternalError());

  dofh.distribute_mg_dofs (fe);

  Assert(dofh.has_active_dofs()==true, ExcInternalError());
  Assert(dofh.has_level_dofs()==true, ExcInternalError());
  if (myid==0)
    {
      deallog << "Levels: " << tr.n_global_levels() << std::endl;
      std::cout << "Levels: " << tr.n_global_levels() << std::endl;
      
      deallog << "n_locally_owned_dofs_per_processor:" << std::endl;
      for (unsigned int i=0; i<dofh.n_locally_owned_dofs_per_processor().size(); ++i)
	deallog << dofh.n_locally_owned_dofs_per_processor()[i] << std::endl;

      deallog << "locally_owned_mg_dofs_per_processor:" << std::endl;
      for (unsigned int lvl=0;lvl<tr.n_global_levels();++lvl)
	{
	  deallog << "level " << lvl << ":" << std::endl;
	  
	  const std::vector<IndexSet> & vec = dofh.locally_owned_mg_dofs_per_processor(lvl);
      
	  for (unsigned int i=0; i<vec.size(); ++i)
	    deallog << vec[i].n_elements() << std::endl;
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

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi(JobIdentifier::base_name(__FILE__)).c_str());
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

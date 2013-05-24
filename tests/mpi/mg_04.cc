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


// distributing FE_Q multigrid DoFs fails
/*
An error occurred in line <1344> of file </scratch/deal-trunk/deal.II/source/dofs/dof_handler_policy.cc> in function
    void dealii::internal::DoFHandler::Policy::{anonymous}::set_mg_dofindices_recursively(const dealii::parallel::distributed::Triangulation<dim, spacedim>&, const typename dealii::internal::p4est::types<dim>::quadrant&, const typename dealii::DoFHandler<dim, spacedim>::level_cell_iterator&, const typename dealii::internal::p4est::types<dim>::quadrant&, unsigned int*, unsigned int) [with int dim = 2; int spacedim = 2; typename dealii::internal::p4est::types<dim>::quadrant = p4est_quadrant; typename dealii::DoFHandler<dim, spacedim>::level_cell_iterator = dealii::TriaIterator<dealii::DoFCellAccessor<dealii::DoFHandler<2>, true> >]
The violated condition was: 
    (dof_indices[i] == (DoFHandler<dim,spacedim>::invalid_dof_index)) || (dof_indices[i]==dofs[i])
The name and call sequence of the exception was:
    ExcInternalError()
Additional Information: 
(none)
*/

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
void output(parallel::distributed::Triangulation<dim> & tr)
{
  const std::string filename = ("mg_04/mesh." +
                                Utilities::int_to_string
                                (tr.locally_owned_subdomain(), 4) +
                                ".vtu");


}

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
					       Triangulation<dim>::none,
					       parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_ball (tr);
  //GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  DoFHandler<dim> dofh(tr);

  output(tr);

  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs (fe);
  dofh.distribute_mg_dofs (fe);

  {
   std::ofstream file((std::string("mg_04/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
          for (unsigned int lvl=0;lvl<tr.n_levels();++lvl)
            {
              file << "level " << lvl << ": ";
              typename DoFHandler<dim>::cell_iterator
              cell = dofh.begin(lvl),
              endc = dofh.end(lvl);

              for (;cell!=endc;++cell)
                {
              std::vector<unsigned int> dofs(fe.n_dofs_per_cell());
              cell->get_mg_dof_indices(dofs);

              for (unsigned int i=0;i<dofs.size();++i)
                if (dofs[i]==numbers::invalid_dof_index)
                  file << "- ";
                else
                  file << dofs[i] << " ";
              file << " | ";
                }
              file << std::endl;
            }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0)
    {
      for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);++i)
        {
          cat_file((std::string("mg_04/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(i)).c_str());
        }
    }      

  if (myid==0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("mg_04").c_str());
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

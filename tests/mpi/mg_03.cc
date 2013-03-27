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


// distribute_dofs on one cpu crashes with a segmentation fault.
/* valgrind says:
==7944== Invalid read of size 8
==7944==    at 0x7A60612: dealii::DoFCellAccessor<dealii::DoFHandler<2, 2> >::DoFCellAccessor(dealii::Triangulation<2, 2> const*, int, int, dealii::DoFCellAccessor<dealii::DoFHandler<2, 2> >::AccessorData const*) (dof_accessor.templates.h:3515)
==7944==    by 0x7AC3618: dealii::TriaRawIterator<dealii::DoFCellAccessor<dealii::DoFHandler<2, 2> > >::TriaRawIterator() (tria_iterator.templates.h:36)
==7944==    by 0x7AC074D: dealii::TriaIterator<dealii::DoFCellAccessor<dealii::DoFHandler<2, 2> > >::TriaIterator() (tria_iterator.templates.h:160)
==7944==    by 0x7B70D2D: dealii::TriaActiveIterator<dealii::DoFCellAccessor<dealii::DoFHandler<2, 2> > >::TriaActiveIterator() (tria_iterator.templates.h:330)
==7944==    by 0x9A434DA: dealii::internal::DoFHandler::Policy::ParallelDistributed<2, 2>::distribute_dofs(unsigned int, dealii::DoFHandler<2, 2>&) const (dof_handler_policy.cc:2054)
==7944==    by 0x99FE383: dealii::DoFHandler<2, 2>::distribute_dofs(dealii::FiniteElement<2, 2> const&, unsigned int) (dof_handler.cc:984)
==7944==    by 0x41C8E7: void test<2>() (mg_03.cc:73)
==7944==    by 0x417949: main (mg_03.cc:134)
==7944==  Address 0x0 is not stack'd, malloc'd or (recently) free'd
==7944==
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
  const std::string filename = ("mg_03/mesh." +
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

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim> dofh(tr);

  output(tr);

  static const FE_DGP<dim> fe(0);
  dofh.distribute_dofs (fe);
  dofh.distribute_mg_dofs (fe);

  {
   std::ofstream file((std::string("mg_03/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
          for (unsigned int lvl=0;lvl<tr.n_levels();++lvl)
            {
              file << "level " << lvl << ": ";
              typename DoFHandler<dim>::cell_iterator
              cell = dofh.begin(lvl),
              endc = dofh.end(lvl);

              for (;cell!=endc;++cell)
                {
              std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell());
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
          cat_file((std::string("mg_03/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(i)).c_str());
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
      std::ofstream logfile(output_file_for_mpi("mg_03").c_str());
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

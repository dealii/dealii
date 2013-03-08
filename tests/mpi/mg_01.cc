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


// check level_subdomain_id for distributed Triangulation

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>

template<int dim>
void output(parallel::distributed::Triangulation<dim> & tr)
{
    const std::string filename = ("mg_01/mesh." +
                                Utilities::int_to_string
                                (tr.locally_owned_subdomain(), 4) +
                                ".fig");
    std::ofstream stream(filename.c_str());

	GridOutFlags::XFig flags;
	flags.color_by = GridOutFlags::XFig::level_subdomain_id;
	GridOut out;
	out.set_flags(flags);

	out.write_xfig(tr, stream);
}

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
					       Triangulation<dim>::none,
					       parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim> dofh(tr);

  output(tr);

  //tr.begin_active()->set_level_subdomain_id(1+myid);

  for (unsigned int i=0;i<10;++i)
  {

      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==i)
        {
          deallog << "ID = " << i << std::endl;
          for (unsigned int lvl=0;lvl<tr.n_levels();++lvl)
            {
              typename Triangulation<dim>::cell_iterator
              cell = tr.begin(lvl),
              endc = tr.end(lvl);

              for (; cell!=endc; ++cell)
                if (cell->level_subdomain_id()!=4294967294)
                  deallog << cell->level_subdomain_id();
                else
                  deallog << "-";
              deallog << std::endl;
            }
          deallog << std::endl;
        }
  }

  if (myid==0)
    deallog << "my levels: " << tr.n_levels() << "<= global levels:" << tr.n_global_levels() << std::endl;
  

//  static const FE_DGP<dim> fe(0);
//  dofh.distribute_dofs (fe);

//  DataOut<dim> data_out;
//  data_out.attach_dof_handler (dofh);
//
//  TrilinosWrappers::MPI::Vector x;
//
//  x.reinit(dofh.locally_owned_dofs(), MPI_COMM_WORLD);
//  x=2.0;
//
//  data_out.add_data_vector (x, "x");
//  data_out.build_patches ();
//
//  if (myid==0)
//    {
//      for (unsigned int i=0; i<dofh.n_locally_owned_dofs_per_processor().size(); ++i)
//	deallog << dofh.n_locally_owned_dofs_per_processor()[i] << std::endl;
//      data_out.write_vtu (deallog.get_file_stream());
//    }
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("mg_01").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

#ifdef DEAL_II_WITH_MPI
  MPI_Finalize();
#endif
}

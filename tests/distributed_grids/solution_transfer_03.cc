//---------------------------------------------------------------------------
//    $Id: 3d_refinement_02.cc 19830 2009-10-12 18:19:51Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// test distributed solution_transfer with one processor

#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <distributed/tria.h>
#include <distributed/solution_transfer.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/intergrid_map.h>

#include <lac/petsc_vector.h>
#include <fe/fe_q.h>

#include <numerics/vectors.h>
#include <fe/mapping_q1.h>
#include <base/function.h>

#include <fstream>
#include <cstdlib>


#include <grid/grid_out.h>

template<int dim>
class MyFunction : public Function<dim>
{
  public:
    MyFunction () : Function<dim>() {};

    virtual double value (const Point<dim>   &p,
                          const unsigned int) const
      { double f=p[0]*2.0+1.0;
        if (dim>1)
          f*=p[1]*3.3-1.0;
        if (dim>2)
          f*=p[2]*5.0;
        return f;
      };
};

template<int dim>
void test(std::ostream& /*out*/)
{
  MyFunction<dim> func;
  MappingQ1<dim> mapping;
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global (3);
  DoFHandler<dim> dofh(tr);
  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs (fe);

  parallel::distributed::SolutionTransfer<dim, Vector<double> > soltrans(dofh);

  for (int i=0;i<4;++i)
    {

      for (typename Triangulation<dim>::cell_iterator
	     cell = tr.begin();
	   cell != tr.end(); ++cell)
	{
	  if (cell->has_children() && !(rand() % 3))
	    for (unsigned int c=0;c<GeometryInfo<dim>::max_children_per_cell; ++c)
	      if (!cell->child(c)->has_children())
		cell->child(c)->set_coarsen_flag();
	}
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = tr.begin_active();
	   cell != tr.end(); ++cell)
	{
	  if (!(rand() % 3))
	    {
	      cell->clear_coarsen_flag();
	      cell->set_refine_flag();
	    }

	}

      tr.prepare_coarsening_and_refinement();

      Vector<double> solution(dofh.n_dofs());
      VectorTools::interpolate (mapping,
				* static_cast<dealii::DoFHandler<dim>* >(&dofh),
				func,
				solution);

      soltrans.prepare_for_coarsening_and_refinement(solution);

      tr.execute_coarsening_and_refinement ();
      if (0)
	{
	  std::ofstream out ("grid-1.eps");
	  GridOut grid_out;
	  grid_out.write_eps (tr, out);
	}

      dofh.distribute_dofs (fe);

      Vector<double> interpolated_solution(dofh.n_dofs());
      soltrans.interpolate(interpolated_solution);

      deallog << "norm: " << interpolated_solution.l2_norm() << std::endl;
      Vector<double> difference(tr.n_global_active_cells());

      VectorTools::integrate_difference(mapping,
					dofh,
					interpolated_solution,
					func,
					difference,
					QGauss<dim>(2),
					VectorTools::L2_norm);
      deallog << "error: " << difference.l2_norm() << std::endl;
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

  std::ofstream logfile("solution_transfer_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}

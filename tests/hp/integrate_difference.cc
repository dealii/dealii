//----------------------------  integrate_difference.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  integrate_difference.cc  ---------------------------


// call VectorTools::integrate_difference with fe's distributed in the
// same random way as in hp/random


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>



template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;
  
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (4-dim);  

  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> q_collection;
  for (unsigned int i=1; i<=4; ++i)
    {
      fe_collection.push_back(FE_Q<dim> (i));
      q_collection.push_back (QGauss<dim> (i+2));
    }
  

  hp::DoFHandler<dim> dof_handler(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (rand() % fe_collection.size());
  
  dof_handler.distribute_dofs(fe_collection);

  Vector<double> vec (dof_handler.n_dofs());
  for (unsigned int i=0; i<vec.size(); ++i)
    vec(i) = i;

  Vector<float> diff (tria.n_active_cells());

  VectorTools::NormType norms[] = 
    {
          VectorTools::mean,
          VectorTools::L1_norm,
          VectorTools::L2_norm,
          VectorTools::Linfty_norm,
          VectorTools::H1_seminorm,
          VectorTools::W1p_seminorm
    };
  for (unsigned int i=0; i<sizeof(norms)/sizeof(norms[0]); ++i)
    {
      VectorTools::integrate_difference (dof_handler,
                                         vec,
                                         Functions::SquareFunction<dim>(),
                                         diff,
                                         q_collection,
                                         norms[i]);
      deallog << "i=" << i << ", diff=" << diff.l2_norm() << std::endl;
    }
}


int main ()
{
  std::ofstream logfile("integrate_difference/output");
  logfile.precision(2);
  deallog << std::setprecision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}

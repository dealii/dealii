//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// like _10, but with only a single non-primitive element. this test exists in
// order to find out why at the time of writing the test the branch
// implementing distributed meshes produced different output for the _10 test
//
// this test uses the FE_RaviartThomas as the single non-primitive element

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim>& tr,
	   const FiniteElement<dim>& fe)
{
  deallog << "FE=" << fe.get_name()
	  << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  Vector<double> fe_function(dof.n_dofs());
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    fe_function(i) = i+1;
  
  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
			   update_values | update_gradients | update_hessians);
  fe_values.reinit (dof.begin_active());

  std::vector<Tensor<2,dim> > scalar_values (quadrature.size());
  std::vector<std::vector<Tensor<2,dim> > >
    vector_values (quadrature.size(),
		   std::vector<Tensor<2,dim> >(fe.n_components()));

  fe_values.get_function_hessians (fe_function, vector_values);
  
  for (unsigned int c=0; c<fe.n_components(); ++c)
    {
      FEValuesExtractors::Scalar single_component (c);
      fe_values[single_component].get_function_hessians (fe_function,
							 scalar_values);
      deallog << "component=" << c << std::endl;
      
      for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
	{
	  deallog << scalar_values[q] << std::endl;
	  Assert ((scalar_values[q] - vector_values[q][c]).norm()
		  <= 1e-12 * scalar_values[q].norm(),
		  ExcInternalError());
	}
    }
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const HyperBallBoundary<dim> boundary;
  tr.set_boundary (0, boundary);

  FE_RaviartThomas<dim> fe(1);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("fe_values_view_10_single_01/output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}

//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// test the FEValues views and extractor classes. this test is for
// get_function_gradients for vector components and a primitive element

#include "../tests.h"
#include <base/logstream.h>
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

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

  std::vector<Tensor<2,dim> > selected_vector_values (quadrature.n_quadrature_points);
  std::vector<std::vector<Tensor<1,dim> > >
    vector_values (quadrature.n_quadrature_points,
		   std::vector<Tensor<1,dim> >(fe.n_components()));

  fe_values.get_function_gradients (fe_function, vector_values);
  
  for (unsigned int c=0; c<fe.n_components(); ++c)
				     // use a vector extractor if there
				     // are sufficiently many components
				     // left after the current component
				     // 'c'
    if (c+dim <= fe.n_components())
      {
	FEValuesExtractors::Vector vector_components (c);
	fe_values[vector_components].get_function_gradients (fe_function,
							    selected_vector_values);
	deallog << "component=" << c << std::endl;
      
	for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
	  for (unsigned int d=0; d<dim; ++d)
	    {
	      deallog << selected_vector_values[q][d] << std::endl;
	      Assert ((selected_vector_values[q][d] - vector_values[q][c+d]).norm()
		      <= 1e-12 * selected_vector_values[q][d].norm(),
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

  FESystem<dim> fe (FE_Q<dim>(1), 1,
		    FE_Q<dim>(2), 2,
		    FE_DGQ<dim>(3), dim);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("fe_values_view_12/output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-8);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}

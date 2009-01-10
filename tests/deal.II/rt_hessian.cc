//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// the Hessian of the RT element was a tensor of NaN's. This doesn't make much
// sense

#include "../tests.h"
#include <base/logstream.h>
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim>& tr,
	   const FiniteElement<dim>& fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature, update_covariant_transformation | update_hessians);

  fe_values.reinit (dof.begin_active());

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      deallog << fe_values.shape_hessian_component (0,0,0)[i][j] << std::endl;
  
				   // compare the hessian with
				   // itself. this fails if the
				   // values are NaN's which the
				   // Hessian consists of at the
				   // time this test is written
  Assert (fe_values.shape_hessian_component (0,0,0)
	  ==
	  fe_values.shape_hessian_component (0,0,0),
	  ExcInternalError());
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  test(tr, FE_RaviartThomas<dim>(1));
}


int main()
{
  std::ofstream logfile ("rt_hessian/output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-8);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}

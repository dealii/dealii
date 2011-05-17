//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// computing the Hessian of the RT element without computing values failed
// because some field wasn't properly initialized

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim>& tr,
	   const FiniteElement<dim>& fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name()
	  << std::endl;

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature, update_hessians);

					   // we used to fail here
  fe_values.reinit (dof.begin_active());

  deallog << "OK" << std::endl;
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
  std::ofstream logfile ("rt_covariant/output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-10);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}

//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// make sure FEValuesExtractors::Vector can be default-constructed but that it
// produces an invalid, unusable object

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim>& tr,
	   const FiniteElement<dim>& fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
			   update_values);
  fe_values.reinit (dof.begin_active());

  FEValuesExtractors::Vector extr; // invalid object
  try
    {
      fe_values[extr];             // invalid access
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
      goto ok;
    }

  Assert (false, ExcMessage ("No exception!?"));

  ok:
  ;
}



template<int dim>
void test()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe (FE_Q<dim>(1), dim);
  test(tr, fe);
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("fe_values_view_invalid_02/output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test<1>();
  test<2>();
  test<3>();
}

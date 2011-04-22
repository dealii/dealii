//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// like _24, but for a simpler mesh for which the output has been
// verified to be correct

#include "../tests.h"
#include <base/logstream.h>
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
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
			   update_values | update_gradients | update_q_points);
  fe_values.reinit (dof.begin_active());

				   // let the FEValues object compute the
				   // divergences at quadrature points
  std::vector<Tensor<1,dim> > divergences (quadrature.size());
  FEValuesExtractors::SymmetricTensor<2> extractor(0);
  fe_values[extractor]
    .get_function_divergences (fe_function, divergences);

				   // now do the same "by hand"
  std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
  dof.begin_active()->get_dof_indices (local_dof_indices);

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      deallog << "i=" << i << std::endl;

      for (unsigned int q=0; q<quadrature.size(); ++q)
	deallog << "  q_point=" << fe_values.quadrature_point(q) << std::endl
		<< "    value= " << fe_values[extractor].value (i,q) << std::endl
      		<< "    div= " << fe_values[extractor].divergence (i,q) << std::endl;
    }
}



template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe (FE_Q<dim>(1),
		    SymmetricTensor<2,dim>::n_independent_components);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("fe_values_view_25/output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_cube<1>();
  test_hyper_cube<2>();
  test_hyper_cube<3>();
}

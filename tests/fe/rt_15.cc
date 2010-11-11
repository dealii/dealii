//----------------------------  rt_15.cc  ---------------------------
//    rt_15.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_15.cc  ---------------------------

// Like rt_13, but use FESubfaceValues
//
// the test used to fail because of the issue with computing the
// normals using FEFaceValue, where FEFaceValue by accident uses the
// *face* mapping, not the *cell* mapping to compute the Piola
// transform (leading to a missing power of h in the determinant)

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/vector_memory.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2


std::ofstream logfile ("rt_15/output");

template<int dim>
void
test (const unsigned int degree)
{
  FE_RaviartThomas<dim> fe_rt(degree);

  deallog << "Degree=" << degree
	  << std::endl;
  
  for (double h=1; h>1./128; h/=2)
    {      
      deallog << "  h=" << h
	      << std::endl;

      Triangulation<dim> tr;
      GridGenerator::hyper_cube(tr, 0., h);

      DoFHandler<dim> dof(tr);
      dof.distribute_dofs(fe_rt);

      QTrapez<dim-1> quadrature;

      FESubfaceValues<dim> fe_values (fe_rt, quadrature, update_gradients);
      fe_values.reinit (dof.begin_active(), 0, 0);
      for (unsigned int q=0; q<quadrature.size(); ++q)
	{
	  deallog << "    Quadrature point " << q << ": ";
	  for (unsigned int i=0; i<fe_rt.dofs_per_cell; ++i)
	    {
	      deallog << '[';
	      for (unsigned int c=0; c<fe_rt.n_components(); ++c)
		deallog << fe_values.shape_grad_component(i,q,c) << ' ';
	      deallog << ']';
	    }
	  deallog << std::endl;
	}
    }
}


int
main()
{
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int i=0; i<4; ++i)
    test<2>(i);
  
  return 0;
}




//----------------------------  interpolate_rt_01.cc  ---------------------------
//    $Id: interpolate_rt_01.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  interpolate_rt_01.cc  ---------------------------


// check that VectorTools::interpolate works for RT elements correctly on
// a uniformly refined mesh for functions of degree q

#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/vectors.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

#include <fstream>
#include <vector>


template <int dim>
class F :  public Function<dim>
{
  public:
    F (const unsigned int q,
       const unsigned int n_components)
		    :
		    Function<dim>(n_components),
		    q(q)
      {}
    
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &v) const
      {
	for (unsigned int c=0; c<v.size(); ++c)
	  {
	    v(c) = 0;
	    for (unsigned int d=0; d<dim; ++d)
	      for (unsigned int i=0; i<=q; ++i)
		v(c) += (d+1)*(i+1)*std::pow (p[d], 1.*i)+c;
	  }
      }

  private:
    const unsigned int q;
};



template <int dim>
void test ()
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  for (unsigned int p=0; p<6-dim; ++p)
    {
      FE_RaviartThomas<dim> fe(p);
      DoFHandler<dim>        dof_handler(triangulation);
      dof_handler.distribute_dofs (fe);

      deallog << "n_dofs=" << dof_handler.n_dofs() << std::endl;
      
      ConstraintMatrix constraints;
      constraints.close ();

      Vector<double> interpolant (dof_handler.n_dofs());
      Vector<float>  error (triangulation.n_active_cells());
      for (unsigned int q=0; q<=p+2; ++q)
	{
					   // interpolate the function
	  VectorTools::project (dof_handler,
				constraints,
				QGauss<dim>(p+2),
				F<dim> (q, fe.n_components()),
				interpolant);
      
					   // then compute the interpolation error
	  VectorTools::integrate_difference (dof_handler,
					     interpolant,
					     F<dim> (q, fe.n_components()),
					     error,
					     QGauss<dim>(std::max(p,q)+1),
					     VectorTools::L2_norm);
	  deallog << fe.get_name() << ", P_" << q
		  << ", rel. error=" << error.l2_norm() / interpolant.l2_norm()
		  << std::endl;
	  
	  if (q<=p)
	    Assert (error.l2_norm() < 1e-12*interpolant.l2_norm(),
		    ExcInternalError());
	}
    }
}



int main ()
{
  std::ofstream logfile("interpolate_rt_01/output");
  logfile.precision (3);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}


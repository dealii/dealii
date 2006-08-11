//----------------------------  project_common.cc  ---------------------------
//    $Id: project_common.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_common.cc  ---------------------------


// common framework to check whether an element of polynomial order p can
// represent functions of order q

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
#include <fe/fe_abf.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgp_monomial.h>
#include <fe/fe_dgp_nonparametric.h>
#include <fe/fe_dgq.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_q.h>
#include <fe/fe_q_hierarchical.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

#include <fstream>
#include <vector>


template <int dim>
void test ();


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
    
    virtual double value (const Point<dim> &p,
			  const unsigned int component) const
      {
	Assert ((component == 0) && (this->n_components == 1),
		ExcInternalError());
	double val = 0;
	for (unsigned int d=0; d<dim; ++d)
	  for (unsigned int i=0; i<=q; ++i)
	    val += (d+1)*(i+1)*std::pow (p[d], 1.*i);
	return val;
      }


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
void test_no_hanging_nodes (const FiniteElement<dim> &fe,
			    const unsigned int        p)
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  DoFHandler<dim>        dof_handler(triangulation);
  dof_handler.distribute_dofs (fe);

  deallog << "n_dofs=" << dof_handler.n_dofs() << std::endl;

				   // there are no constraints here
  ConstraintMatrix constraints;
  constraints.close ();

  Vector<double> projection (dof_handler.n_dofs());
  Vector<float>  error (triangulation.n_active_cells());
  for (unsigned int q=0; q<=p+2; ++q)
    {
				       // project the function
      VectorTools::project (dof_handler,
			    constraints,
			    QGauss<dim>(p+2),
			    F<dim> (q, fe.n_components()),
			    projection);
      
				       // then compute the interpolation error
      VectorTools::integrate_difference (dof_handler,
					 projection,
					 F<dim> (q, fe.n_components()),
					 error,
					 QGauss<dim>(std::max(p,q)+1),
					 VectorTools::L2_norm);
      deallog << fe.get_name() << ", P_" << q
	      << ", rel. error=" << error.l2_norm() / projection.l2_norm()
	      << std::endl;
	  
      if (q<=p)
	Assert (error.l2_norm() <= 1e-12*projection.l2_norm(),
		ExcInternalError());
    }
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}


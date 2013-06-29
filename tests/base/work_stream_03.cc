//----------------------------  work_stream_03.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  work_stream_03.cc  ---------------------------


// Moritz originally implemented thread local scratch objects for
// WorkStream in r24748 but it led to failures in the testsuite. what
// exactly went on was a mystery and this test is a first step in
// figuring out what happens by running a simplified version of one of
// the failing tests (deal.II/project_q_01) multiple times and
// verifying that it indeed works

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <vector>


char logname[] = "work_stream_03/output";



template <int dim>
class F :  public Function<dim>
{
  public:
    F (const unsigned int q)
		    :
		    q(q)
      {}

    virtual double value (const Point<dim> &p,
			  const unsigned int component) const
      {
	Assert (component == 0, ExcInternalError());
	double val = 0;
	for (unsigned int d=0; d<dim; ++d)
	  for (unsigned int i=0; i<=q; ++i)
	    val += (d+1)*(i+1)*std::pow (p[d], 1.*i);
	return val;
      }

  private:
    const unsigned int q;
};


template <int dim>
void do_project (const Triangulation<dim> &triangulation,
		 const unsigned int        p)
{
  std::cout << "Start: " << __PRETTY_FUNCTION__ << ' ' << p << std::endl;
  FE_Q<dim> fe(p);
  DoFHandler<dim>        dof_handler(triangulation);
  dof_handler.distribute_dofs (fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler,
					   constraints);
  constraints.close ();

  Vector<double> projection (dof_handler.n_dofs());
  Vector<float>  error (triangulation.n_active_cells());
  for (unsigned int q=0; q<=p; ++q)
    {
				       // project the function
      VectorTools::project (dof_handler,
			    constraints,
			    QGauss<dim>(p+2),
			    F<dim> (q),
			    projection);
				       // just to make sure it doesn't get
				       // forgotten: handle hanging node
				       // constraints
      constraints.distribute (projection);

				       // then compute the interpolation error
      VectorTools::integrate_difference (dof_handler,
					 projection,
					 F<dim> (q),
					 error,
					 QGauss<dim>(std::max(p,q)+1),
					 VectorTools::L2_norm);
      Assert (error.l2_norm() <= 1e-12*projection.l2_norm(),
	      ExcInternalError());
    }
  std::cout << "Done: " << __PRETTY_FUNCTION__ << ' ' << p << std::endl;
}




template <int dim>
void test ()
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  Threads::TaskGroup<> g;
  for (unsigned int p=1; p<4; ++p)
    g += Threads::new_task (&do_project<dim>, triangulation, p);
  g.join_all ();
}


int main ()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}


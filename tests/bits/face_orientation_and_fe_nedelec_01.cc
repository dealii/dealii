//----------------------------  face_orientation_and_fe_nedelec_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  face_orientation_and_fe_nedelec_01.cc  ---------------------------


// make sure that we treat FE_Nedelec elements correctly if
// face_orientation==false. for nedelec elements, there is actually nothing to
// do because these elements only have dofs on edges, not on faces, but we
// should test anyway that there is nothing that actively goes wrong.

char logname[] = "face_orientation_and_fe_nedelec_01/output";


#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <lac/constraint_matrix.h>
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



DeclException1 (ExcFailedProjection,
		double,
		<< "The projection was supposed to exactly represent the "
		<< "original function, but the relative residual is "
		<< arg1);


template <int dim>
void do_project (const Triangulation<dim> &triangulation,
		 const FiniteElement<dim> &fe,
		 const unsigned int        p,
		 const unsigned int        order_difference)
{  
  DoFHandler<dim>        dof_handler(triangulation);
  dof_handler.distribute_dofs (fe);

  deallog << "n_dofs=" << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler,
					   constraints);
  constraints.close ();

  Vector<double> projection (dof_handler.n_dofs());
  Vector<float>  error (triangulation.n_active_cells());
  for (unsigned int q=0; q<=p+2-order_difference; ++q)
    {
				       // project the function
      VectorTools::project (dof_handler,
			    constraints,
			    QGauss<dim>(p+2),
			    F<dim> (q, fe.n_components()),
			    projection);
				       // just to make sure it doesn't get
				       // forgotten: handle hanging node
				       // constraints
      constraints.distribute (projection);
      
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
	  
      if (q<=p-order_difference)
	Assert (error.l2_norm() <= 1e-10*projection.l2_norm(),
		ExcFailedProjection(error.l2_norm() / projection.l2_norm()));
    }
}



// test with a 3d grid that has cells with face_orientation==false and hanging
// nodes. this trips up all sorts of pieces of code, for example there was a
// crash when computing hanging node constraints on such faces (see
// bits/face_orientation_crash), and it triggers all sorts of other
// assumptions that may be hidden in places
//
// the mesh we use is the 7 cells of the hyperball mesh in 3d, with each of
// the cells refined in turn. that then makes 7 meshes with 14 active cells
// each. this also cycles through all possibilities of coarser or finer cell
// having face_orientation==false
template <int dim>
void test_with_wrong_face_orientation (const FiniteElement<dim> &fe,
				       const unsigned int        p,
				       const unsigned int        order_difference = 0)
{
  if (dim != 3)
    return;
  
  for (unsigned int i=0; i<7; ++i)
    {
      Triangulation<dim>     triangulation;
      GridGenerator::hyper_ball (triangulation);
      typename Triangulation<dim>::active_cell_iterator
	cell = triangulation.begin_active();
      std::advance (cell, i);
      cell->set_refine_flag ();
      triangulation.execute_coarsening_and_refinement ();
  
      do_project (triangulation, fe, p, order_difference);
    }
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
}



template <int dim>
void test ()
{
  if (dim > 1)
    test_with_wrong_face_orientation (FE_Nedelec<dim>(1), 0);
}

//----------------------------  vectors_boundary_rhs.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vectors_boundary_rhs.cc  ---------------------------


// check VectorTools::create_boundary_right_hand_side


#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <lac/constraint_matrix.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <hp/fe_collection.h>
#include <hp/q_collection.h>
#include <fe/mapping_q.h>
#include <hp/mapping_collection.h>
#include <numerics/vectors.h>

#include <fstream>


template<int dim>
class MySquareFunction : public Function<dim>
{
  public:
    MySquareFunction () : Function<dim>(2) {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const
      {	return (component+1)*p.square(); }
    
    virtual void   vector_value (const Point<dim>   &p,
				 Vector<double>     &values) const
      {
	for (unsigned int d=0; d<dim; ++d)
	  values(d) = value(p,d);
      }
};


template <int dim>
void
check ()
{
  Triangulation<dim> tr;  
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

				   // create a system element composed
				   // of one Q1 and one Q2 element
  FESystem<dim> element(FE_Q<dim>(1), 1,
			FE_Q<dim>(2), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

				   // use a more complicated mapping
				   // of the domain and a quadrature
				   // formula suited to the elements
				   // we have here
  MappingQ<dim> mapping(3);

  QGauss<dim-1> quadrature(3);

  Vector<double> rhs (dof.n_dofs());
  VectorTools::create_boundary_right_hand_side (dof, quadrature,
				       MySquareFunction<dim>(),
				       rhs);
  for (unsigned int i=0; i<rhs.size(); ++i)
    deallog << rhs(i) << std::endl;
}



int main ()
{
  std::ofstream logfile ("vectors_boundary_rhs_01/output");
  deallog << std::setprecision (4);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}

//----------------------------  derivative_approximation.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  derivative_approximation.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include <base/logstream.h>
#include <base/function_lib.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/mapping_q.h>
#include <numerics/vectors.h>
#include <numerics/error_estimator.h>

#include <fstream>


template<int dim>
class MySquareFunction : public Function<dim>
{
  public:
    MySquareFunction () : Function<dim>(2) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const
      {	return (component+1)*p.square(); };
    
    virtual void   vector_value (const Point<dim>   &p,
				 Vector<double>     &values) const
      { values(0) = value(p,0);
	values(1) = value(p,1); };
};



template <int dim>
Quadrature<dim-1> &
get_q_face (Function<dim>&)
{
  static QGauss4<dim-1> q;
  return q;
};


Quadrature<0> &
get_q_face (Function<1>&)
{
  return *static_cast<Quadrature<0>*>(0);
};




template <int dim>
void
check ()
{
  CosineFunction<dim> function;
  
  Triangulation<dim> tr;  
  if (dim==2)
    GridGenerator::hyper_ball(tr);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);
  
  FE_Q<dim> element(3);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  MappingQ<dim> mapping(3);
  Quadrature<dim-1> &q_face = get_q_face(function);

  std::map<unsigned char,const Function<dim>*> neumann_bc;
  neumann_bc[0] = &function;
  
  Vector<double> v (dof.n_dofs());
  VectorTools::interpolate (mapping, dof, function, v);

  Vector<float> error (tr.n_active_cells());

  KellyErrorEstimator<dim>::estimate (mapping, dof, q_face, neumann_bc,
				      v, error);

  deallog << "Estimated error:" << std::endl;
  for (unsigned int i=0; i<error.size(); ++i)
    deallog << error(i)*100 << std::endl;
}


int main ()
{
  ofstream logfile ("error_estimator.output");
  logfile.precision (2);
  logfile.setf(ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}

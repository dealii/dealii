//----------------------------  boundaries.cc  ---------------------------
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
//----------------------------  boundaries.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */
/* Purpose: check interpolation and projection of boundary values. */



#include <base/logstream.h>
#include <base/function_lib.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/mapping_q.h>
#include <numerics/vectors.h>

#include <fstream>


template<int dim>
class MySquareFunction : public Function<dim>
{
  public:
    MySquareFunction () : Function<dim>(2) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const
      {	return 100*(component+1)*p.square()*sin(p.square()); };
    
    virtual void   vector_value (const Point<dim>   &p,
				 Vector<double>     &values) const
      { for (unsigned int d=0; d<n_components; ++d) values(d) = value(p,d); };
};


template <int dim>
const Quadrature<dim-1> &
boundary_q (const DoFHandler<dim> &)
{
  static const QGauss4<dim-1> q;
  return q;
};


const Quadrature<0> &
boundary_q (const DoFHandler<1> &)
{
  return *static_cast<const Quadrature<0>*>(0);
};


void write_map (const std::map<unsigned int,double> &bv)
{
  for (std::map<unsigned int,double>::const_iterator
	 i=bv.begin(); i!=bv.end(); ++i)
    deallog << i->first << ' ' << i->second <<std::endl;
};

      


template <int dim>
void
check ()
{
  Triangulation<dim> tr;  
  if (dim==2)
    {
      GridGenerator::hyper_ball(tr, Point<dim>(), 1);
    }
  else
    GridGenerator::hyper_cube(tr, -1./sqrt(dim),1./sqrt(dim));
  if (dim != 1)
    {
      static const HyperBallBoundary<dim> boundary;
      tr.set_boundary (0, boundary);
    };
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

				   // use a cubic mapping to make
				   // things a little more complicated
  MappingQ<dim> mapping(3);

  
				   // list of finite elements for
				   // which we want check, and
				   // associated list of boundary
				   // value functions
  std::vector<const FiniteElement<dim>*> fe_list;
  std::vector<const Function<dim>*> function_list;

				   // FE1: a system of a quadratic and
				   // a linear element
  fe_list.push_back (new FESystem<dim> (FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));
  function_list.push_back (new MySquareFunction<dim>());

				   // FE2: a linear element, to make
				   // things simple
  fe_list.push_back (new FE_Q<dim> (1));
  function_list.push_back (new SquareFunction<dim>());
  
				   // check all of them
  for (unsigned int i=0; i<fe_list.size(); ++i)
    {
      const FiniteElement<dim> &fe = *fe_list[i];
      
      DoFHandler<dim> dof(tr);
      dof.distribute_dofs(fe);

      typename FunctionMap<dim>::type function_map;
      function_map[0] = function_list[i];

				       // interpolate boundary values
      deallog << "Interpolated boundary values" << std::endl;
      std::map<unsigned int,double> interpolated_bv;
      VectorTools::interpolate_boundary_values (mapping, dof, function_map,
						interpolated_bv, std::vector<bool>());
      write_map (interpolated_bv);

				       // project boundary values
				       // presently this is not
				       // implemented for 3d
      if (dim != 3)
	{
	  deallog << "Projected boundary values" << std::endl;
	  std::map<unsigned int,double> projected_bv;
	  VectorTools::project_boundary_values (mapping, dof, function_map,
						boundary_q(dof), projected_bv);
	  write_map (projected_bv);
	};
    };
  
      
				   // delete objects now no more needed
  for (unsigned int i=0; i<fe_list.size(); ++i)
    {
      delete fe_list[i];
      delete function_list[i];
    };
}


int main ()
{
  std::ofstream logfile ("boundaries.output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
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

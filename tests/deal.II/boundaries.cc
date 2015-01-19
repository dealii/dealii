// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */
/* Purpose: check interpolation and projection of boundary values. */



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>


template<int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction () : Function<dim>(2) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return 100*(component+1)*p.square()*std::sin(p.square());
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    for (unsigned int d=0; d<this->n_components; ++d) values(d) = value(p,d);
  }
};


template <int dim>
const Quadrature<dim-1> &
boundary_q (const DoFHandler<dim> &)
{
  static const QGauss<dim-1> q(4);
  return q;
}


const Quadrature<0> &
boundary_q (const DoFHandler<1> &)
{
  static const Quadrature<0> *q = 0;
  return *q;
}


void write_map (const std::map<types::global_dof_index,double> &bv)
{
  for (std::map<types::global_dof_index,double>::const_iterator
       i=bv.begin(); i!=bv.end(); ++i)
    deallog << i->first << ' ' << i->second <<std::endl;
}




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
    GridGenerator::hyper_cube(tr, -1./std::sqrt(static_cast<double>(dim)),1./std::sqrt(static_cast<double>(dim)));
  static const HyperBallBoundary<dim> boundary;
  if (dim != 1)
    {
      tr.set_boundary (0, boundary);
    }
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
  function_list.push_back (new Functions::SquareFunction<dim>());

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
      std::map<types::global_dof_index,double> interpolated_bv;
      VectorTools::interpolate_boundary_values (mapping, dof, function_map,
                                                interpolated_bv, std::vector<bool>());
      write_map (interpolated_bv);

      // project boundary values
      // presently this is not
      // implemented for 3d
      if (dim != 3)
        {
          deallog << "Projected boundary values" << std::endl;
          std::map<types::global_dof_index,double> projected_bv;
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
  initlog();
  deallog << std::setprecision (2);
  deallog << std::fixed;

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

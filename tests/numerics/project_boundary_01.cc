// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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



/* Purpose: check projection of boundary values for really small functions,
   otherwise similar to boundaries.cc. */



#include "../tests.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>



template <int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction (const unsigned int n_components) :
    Function<dim>(n_components),
    scaling(1)
  {}

  void set_scaling (const double scaling)
  {
    this->scaling = scaling;
  }

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return 100*(component+1)*p.square()*std::sin(p.square()) * scaling;
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    for (unsigned int d=0; d<this->n_components; ++d) values(d) = value(p,d);
  }

private:
  double scaling;
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
  static const Quadrature<0> *q = nullptr;
  return *q;
}


void write_map (const std::map<types::global_dof_index,double> &bv)
{
  for (std::map<types::global_dof_index,double>::const_iterator
       i=bv.begin(); i!=bv.end(); ++i)
    // also output log of value to also display small numbers
    deallog << i->first << ' ' << i->second << " "
            << (std::abs(i->second)>0 ? std::log(std::abs(i->second)) : -10000 )
            << std::endl;
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
  GridTools::copy_boundary_to_manifold_id(tr);

  static const SphericalManifold<dim> boundary;
  if (dim != 1)
    {
      tr.set_manifold (0, boundary);
    }
  tr.refine_global (1);

  MappingQ<dim> mapping(3);


  // list of finite elements for which we want check, and associated list of
  // boundary value functions
  std::vector<const FiniteElement<dim>*> fe_list;
  std::vector<MySquareFunction<dim>*> function_list;

  // FE1: a system of a quadratic and a linear element
  fe_list.push_back (new FESystem<dim> (FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));
  function_list.push_back (new MySquareFunction<dim>(2));

  // FE2: a linear element, to make things simple
  fe_list.push_back (new FE_Q<dim> (1));
  function_list.push_back (new MySquareFunction<dim>(1));

  // test four different cases for the parameter: 1, 1e-20, 1e-170, 1e-800 (= 0)
  double factors [] = {1., 1e-40, 1e-170, 1e-800};
  for (unsigned int it=0; it<4; ++it)
    {
      deallog.push(Utilities::int_to_string(it,1));

      // check all of them
      for (unsigned int i=0; i<fe_list.size(); ++i)
        {
          function_list[i]->set_scaling(factors[it]);
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

          deallog << "Projected boundary values" << std::endl;
          std::map<types::global_dof_index,double> projected_bv;
          VectorTools::project_boundary_values (mapping, dof, function_map,
                                                boundary_q(dof), projected_bv);
          write_map (projected_bv);
        };
      deallog.pop();
    }

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
  deallog << std::setprecision (3);
  deallog << std::fixed;

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}

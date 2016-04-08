// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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


// check that the kelly error estimator returns the same result whether we
// compute all indicators at once, or on different subdomains separately



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>


template<int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction () : Function<dim>(2) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return (component+1)*p.square();
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    values(0) = value(p,0);
    values(1) = value(p,1);
  }
};



template <int dim>
Quadrature<dim-1> &
get_q_face (Function<dim> &)
{
  static QGauss<dim-1> q(4);
  return q;
}


Quadrature<0> &
get_q_face (Function<1> &)
{
  Quadrature<0> *q = 0;
  return *q;
}



template <int dim>
void make_mesh (Triangulation<dim> &tria)
{

  GridGenerator::hyper_cube(tria, -1, 1);

  // refine the mesh in a random way so as to
  // generate as many cells with
  // hanging nodes as possible
  tria.refine_global (4-dim);
  const double steps[4] = { /*d=0*/ 0, 7, 3, 3 };
  for (unsigned int i=0; i<steps[dim]; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active();
      for (unsigned int index=0; cell != tria.end(); ++cell, ++index)
        if (index % (3*dim) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement ();
    }

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active (),
  endc = tria.end ();
  for (; cell!=endc; ++cell)
    {
      unsigned int subdomain = 0;
      for (unsigned int d=0; d<dim; ++d)
        if (cell->center()(d) > 0)
          subdomain |= (1<<d);
      AssertThrow (subdomain < (1<<dim), ExcInternalError());

      cell->set_subdomain_id (subdomain);
    }
}




template <int dim>
void
check ()
{
  Functions::CosineFunction<dim> function;

  Triangulation<dim> tria;
  make_mesh (tria);

  FE_Q<dim> element(QIterated<1>(QTrapez<1>(),3));
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(element);

  MappingQ<dim> mapping(3);
  Quadrature<dim-1> &q_face = get_q_face(function);

  std::map<types::boundary_id,const Function<dim>*> neumann_bc;
  neumann_bc[0] = &function;

  Vector<double> v (dof.n_dofs());
  VectorTools::interpolate (mapping, dof, function, v);

  Vector<float> error1 (tria.n_active_cells());
  Vector<float> error2 (tria.n_active_cells());

  // compute error by looking at all cells at
  // once and output this as the base line
  // results. scale results so that they show
  // up in the output file as a reasonable
  // number
  KellyErrorEstimator<dim>::estimate (mapping, dof, q_face, neumann_bc,
                                      v, error1);
  const double scaling_factor = 500000./error1.linfty_norm();
  error1 *= scaling_factor;

  deallog << "Estimated error indicators:" << std::endl;
  for (unsigned int i=0; i<error1.size(); ++i)
    deallog << error1(i) << std::endl;

  // then do the same with different
  // subdomain ids and add up the result
  for (unsigned int subdomain=0; subdomain<(1<<dim); ++subdomain)
    {
      deallog << "Subdomain id=" << subdomain << std::endl;

      Vector<float> this_error (tria.n_active_cells());
      KellyErrorEstimator<dim>::estimate (mapping, dof, q_face, neumann_bc,
                                          v, this_error,
                                          std::vector<bool>(), 0,
                                          MultithreadInfo::n_threads(),
                                          subdomain);
      this_error *= scaling_factor;

      // copy the result into error2. since
      // every invokation of the kelly
      // estimator should only operate on
      // the cells of one subdomain,
      // whenever there is something in
      // this_error, the corresponding
      // entry in error2 should still be
      // empty
      for (unsigned int i=0; i<this_error.size(); ++i)
        {
          deallog << i << ' ' << this_error(i) << std::endl;

          Assert ((this_error(i)==0) || (error2(i)==0),
                  ExcInternalError());
          if (this_error(i) != 0)
            error2(i) = this_error(i);
        }
    }

  // now compare the results of the two
  // computations
  AssertThrow (error1 == error2, ExcInternalError());

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);

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

//----------------------------  rt_6.cc  ---------------------------
//    rt_6.cc,v 1.1 2003/06/09 15:59:07 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_6.cc  ---------------------------

// adaptation of up_and_down for RT elements

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 4


template <int dim>
Point<dim> transform (const Point<dim> p)
{
  switch (dim) 
    {
      case 1:
            return p;
      case 2:
            return Point<dim>(p(0)*(1+p(1)), p(1)*(1+p(0)));
      case 3:
            return Point<dim>(p(0)*(1+p(1))*(1+p(2)),
                              p(1)*(1+p(0))*(1+p(2)),
                              p(2)*(1+p(0))*(1+p(1)));
      default:
            Assert (false, ExcNotImplemented());
            return Point<dim>();
    };
}


template <int dim>
void check_element (const Triangulation<dim> &tr,
                    const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof_handler(tr);
  dof_handler.distribute_dofs (fe);

                                   // create a mostly arbitrary
                                   // function plus a trend on this
                                   // grid
  Vector<double> tmp(dof_handler.n_dofs());
  for (unsigned int i=0; i<tmp.size(); ++i)
    tmp(i) = (i + 13*i%17);

                                   // restrict this function to the
                                   // next coarser level and
                                   // distribute it again to the
                                   // higher level
  Vector<double> x(tmp.size());
  Vector<double> v(fe.dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin();
       cell!=dof_handler.end(); ++cell)
    if (cell->has_children() &&
        cell->child(0)->active())
      {
                                         // first make sure that what
                                         // we do is reasonable. for
                                         // this, _all_ children have
                                         // to be active, not only
                                         // some of them
        for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
          Assert (cell->child(c)->active(), ExcInternalError());

                                         // then restrict and prolong
        cell->get_interpolated_dof_values (tmp, v);
        cell->set_dof_values_by_interpolation (v, x);
      };

                                   // now x is a function on the fine
                                   // grid that is representable on
                                   // the coarse grid. so another
                                   // cycle should not alter it any
                                   // more:
  Vector<double> x2(x.size());
  for (typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin();
       cell!=dof_handler.end(); ++cell)
    if (cell->has_children() &&
        cell->child(0)->active())
      {
        cell->get_interpolated_dof_values (x, v);
        cell->set_dof_values_by_interpolation (v, x2);
      };

                                   // then check that this is so:
  x2 -= x;
  const double relative_residual = (x2.l2_norm() / x.l2_norm());

  const double threshold = 1e-6;
  deallog << ", dofs_per_cell=" << fe.dofs_per_cell
          << "; relative residual: "
          << (relative_residual<threshold ? "ok" : "botched up!")
	  << " (relative residual=" << relative_residual << ")"
          << std::endl;
  
  Assert (relative_residual < threshold*x.l2_norm(), ExcInternalError());
}


template <int dim>
void test () 
{
                                   // make a coarse triangulation as a
                                   // hypercube. if in more than 1d,
                                   // distort it so that it is no more
                                   // an affine image of the
                                   // hypercube, to make things more
                                   // difficult. then refine it twice
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  Point<dim> (*p) (Point<dim>) = &transform<dim>;
  GridTools::transform (p, tr);
  tr.refine_global (2);

                                   // now for a list of finite
                                   // elements, for which we want to
                                   // test. we happily waste tons of
                                   // memory here, but who cares...
  const FiniteElement<dim> *fe_list[]
    = 
    {
          new FE_RaviartThomas<dim>(0),
          new FE_RaviartThomas<dim>(1)
    };
  
  for (unsigned int i=0; i<sizeof(fe_list)/sizeof(fe_list[0]); ++i)
    if (fe_list[i] != 0)
      {
        deallog << dim << "d, uniform grid, fe #" << i;
        check_element (tr, *fe_list[i]);
      }
}



int
main()
{
  std::ofstream logfile ("rt_6.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  test<2>();
  
  return 0;
}




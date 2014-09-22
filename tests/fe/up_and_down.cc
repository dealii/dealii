// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2


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
    tmp(i) = i;//(i + 13*i%17);

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
        for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
          Assert (cell->child(c)->active(), ExcInternalError());

        // then restrict and prolongate
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
          << std::endl;

//TODO:[WB] Why this exception with a value different from above. Output of the error should be sufficient!
//  Assert (relative_residual < threshold*x.l2_norm(), ExcInternalError());
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
    // first for some scalar
    // elements:

    // FE_Q
    new FE_Q<dim>(1),
    new FE_Q<dim>(2),
    (dim<3 ? new FE_Q<dim>(3) : 0),
    (dim<2 ? new FE_Q<dim>(4) : 0),

    // FE_DGQ
    new FE_DGQ<dim>(0),
    new FE_DGQ<dim>(1),
    new FE_DGQ<dim>(2),
    (dim<3 ? new FE_DGQ<dim>(3) : 0),
    (dim<3 ? new FE_DGQ<dim>(4) : 0),

    // FE_DGP does not
    // presently have the
    // restriction matrices
    // implemented, so comment
    // out the following
    // tests...
//           new FE_DGP<dim>(0),
//           new FE_DGP<dim>(1),
//           new FE_DGP<dim>(2),
//           new FE_DGP<dim>(3),

    // a non-primitive FE
    (dim != 1 ? new FE_Nedelec<dim>(0) : 0),
    (dim != 1 ? new FE_Nedelec<dim>(1) : 0),

    // some composed elements
    // of increasing
    // complexity, to check the
    // logics by which the
    // matrices of the composed
    // elements are assembled
    // from those of the base
    // elements. note that some
    // of the base elements are
    // additive, some not, so
    // the result will be an
    // element that is mixed in
    // this respect
    new FESystem<dim> (FE_Q<dim>(2), 2),
    new FESystem<dim> (FE_Q<dim>(1), 2,
    FE_DGQ<dim>(2), 2),
    new FESystem<dim> (FE_Q<dim>(1), 2,
    FE_DGQ<dim>(2), 2,
    FE_DGQ<dim>(0), 1),
    new FESystem<dim> (FE_Q<dim>(1), 2,
    FESystem<dim> (FE_Q<dim>(1), 2,
    FE_DGQ<dim>(2), 2,
    FE_DGQ<dim>(2), 1), 2,
    FE_DGQ<dim>(0), 1),
    new FESystem<dim> (FE_Q<dim>(1), 2,
    FESystem<dim> (FE_Q<dim>(1), 2,
    FE_DGQ<dim>(2), 2,
    FESystem<dim>(FE_DGQ<dim>(0),
    3),
    1), 2,
    FE_DGQ<dim>(0), 1),

    // finally mixed elements,
    // with scalar and Nedelec
    // elements, to make things
    // really difficult
    (dim != 1 ?
    new FESystem<dim>(FE_Nedelec<dim>(0), 2)
    : 0),
    (dim != 1 ?
    new FESystem<dim>(FE_Nedelec<dim>(0), 2,
    FE_Q<dim> (2), 2)
    : 0),
    (dim != 1 ?
    new FESystem<dim>(FE_Nedelec<dim>(0), 2,
    FE_DGQ<dim> (2), 2,
    FESystem<dim>(FE_Nedelec<dim>(0), 2,
    FE_Q<dim> (2), 2), 2)
    : 0),
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
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}




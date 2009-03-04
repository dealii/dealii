//----------------------------  project_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_01.cc  ---------------------------

// check VectorTools::project for Vector<double> arguments


#include "../tests.h"
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <hp/dof_handler.h>
#include <lac/constraint_matrix.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>
#include <dofs/dof_accessor.h>

#include <fstream>


// define the multi-linear function x or x*y or x*y*z that we will
// subsequently project onto the ansatz space
template <int dim>
class F : public Function<dim>
{
  public:
    virtual double value (const Point<dim> &p,
                          const unsigned int = 0) const
      {
        double s = 1;
        for (unsigned int i=0; i<dim; ++i)
          s *= p[i];
        return s;
      }
};


template<int dim>
void test()
{
                                   // create 2 triangulations with the
                                   // same coarse grid, and refine
                                   // them differently
  Triangulation<dim> tria;

  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dh (tria);
  dh.distribute_dofs (fe);

  Vector<double> v(dh.n_dofs());

  ConstraintMatrix cm;
  cm.close ();
  VectorTools::project (dh, cm, QGauss<dim>(3), F<dim>(),
                        v);

  for (typename DoFHandler<dim>::active_cell_iterator cell=dh.begin_active();
       cell != dh.end(); ++cell)
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      {
                                         // check that the error is
                                         // somewhat small. it won't
                                         // be zero since we project
                                         // and do not interpolate
        Assert (std::fabs (v(cell->vertex_dof_index(i,0)) -
                           F<dim>().value (cell->vertex(i)))
                < 1e-4,
                ExcInternalError());
        deallog << cell->vertex(i) << ' ' << v(cell->vertex_dof_index(i,0))
                << std::endl;
      }
}


int main()
{
  std::ofstream logfile ("project_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

// 1d is presently not implemented  
//  test<1>();
  test<2>();
  test<3>();
}


//----------------------------  mesh_3d_25.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2008 by Timo Heister and the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_25.cc  ---------------------------


// an adaptation of the mesh_3d_22 test.  check that
// VectorTools::interpolate works for FE_Q(p) elements correctly on a mesh
// with a cell that has a wrong face rotation. when using a MappingQ(3), we
// calculated interpolation points wrongly

#include "../tests.h"
#include "../bits/mesh_3d.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>
#include <cmath>
#include <vector>

using namespace dealii;
using namespace std;

template <int dim>
class F :  public Function<dim>
{
  public:
    F (const unsigned int q) : q(q) {}

    virtual double value (const Point<dim> &p,
			  const unsigned int) const
      {
	double v=0;
	for (unsigned int d=0; d<dim; ++d)
	  for (unsigned int i=0; i<=q; ++i)
	    v += (d+1)*(i+1)*std::pow (p[d], 1.*i);
	return v;
      }

  private:
    const unsigned int q;
};



template <int dim>
void test (Triangulation<dim>& triangulation)
{
  MappingQ<3> mapping(3);
  for (unsigned int p=1; p<7-dim; ++p)
    {
      FE_Q<dim>              fe(p);
      DoFHandler<dim>        dof_handler(triangulation);
      dof_handler.distribute_dofs (fe);

      Vector<double> interpolant (dof_handler.n_dofs());
      Vector<float>  error (triangulation.n_active_cells());
      for (unsigned int q=0; q<=p+2; ++q)
	{
					   // interpolate the function
	  VectorTools::interpolate (mapping, dof_handler,
				    F<dim> (q),
				    interpolant);

					   // then compute the interpolation error
	  VectorTools::integrate_difference (mapping, dof_handler,
					     interpolant,
					     F<dim> (q),
					     error,
					     QGauss<dim>(q+2),
					     VectorTools::L2_norm);
	  deallog << fe.get_name() << ", P_" << q
		  << ", rel. error=" << error.l2_norm() / interpolant.l2_norm()
		  << std::endl;
	  if (q<=p)
	    Assert (error.l2_norm() < 1e-12*interpolant.l2_norm(),
		    ExcInternalError());
	}
    }
}




int main ()
{
  std::ofstream logfile ("mesh_3d_25/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.0e-10);

  Triangulation<3> triangulation;

  create_two_cubes_rotation(triangulation,1);
  test<3>(triangulation);
}


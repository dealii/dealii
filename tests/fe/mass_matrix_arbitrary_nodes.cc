//----------------------------  mass_matrix_arbitrary_nodes.cc  ---------------------------
//    mass_matrix_arbitrary_nodes.cc,v 1.2 2004/01/23 16:34:24 wolf Exp
//    Version:
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mass_matrix_arbitrary_nodes.cc  ---------------------------

// verify that we get a diagonal matrix when using FE_DGQArbitraryNodes with
// the same quadrature formula for integration as for construction of the FE.

#include "../tests.h"
#include <iostream>
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/numerics/matrix_tools.h>




template <int dim>
void check ()
{
  deallog << dim << 'D' << std::endl;
  for (unsigned int q_points=1; q_points<7-dim; ++q_points)
    {
      deallog << "q_points=" << q_points << std::endl;

      Triangulation<dim> tr;
      if (dim == 1)
	GridGenerator::hyper_cube (tr);
      else
      	GridGenerator::hyper_ball (tr);
      tr.refine_global (3-dim);

      FE_DGQArbitraryNodes<dim> fe ((QGauss<1>(q_points)));
      DoFHandler<dim> dh (tr);
      dh.distribute_dofs (fe);

      SparsityPattern sp (dh.n_dofs(),
			  dh.n_dofs(),
			  fe.dofs_per_cell);
      DoFTools::make_sparsity_pattern (dh, sp);
      sp.compress();

      SparseMatrix<double> mass_matrix (sp);
      MatrixTools::create_mass_matrix (dh,
				       QGauss<dim>(q_points),
				       mass_matrix);

				       // verify that the matrix is diagonal
      for (unsigned int i=0; i<dh.n_dofs(); ++i)
	{
	  deallog << mass_matrix.el(i,i) << std::endl;
	  for (unsigned int j=0; j<dh.n_dofs(); ++j)
	    if (i != j)
	      Assert (std::fabs(mass_matrix.el(i,j) / mass_matrix.el(i,i)) < 1e-14,
		      ExcInternalError());
	}
    }
}




int main ()
{
  std::ofstream logfile("mass_matrix_arbitrary_nodes/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}


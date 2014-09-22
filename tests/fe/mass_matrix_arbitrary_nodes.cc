// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}


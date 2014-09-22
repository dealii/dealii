// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// this is the root cause for solution_tranfer_01: the prolongation
// matrices for FE_DGQ<dim-1,dim> were not computed at all

#include "../tests.h"
#include <deal.II/fe/fe_dgq.h>


#include <fstream>

using namespace dealii;


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int spacedim = 2;
  const unsigned int dim = spacedim-1;

  for (unsigned int degree=0; degree<3; ++degree)
    {
      deallog << "Degree=" << degree << std::endl;
      FE_DGQ<dim,spacedim> fe(degree);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          deallog << fe.get_prolongation_matrix(0)(i,j) << std::endl;
    }

  return 0;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Test n_dofs_per-methods of FE_SimplexP  and FE_SimplexDGP.


#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

using namespace dealii;

template <int dim, int spacedim>
void
test(const FiniteElement<dim, spacedim> &fe)
{
  deallog << fe.get_name() << ": " << std::endl;

  deallog << "  n_dofs_per_vertex(): " << fe.n_dofs_per_vertex() << std::endl;
  deallog << "  n_dofs_per_line():   " << fe.n_dofs_per_line() << std::endl;

  deallog << "  n_dofs_per_quad():   ";
  for (unsigned int i = 0; i < (dim == 2 ? 1 : fe.reference_cell().n_faces());
       ++i)
    deallog << fe.n_dofs_per_quad(i) << ' ';
  deallog << std::endl;

  deallog << "  n_dofs_per_hex():    " << fe.n_dofs_per_hex() << std::endl;

  deallog << "  n_dofs_per_face():   ";
  for (unsigned int i = 0; i < fe.reference_cell().n_faces(); ++i)
    deallog << fe.n_dofs_per_face(i) << ' ';
  deallog << std::endl;

  deallog << "  n_dofs_per_cell():   " << fe.n_dofs_per_cell() << std::endl;
  deallog << "  tensor_degree():     " << fe.tensor_degree() << std::endl;

  deallog << std::endl;
}

int
main()
{
  initlog();

  test(FE_SimplexP<2>(1));
  test(FE_SimplexP<2>(2));
  test(FE_SimplexP<3>(1));
  test(FE_SimplexP<3>(2));

  test(FE_SimplexDGP<2>(1));
  test(FE_SimplexDGP<2>(2));
  test(FE_SimplexDGP<3>(1));
  test(FE_SimplexDGP<3>(2));

  test(FE_WedgeP<3>(1));
  test(FE_WedgeP<3>(2));

  test(FE_WedgeDGP<3>(1));
  test(FE_WedgeDGP<3>(2));

  test(FE_PyramidP<3>(1));

  test(FE_PyramidDGP<3>(1));
}

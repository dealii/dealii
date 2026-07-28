// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that DoFs on FESimplexP are distributed correctly on the faces in
// higher order


#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <string>
#include <vector>

#include "../tests.h"


template <int dim>
void
test_interpolation(const unsigned int        degree,
                   const Triangulation<dim> &triangulation)
{
  const FE_SimplexP<dim> fe(degree);

  deallog << "Testing interpolation " << fe.get_name() << std::endl;

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  Vector<float>  error(triangulation.n_active_cells());

  // interpolate monomials of the form x1^a*x2^b*x3^c
  // tets can interpolate a + b + c <= degree exactly
  std::vector<Tensor<1, dim>> exponents;
  for (unsigned int a = 0; a <= degree; ++a)
    for (unsigned int b = 0; b <= degree - a; ++b)
      for (unsigned int c = 0; c <= degree - a - b; ++c)
        {
          Tensor<1, dim> exponent;

          exponent[0] = a;
          if constexpr (dim > 1)
            exponent[1] = b;
          if constexpr (dim > 2)
            exponent[2] = c;

          exponents.emplace_back(exponent);
        }

  deallog << "  Relative interpolation errors:";
  for (const auto &exponent : exponents)
    {
      const Functions::Monomial<dim> test_function(exponent, fe.n_components());

      MappingFE<dim> mapping(FE_SimplexP<dim>(1));
      // interpolate the function
      VectorTools::interpolate(mapping,
                               dof_handler,
                               test_function,
                               interpolant);

      // then compute the interpolation error
      QGaussSimplex<dim> quad(4);
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        interpolant,
                                        test_function,
                                        error,
                                        quad,
                                        VectorTools::L2_norm);

      if (error.l2_norm() / interpolant.l2_norm() < 1e-12)
        deallog << " ok";

      Assert(error.l2_norm() < 1e-12 * interpolant.l2_norm(),
             ExcInternalError());
    }
  deallog << std::endl;
}



template <int dim>
void
create_tria(const unsigned int  face_no,
            unsigned int        nf,
            unsigned int        orientation,
            Triangulation<dim> &tria)
{
  Triangulation<3> dummy;
  GridGenerator::reference_cell(dummy, ReferenceCells::Tetrahedron);

  auto vertices = dummy.get_vertices();

  std::vector<CellData<3>> cells;

  {
    CellData<3> cell;
    cell.vertices    = {0, 1, 2, 3};
    cell.material_id = 0;
    cells.push_back(cell);
  }

  {
    const auto                  face = dummy.begin()->face(face_no);
    std::array<unsigned int, 3> v{
      {face->vertex_index(0), face->vertex_index(1), face->vertex_index(2)}};
    ArrayView<const unsigned int> view(v.cbegin(), v.size());
    ;
    const auto permuted =
      ReferenceCells::Triangle.permute_by_combined_orientation(view,
                                                               orientation);

    auto direction =
      cross_product_3d(vertices[permuted[1]] - vertices[permuted[0]],
                       vertices[permuted[2]] - vertices[permuted[0]]);
    direction = direction / direction.norm();

    vertices.push_back(face->center() + direction);

    CellData<3> cell;
    if (nf == 0)
      cell.vertices = {permuted[0], permuted[1], permuted[2], 4};
    else if (nf == 1)
      cell.vertices = {permuted[1], permuted[0], 4, permuted[2]};
    else if (nf == 2)
      cell.vertices = {permuted[0], 4, permuted[1], permuted[2]};
    else if (nf == 3)
      cell.vertices = {4, permuted[1], permuted[0], permuted[2]};

    cell.material_id = 1;
    cells.push_back(cell);
  }

  tria.create_triangulation(vertices, cells, {});

  const auto cell = tria.begin();

  const auto face = cell->face(face_no);

  auto ncell = tria.begin();
  ncell++;

  unsigned int nf_true = 0;
  for (; nf_true < 4; ++nf_true)
    if (ncell->face(nf_true) == face)
      break;

  deallog << "Face number " << face_no << " neighbor face number " << nf_true
          << " orientation " << int(ncell->combined_face_orientation(nf_true))
          << std::endl;
}



int
main()
{
  initlog();

  // degree, face, face orientation
  //  just test some
  // for (unsigned int f = 0; f < 4; ++f)
  // for (unsigned int nf = 0; nf < 4; ++nf)
  // for (unsigned int r :{0,1,2,3,4,5})
  for (unsigned int f : {0, 1})
    for (unsigned int nf : {2, 3})
      for (unsigned int r : {1, 2})
        {
          Triangulation<3> tria;
          create_tria<3>(f, nf, r, tria);
          for (unsigned int i = 3; i < 4; ++i)
            test_interpolation<3>(i, tria);
        }

  return 0;
}

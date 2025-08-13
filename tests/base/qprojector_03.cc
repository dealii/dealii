// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Verify that QProjector::project_to_all_subfaces() produces consistent output
// before and after replacing QProjector's element-specific code with calls to
// ReferenceCell functions.

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

template <int dim>
void
test(const ReferenceCell reference_cell)
{
  deallog.push(reference_cell.to_string());

  std::vector<Point<dim - 1>> points(dim == 1 ? 1 : 2);
  std::vector<double>         weights(dim == 1 ? 1 : 2,
                              reference_cell.face_reference_cell(0).volume() /
                                (dim == 1 ? 1.0 : 2.0));
  if constexpr (dim > 1)
    for (unsigned int d = 0; d < dim - 1; ++d)
      {
        points[0][d] = 0.05 + 0.05 * d;
        points[1][d] = 0.15 - 0.075 * d;
      }

  const Quadrature<dim - 1> face_quadrature(points, weights);
  const auto                quadratures =
    QProjector<dim>::project_to_all_subfaces(reference_cell, face_quadrature);

  std::vector<internal::SubfaceCase<dim>> subface_cases;
  if (reference_cell.is_hyper_cube())
    {
      if constexpr (dim == 1)
        {
          subface_cases.push_back(
            internal::SubfacePossibilities<1>::case_isotropic);
        }
      else if constexpr (dim == 2)
        {
          subface_cases.push_back(
            internal::SubfacePossibilities<2>::case_isotropic);
        }
      else
        {
          subface_cases.push_back(internal::SubfacePossibilities<3>::case_x);
          subface_cases.push_back(internal::SubfacePossibilities<3>::case_x1y);
          subface_cases.push_back(internal::SubfacePossibilities<3>::case_x2y);
          subface_cases.push_back(
            internal::SubfacePossibilities<3>::case_x1y2y);
          subface_cases.push_back(internal::SubfacePossibilities<3>::case_y);
          subface_cases.push_back(internal::SubfacePossibilities<3>::case_y1x);
          subface_cases.push_back(internal::SubfacePossibilities<3>::case_y2x);
          subface_cases.push_back(
            internal::SubfacePossibilities<3>::case_y1x2x);
          subface_cases.push_back(
            internal::SubfacePossibilities<3>::case_isotropic);
        }
    }
  else
    {
      subface_cases.push_back(
        internal::SubfacePossibilities<dim>::case_isotropic);
    }

  for (unsigned int face_no = 0; face_no < reference_cell.n_faces(); ++face_no)
    {
      deallog << "face_no = " << face_no << std::endl;
      for (const auto &subface_case : subface_cases)
        {
          const auto n_subfaces =
            reference_cell.is_hyper_cube() ?
              (dim == 1 ? 1 : GeometryInfo<dim>::n_subfaces(subface_case)) :
              reference_cell.face_reference_cell(face_no)
                .n_isotropic_children();
          for (types::geometric_orientation combined_orientation = 0;
               combined_orientation <
               reference_cell.n_face_orientations(face_no);
               ++combined_orientation)
            for (unsigned int subface_no = 0; subface_no < n_subfaces;
                 ++subface_no)
              {
                const auto [final_subface_no, final_refinement_case] =
                  reference_cell.equivalent_refinement_case(
                    combined_orientation, subface_case, subface_no);
                std::vector<Point<dim>> subface_vertices;
                // at the present time subfaces must have the same ReferenceCell
                // has the parent face
                for (const unsigned int &vertex_no :
                     reference_cell.face_reference_cell(face_no)
                       .vertex_indices())
                  subface_vertices.push_back(
                    reference_cell.subface_vertex_location<dim>(
                      face_no,
                      final_subface_no,
                      vertex_no,
                      final_refinement_case));
                const auto subface_bbox =
                  BoundingBox<dim>(subface_vertices).create_extended(1.0e-6);

                const auto index = QProjector<dim>::DataSetDescriptor::subface(
                  reference_cell,
                  face_no,
                  subface_no,
                  combined_orientation,
                  face_quadrature.size(),
                  subface_case);
                for (unsigned int qp_n = 0; qp_n < face_quadrature.size();
                     ++qp_n)
                  {
                    deallog << quadratures.point(qp_n + index) << std::endl;
                    // Further verify that the quadrature points are on the
                    // correct face
                    AssertThrow(subface_bbox.point_inside(
                                  quadratures.point(qp_n + index)),
                                ExcInternalError());
                  }
              }
        }
    }

  deallog.pop();
}

int
main()
{
  initlog();

  test<1>(ReferenceCells::Line);

  test<2>(ReferenceCells::Triangle);
  test<2>(ReferenceCells::Quadrilateral);

  // test<2>(ReferenceCells::Tetrahedron);
  test<3>(ReferenceCells::Hexahedron);
}

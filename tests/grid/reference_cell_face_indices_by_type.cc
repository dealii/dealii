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


// Test ReferenceCell::face_indices_by_type()

#include <deal.II/base/exceptions.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// general test case
template <int dim>
void
execute_test(const std::initializer_list<ReferenceCell> &cell_types,
             const std::initializer_list<ReferenceCell> &face_types)
{
  Triangulation<dim> triangulation;

  for (auto tested_cell_type : cell_types)
    {
      // Clear triangulation and fill with one new reference cell.
      triangulation.clear();
      GridGenerator::reference_cell(triangulation, tested_cell_type);

      // get the first (and only cell)
      auto cell = *triangulation.active_cell_iterators().begin();

      // Store the information of a face that we need for testing
      struct FaceInfo
      {
        unsigned int  face_id;
        ReferenceCell face_type;
      };
      std::list<FaceInfo> returned_face_indices;

      // Extract the face information via all possible iterators that
      // face_indices_by_type could return
      for (auto tested_face_type : face_types)
        {
          // now iterate over the faces for the 'tested_face_type' and store the
          // information
          for (auto actual_ref_face_id_for_given_type :
               cell->reference_cell().face_indices_by_type(tested_face_type))
            {
              // Store the information
              returned_face_indices.push_front(
                {actual_ref_face_id_for_given_type, tested_face_type});
            }
        }

      // Some commandline output
      deallog << "Testing " << tested_cell_type.to_string() << std::endl;
      deallog << "received data face info" << std::endl
              << "id type" << std::endl;
      for (auto face_info : returned_face_indices)
        {
          deallog << " " << face_info.face_id << " "
                  << face_info.face_type.to_string() << std::endl;
        }

      // 1. check list for length (must be exactly of n_faces())
      AssertThrow(returned_face_indices.size() == cell->n_faces(),
                  ExcMessage("Too many or not enough face indices were given"));

      // 2. Check if every face index is contained exactly once
      for (unsigned int current_index = 0; current_index < cell->n_faces();
           ++current_index)
        {
          bool found_index = false;
          // search all entries for the index (don't stop if found to check for
          // duplicates)
          for (auto face_info : returned_face_indices)
            if (face_info.face_id == current_index)
              {
                // index mustn't be found twice
                AssertThrow(
                  !found_index,
                  ExcMessage(
                    "Face index has been given on multiple occasions."));
                found_index = true;
              }

          // Check whether we found the index
          AssertThrow(found_index, ExcMessage("A face index is missing."));
        }

      // 3. Check whether face indices were filtered correctly
      for (auto face_info : returned_face_indices)
        AssertThrow(
          cell->face(face_info.face_id)->reference_cell() ==
            face_info.face_type,
          ExcMessage(
            "A face index was reported under a false face filter. This means the index of a face whose reference cell doesn't match the filter was given (eg. filter for tet but got index of face that is actually quad)."));

      // Some more output
      for (unsigned int i = 0; i < tested_cell_type.n_faces(); ++i)
        {
          deallog << "face: " << i << " of type "
                  << cell->face(i)->reference_cell().to_string() << std::endl;
        }
    }
}

void
test_1d()
{
  // All known 3d cell types
  auto cell_types = {ReferenceCells::Line};
  // All known face types
  auto face_types = {ReferenceCells::Vertex,
                     ReferenceCells::Line,
                     ReferenceCells::Triangle,
                     ReferenceCells::Quadrilateral,
                     ReferenceCells::Tetrahedron,
                     ReferenceCells::Pyramid,
                     ReferenceCells::Wedge,
                     ReferenceCells::Hexahedron};

  execute_test<1>(cell_types, face_types);
}

void
test_2d()
{
  // All known 3d cell types
  auto cell_types = {ReferenceCells::Triangle, ReferenceCells::Quadrilateral};
  // All known face types
  auto face_types = {ReferenceCells::Vertex,
                     ReferenceCells::Line,
                     ReferenceCells::Triangle,
                     ReferenceCells::Quadrilateral,
                     ReferenceCells::Tetrahedron,
                     ReferenceCells::Pyramid,
                     ReferenceCells::Wedge,
                     ReferenceCells::Hexahedron};

  execute_test<2>(cell_types, face_types);
}

void
test_3d()
{
  // All known 3d cell types
  auto cell_types = {ReferenceCells::Tetrahedron,
                     ReferenceCells::Pyramid,
                     ReferenceCells::Wedge,
                     ReferenceCells::Hexahedron};
  // All known face types
  auto face_types = {ReferenceCells::Vertex,
                     ReferenceCells::Line,
                     ReferenceCells::Triangle,
                     ReferenceCells::Quadrilateral,
                     ReferenceCells::Tetrahedron,
                     ReferenceCells::Pyramid,
                     ReferenceCells::Wedge,
                     ReferenceCells::Hexahedron};

  execute_test<3>(cell_types, face_types);
}

int
main()
{
  initlog();

  // TODO: Add test for 0D (or skip it...)
  test_1d();
  test_2d();
  test_3d();
}

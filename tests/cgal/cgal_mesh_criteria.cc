// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Check that different `cell_size` parameter gives different number of cells.

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/io.h>
#include <deal.II/cgal/utilities.h>

#include "../tests.h"

using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint = CGAL::Point_3<K>;
using namespace CGALWrappers;
using Mesh_domain =
  CGAL::Polyhedral_mesh_domain_with_features_3<K,
                                               CGAL::Surface_mesh<CGALPoint>>;
using Tr = typename CGAL::
  Mesh_triangulation_3<Mesh_domain, CGAL::Default, ConcurrencyTag>::type;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                                     Mesh_domain::Corner_index,
                                                     Mesh_domain::Curve_index>;

void
test()
{
  CGAL::Surface_mesh<CGALPoint> sm;
  C3t3                          tria;
  std::ifstream                 input(SOURCE_DIR "/input_grids/cube.off");
  input >> sm;
  AdditionalData<3> data;
  data.cell_size = .1;
  cgal_surface_mesh_to_cgal_triangulation(sm, tria, data);
  const unsigned int n_intial_facets = tria.number_of_facets_in_complex();
  tria.clear();
  data.cell_size = .5;
  cgal_surface_mesh_to_cgal_triangulation(sm, tria, data);

  Assert(
    n_intial_facets > tria.number_of_facets_in_complex(),
    ExcMessage(
      "The number of facets in the finer mesh must be greater than the number of facets in the coarse mesh."));
  deallog << std::boolalpha
          << (n_intial_facets > tria.number_of_facets_in_complex())
          << std::endl;
}

int
main()
{
  initlog();
  test();
}

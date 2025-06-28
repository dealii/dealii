// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Read a Surface_mesh from an .off file, then create a coarse mesh filled with
// tets

#include <deal.II/base/config.h>

#include <deal.II/cgal/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/io.h>

#include "../tests.h"

// Create a Surface_mesh from an .off file, then fill it with tets and
// check that we succeeded.

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
  // This test might trigger spurious floating point exception despite
  // functioning properly. Simply disable floating point exceptions again
  // (after they had been enabled int tests.h)
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  {
    const int current_fe_except = fegetexcept();
    fedisableexcept(current_fe_except);
  }
#endif

  const std::vector<std::string> fnames{SOURCE_DIR "/input_grids/cube.off",
                                        SOURCE_DIR
                                        "/input_grids/tetrahedron.off"};
  CGAL::Surface_mesh<CGALPoint>  sm;
  C3t3                           tria;

  for (const auto &name : fnames)
    {
      std::ifstream input(name);
      input >> sm;
      cgal_surface_mesh_to_cgal_triangulation(sm, tria);
      Assert(tria.is_valid(), ExcMessage("Result is not valid."));
      sm.clear(); // reset surface
      tria.clear();
      deallog << "OK" << std::endl;
    }
}

int
main()
{
  initlog();
  test();
}

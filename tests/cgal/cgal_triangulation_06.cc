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

// Read a surface mesh, make a coarse CGAL triangulation out of it, and
// translate the result into a deal.II Triangulation.

#include <deal.II/base/config.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/triangulation.h>
#include <deal.II/cgal/utilities.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/io.h>

#include "../tests.h"

using namespace CGALWrappers;

using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint = K::Point_3;

using Mesh_domain =
  CGAL::Polyhedral_mesh_domain_with_features_3<K,
                                               CGAL::Surface_mesh<CGALPoint>>;
using Tr = typename CGAL::
  Mesh_triangulation_3<Mesh_domain, CGAL::Default, ConcurrencyTag>::type;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;
using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

int
main()
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

  initlog();
  // Build a deal.II triangulation
  CGAL::Surface_mesh<CGALPoint> cgal_surface_mesh;
  {
    std::ifstream input_mesh(SOURCE_DIR "/input_grids/tetrahedron.off");
    input_mesh >> cgal_surface_mesh;
  }

  C3t3 cgal_triangulation;
  cgal_surface_mesh_to_cgal_triangulation(cgal_surface_mesh,
                                          cgal_triangulation);

  Triangulation<3> tria;
  cgal_triangulation_to_dealii_triangulation(cgal_triangulation, tria);

  deallog << "OK" << std::endl;
}

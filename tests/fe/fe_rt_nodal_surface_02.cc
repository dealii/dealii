// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check that FE_RaviartThomasNodal<2, 3> shape values are correct at
// quadrature points that are processed by the scalar remainder path of
// MappingQ (i.e., for quadrature rules whose number of points leaves a
// remainder with respect to the SIMD vectorization width; a 1-point rule is
// the simplest such case). The Piola transform divides by the volume element
// sqrt(det(J^T J)), and the remainder path used to skip filling
// volume_elements for dim != spacedim, so such quadrature points divided by
// an uninitialized value, yielding NaN/inf (or garbage) shape values. Guard
// against this by checking that the values from a 1-point rule are finite
// and agree with the values at the same point from an 8-point rule, which is
// processed entirely by the vectorized path for all common SIMD widths.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int dim = 2, spacedim = 3;

  // A single quad lying in a tilted but flat plane in 3d.
  Triangulation<dim, spacedim> tria;
  {
    std::vector<Point<spacedim>> vertices = {{0, 0, 0},
                                             {1, 0, 0},
                                             {0, 1, 1},
                                             {1, 1, 1}};
    std::vector<CellData<dim>>   cells(1);
    cells[0].vertices = {0, 1, 2, 3};
    tria.create_triangulation(vertices, cells, SubCellData());
  }

  FE_RaviartThomasNodal<dim, spacedim> fe(0);
  deallog << fe.get_name() << std::endl;

  DoFHandler<dim, spacedim> dh(tria);
  dh.distribute_dofs(fe);

  MappingQ<dim, spacedim> mapping(1);

  const Point<dim> midpoint(0.5, 0.5);

  // Evaluate the shape values at the cell midpoint through a 1-point rule
  // (processed by the scalar remainder path for any SIMD width > 1) and
  // through a rule of 8 identical points (processed entirely by the
  // vectorized path for SIMD widths 2, 4, and 8), and compare.
  const Quadrature<dim> quad_remainder(midpoint);
  const Quadrature<dim> quad_vectorized(std::vector<Point<dim>>(8, midpoint),
                                        std::vector<double>(8, 1. / 8.));

  FEValues<dim, spacedim> fev_remainder(mapping,
                                        fe,
                                        quad_remainder,
                                        update_values);
  FEValues<dim, spacedim> fev_vectorized(mapping,
                                         fe,
                                         quad_vectorized,
                                         update_values);
  fev_remainder.reinit(tria.begin_active());
  fev_vectorized.reinit(tria.begin_active());

  const FEValuesExtractors::Vector field(0);

  bool   all_finite     = true;
  double max_difference = 0.0;
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      const Tensor<1, spacedim> v = fev_remainder[field].value(i, 0);
      for (unsigned int d = 0; d < spacedim; ++d)
        if (!std::isfinite(v[d]))
          all_finite = false;

      for (unsigned int q = 0; q < quad_vectorized.size(); ++q)
        max_difference =
          std::max(max_difference,
                   (v - fev_vectorized[field].value(i, q)).norm());
    }

  deallog << "  all shape values finite: " << (all_finite ? "yes" : "no")
          << std::endl;
  deallog << "  remainder and vectorized paths agree (diff < 1e-12): "
          << (max_difference < 1e-12 ? "yes" : "no") << std::endl;
}

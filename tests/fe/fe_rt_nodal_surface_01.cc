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


// Check FE_RaviartThomasNodal<2, 3> on a (flat, tilted) surface mesh embedded
// in 3d. The element represents a 3d vector field that must be tangential to
// the surface, so the shape function values dotted with the surface normal
// have to vanish. Also check that the element reports spacedim components.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../tests.h"


void
test(const unsigned int degree)
{
  const unsigned int dim = 2, spacedim = 3;

  // A single quad lying in a tilted but flat plane in 3d, so that the mapping
  // is affine and the surface Piola transform is exact.
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

  FE_RaviartThomasNodal<dim, spacedim> fe(degree);
  deallog << fe.get_name() << std::endl;
  deallog << "  n_components  = " << fe.n_components() << std::endl;
  deallog << "  dofs_per_cell = " << fe.dofs_per_cell << std::endl;

  DoFHandler<dim, spacedim> dh(tria);
  dh.distribute_dofs(fe);

  MappingQ<dim, spacedim> mapping(1);
  QGauss<dim>             quad(degree + 2);
  FEValues<dim, spacedim> fev(mapping,
                              fe,
                              quad,
                              update_values | update_normal_vectors);
  fev.reinit(tria.begin_active());

  const FEValuesExtractors::Vector field(0);

  double max_normal_component = 0.0;
  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      const Tensor<1, spacedim> n = fev.normal_vector(q);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        {
          const Tensor<1, spacedim> v = fev[field].value(i, q);
          max_normal_component =
            std::max(max_normal_component, std::abs(v * n));
        }
    }

  // The exact value is zero; print a tolerance-based verdict so the test is
  // portable across platforms.
  deallog << "  tangential (max|shape.n| < 1e-10): "
          << (max_normal_component < 1e-10 ? "yes" : "no") << std::endl;
}


int
main()
{
  initlog();

  test(0);
  test(1);
}

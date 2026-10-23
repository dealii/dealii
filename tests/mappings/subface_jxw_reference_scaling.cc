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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

void
check(const Mapping<3>       &mapping,
      const RefinementCase<3> refinement,
      const double            reference_ratio,
      const unsigned int      n_subfaces)
{
  Triangulation<3> tria;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            std::vector<unsigned int>{1, 1, 2},
                                            Point<3>(0, 0, 0),
                                            Point<3>(1, 1, 2));
  // This trilinear map is represented exactly by MappingQ(1). The top face
  // is nonplanar and has surface Jacobian sqrt(1+x*x+y*y).
  GridTools::transform(
    [](const Point<3> &p) { return Point<3>(p[0], p[1], p[2] + p[0] * p[1]); },
    tria);
  // Refine only the upper neighbor; FESubfaceValues requires the lower
  // coarse cell to remain active while its top face has hanging nodes.
  std::next(tria.begin_active())->set_refine_flag(refinement);
  tria.execute_coarsening_and_refinement();
  const auto coarse = tria.begin_active();
  AssertThrow(!coarse->has_children(), ExcInternalError());
  AssertThrow(!coarse->at_boundary(5), ExcInternalError());
  AssertThrow(coarse->neighbor(5)->has_children(), ExcInternalError());
  AssertThrow(coarse->face(5)->n_children() == n_subfaces, ExcInternalError());

  const FE_Q<3>       fe(1);
  const QGauss<2>     quadrature(8);
  FESubfaceValues<3>  subface(mapping,
                             fe,
                             quadrature,
                             update_JxW_values | update_quadrature_points |
                               update_normal_vectors);
  std::vector<double> areas(n_subfaces, 0.);
  for (unsigned int s = 0; s < n_subfaces; ++s)
    {
      subface.reinit(coarse, 5, s);
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          const auto    &unit = quadrature.point(q);
          const double   x    = .5 * (unit[0] + s % 2);
          const double   y = n_subfaces == 2 ? unit[1] : .5 * (unit[1] + s / 2);
          const Point<3> expected_point(x, y, 1. + x * y);
          AssertThrow(subface.quadrature_point(q).distance(expected_point) <
                        1e-12,
                      ExcInternalError());
          Tensor<1, 3> normal;
          normal[0]             = -y;
          normal[1]             = -x;
          normal[2]             = 1.;
          const double jacobian = normal.norm();
          normal /= jacobian;
          AssertThrow(
            std::abs(subface.JxW(q) -
                     reference_ratio * quadrature.weight(q) * jacobian) < 1e-12,
            ExcMessage("Subface JxW must use the reference measure ratio"));
          AssertThrow((subface.normal_vector(q) - normal).norm() < 1e-12,
                      ExcInternalError());
          areas[s] += subface.JxW(q);
        }
    }
  // The factor is not the ratio of physical areas: these children have
  // the same reference area but different physical surface areas.
  AssertThrow(std::abs(areas.front() - areas.back()) > 1e-2,
              ExcInternalError());
  deallog << n_subfaces
          << " subfaces: analytic JxW and normals OK; physical areas unequal"
          << std::endl;
}

int
main()
{
  initlog();
  const MappingQ<3>        q1(1);
  const MappingManifold<3> manifold;
  for (const Mapping<3> *mapping : {static_cast<const Mapping<3> *>(&q1),
                                    static_cast<const Mapping<3> *>(&manifold)})
    {
      deallog << (mapping == &q1 ? "MappingQ" : "MappingManifold") << std::endl;
      check(*mapping, RefinementCase<3>::isotropic_refinement, .25, 4);
      check(*mapping, RefinementCase<3>::cut_x, .5, 2);
      check(*mapping, RefinementCase<3>::cut_xy, .25, 4);
    }
}

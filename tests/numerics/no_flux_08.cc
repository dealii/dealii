// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// we used to get this crash:
//
// An error occurred in line <4646> of file
// </w/heister/deal-trunk/deal.II/include/deal.II/numerics/vectors.templates.h>
// in function
//     static void dealii::VectorTools::compute_no_normal_flux_constraints(const
//     DoFHandlerType<dim, spacedim>&, unsigned int, const
//     std::set<types::boundary_id>&, dealii::AffineConstraints<double>&, const
//     dealii::Mapping<dim, spacedim>&) [with int dim = 3, DoFHandlerType =
//     dealii::DoFHandler, int spacedim = 3]
// The violated condition was:
//     std::fabs (tangent.norm()-1) < 1e-12
// The name and call sequence of the exception was:
//     ExcInternalError()
//
// quarter_hyper_shell works, even though it is a very similar mesh.
//
// this was fixed with r24044 together with the no_flux_07 test that
// reduces it to its essence

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
check()
{
  Triangulation<dim> tr;
  GridGenerator::half_hyper_shell(tr, Point<dim>(), 0.5, 1, 0);
  tr.reset_manifold(0);

  AffineConstraints<double> cm;
  MappingQ<dim>             mapping(1);

  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dofh(tr);

  dofh.distribute_dofs(fe);

  const std::set<types::boundary_id> no_normal_flux_boundaries = {0};
  VectorTools::compute_no_normal_flux_constraints(
    dofh, 0, no_normal_flux_boundaries, cm, mapping);

  cm.print(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(4);
  deallog.get_file_stream().setf(std::ios::fixed);

  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}

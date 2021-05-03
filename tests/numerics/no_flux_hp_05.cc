// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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



#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr, const hp::FECollection<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  DoFRenumbering::component_wise(dof);

  for (const unsigned int i : GeometryInfo<dim>::face_indices())
    {
      deallog << "FE=" << fe[0].get_name() << ", case=" << i << std::endl;

      std::set<types::boundary_id> boundary_ids;
      for (unsigned int j = 0; j <= i; ++j)
        boundary_ids.insert(j);

      AffineConstraints<double> cm;
      VectorTools::compute_no_normal_flux_constraints(dof, 0, boundary_ids, cm);

      cm.print(deallog.get_file_stream());
    }
}


template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const SphericalManifold<dim> boundary;
  tr.set_manifold(0, boundary);

  tr.refine_global(1);

  for (unsigned int degree = 1; degree < 4; ++degree)
    {
      hp::FECollection<dim> fe(
        FESystem<dim>(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), degree)),
                      dim,
                      FE_DGP<dim>(degree + 1),
                      1));
      test(tr, fe);
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


// Check how things go when there are faces with odd face
// orientations. If the quadrature points were ordered differently on
// the two sides of an interface, i.e., if each FEFaceValues used the
// coordinate system of the cell it's on rather than the common
// coordinate system of the face, then one would get different values
// for the same shape function on the q'th quadrature point from both
// sides. We can check that by using a continuous element and
// computing the jump -- it ought to be zero at every point.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



template <int dim>
void
test()
{
  for (unsigned int twist = 0; twist < 4; ++twist)
    {
      Triangulation<dim> triangulation;
      GridGenerator::moebius(triangulation, 7, twist, 1.0, 0.2);

      DoFHandler<dim> dofh(triangulation);
      FE_Q<dim>       fe(1);
      dofh.distribute_dofs(fe);

      MappingQ<dim> mapping(1);
      UpdateFlags   update_flags = update_values;

      FEInterfaceValues<dim> fiv(mapping,
                                 fe,
                                 QGauss<dim - 1>(fe.degree + 1),
                                 update_flags);


      for (auto cell : dofh.active_cell_iterators())
        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          if (cell->at_boundary(f) == false)
            {
              fiv.reinit(cell,
                         f,
                         numbers::invalid_unsigned_int,
                         cell->neighbor(f),
                         cell->neighbor_of_neighbor(f),
                         numbers::invalid_unsigned_int);

              Assert(fiv.get_fe_face_values(0).get_cell() == cell,
                     ExcInternalError());
              Assert(fiv.get_fe_face_values(1).get_cell() == cell->neighbor(f),
                     ExcInternalError());

              // We have a continuous element, so the set of DoFs that
              // are active on the interface is twice the number of
              // DoFs per cell (namely, the two adjacent cells) minus
              // the common set of DoFs on the face itself.
              AssertDimension(fiv.n_current_interface_dofs(),
                              2 * fe.dofs_per_cell - fe.dofs_per_face);
              Assert(!fiv.at_boundary(), ExcInternalError());

              for (unsigned int q = 0; q < fiv.get_quadrature().size(); ++q)
                for (unsigned int i = 0; i < fiv.n_current_interface_dofs();
                     ++i)
                  Assert(std::abs(fiv.jump(i, q)) <= 1e-12, ExcInternalError());
            }

      deallog << "Twist " << twist << "*pi/2 seems to work just fine..."
              << std::endl;
    }
}



int
main()
{
  initlog();
  test<3>();
}

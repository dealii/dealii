// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// check serialization for CellData<structdim>

#include <deal.II/grid/tria.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"


template <int structdim>
void
test()
{
  dealii::CellData<structdim> t1;

  for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
    t1.vertices[i] = i;
  t1.material_id = 0;
  t1.boundary_id = 1;
  t1.manifold_id = 2;

  {
    dealii::CellData<structdim> t2 = t1;
    verify(t1, t2);
  }

  // test vertices
  for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
    {
      dealii::CellData<structdim> t2 = t1;
      t2.vertices[i]++;
      verify(t1, t2);
    }

  // test material_id
  {
    dealii::CellData<structdim> t2 = t1;
    t2.material_id++;
    verify(t1, t2);
  }

  // test boundary_id
  {
    dealii::CellData<structdim> t2 = t1;
    t2.boundary_id++;
    verify(t1, t2);
  }

  // test manifold_id
  {
    dealii::CellData<structdim> t2 = t1;
    t2.manifold_id++;
    verify(t1, t2);
  }
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}

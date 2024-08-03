// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2024 by the deal.II authors
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


// Test basic properties of FECouplingValues

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_coupling_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "../test_grids.h"


template <typename FCV>
void
inspect_fcv(const FCV                                  &fcv,
            const std::vector<types::global_dof_index> &v1,
            const std::vector<types::global_dof_index> &v2)
{
  deallog << "n_coupling_dofs: " << fcv.n_coupling_dofs() << std::endl;

  auto indices = fcv.get_coupling_dof_indices(v1, v2);
  Assert(indices.size() == fcv.n_coupling_dofs(), ExcInternalError());

  deallog << "coupling_dof_indices: ";
  for (const auto &i : indices)
    deallog << i << ' ';
  deallog << std::endl;

  unsigned int idx = 0;
  for (const auto &v : indices)
    {
      const auto pair = fcv.coupling_dof_to_dof_indices(idx);

      deallog << "  index " << idx << " global_dof_index: " << v
              << ", index 1: " << static_cast<int>(pair.first)
              << ", index 2: " << static_cast<int>(pair.second) << std::endl;
      ++idx;
    }

  deallog << "n_quadrature_points: " << fcv.n_quadrature_points() << std::endl;

  for (const auto &q : fcv.quadrature_point_indices())
    {
      const auto &[id1, id2] = fcv.coupling_quadrature_to_quadrature_indices(q);
      const auto &[q1, q2]   = fcv.quadrature_point(q);
      deallog << "  q index " << q << ", index 1: " << id1
              << ", index 2: " << id2 << ", q1: " << q1 << ", q2: " << q2
              << std::endl;
    }

  deallog << std::endl;
}



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  FE_Q<dim>       fe(1);
  dofh.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags = update_quadrature_points | update_JxW_values;

  FEValues<dim> fv1(mapping, fe, QGauss<dim>(fe.degree + 1), update_flags);
  FEValues<dim> fv2(mapping, fe, QGauss<dim>(fe.degree + 1), update_flags);

  FEValues<dim> fv3(mapping, fe, QTrapezoid<dim>(), update_flags);
  FEValues<dim> fv4(mapping, fe, QTrapezoid<dim>(), update_flags);

  FEFaceValues<dim> fv5(mapping, fe, QTrapezoid<dim - 1>(), update_flags);
  FEFaceValues<dim> fv6(mapping, fe, QTrapezoid<dim - 1>(), update_flags);

  std::vector<types::global_dof_index> v1(fe.dofs_per_cell);
  std::vector<types::global_dof_index> v2(fe.dofs_per_cell);

  auto cell1 = dofh.begin();
  fv1.reinit(cell1);
  cell1->get_dof_indices(v1);

  unsigned int boundary_face = numbers::invalid_unsigned_int;
  unsigned int internal_face = numbers::invalid_unsigned_int;

  for (const unsigned int f : cell1->face_indices())
    {
      if (cell1->at_boundary(f))
        {
          boundary_face = f;
        }
      else if (!cell1->at_boundary(f))
        {
          internal_face = f;
        }
      if (boundary_face != numbers::invalid_unsigned_int &&
          internal_face != numbers::invalid_unsigned_int)
        break;
    }

  auto cell2 = cell1->neighbor(internal_face);
  fv2.reinit(cell2);
  cell2->get_dof_indices(v2);

  {
    deallog << "** coupling between cell1 and cell2 (tensor) **" << std::endl;
    FECouplingValues<dim> fcv(fv1, fv2, DoFCouplingType::contiguous);
    inspect_fcv(fcv, v1, v2);
  }
  {
    deallog << "** coupling between cell1 and cell2 (unrolled) **" << std::endl;
    FECouplingValues<dim> fcv(fv1,
                              fv2,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::unrolled);
    inspect_fcv(fcv, v1, v2);
  }
  {
    deallog << "** coupling between cell1 and cell2 Trapez (overlapping) **"
            << std::endl;
    fv3.reinit(cell1);
    fv4.reinit(cell2);
    FECouplingValues<dim> fcv(fv3,
                              fv4,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::overlapping);
    inspect_fcv(fcv, v1, v2);
  }
  {
    deallog
      << "** coupling between cell1 and boundary face, trapez (overlapping) **"
      << std::endl;
    fv3.reinit(cell1);
    fv5.reinit(cell1, boundary_face);
    FECouplingValues<dim> fcv(fv3,
                              fv5,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::overlapping);
    inspect_fcv(fcv, v1, v1);
  }
  {
    deallog
      << "** coupling between cell1 and cell1->internal_face, trapez (overlapping) **"
      << std::endl;
    fv3.reinit(cell1);
    fv5.reinit(cell1, internal_face);
    FECouplingValues<dim> fcv(fv3,
                              fv5,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::overlapping);
    inspect_fcv(fcv, v1, v2);
  }
  {
    deallog
      << "** coupling between cell1 and cell2->internal_face, trapez (overlapping) **"
      << std::endl;
    fv3.reinit(cell1);
    fv5.reinit(cell2, cell1->neighbor_of_neighbor(internal_face));
    FECouplingValues<dim> fcv(fv3,
                              fv5,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::overlapping);
    inspect_fcv(fcv, v1, v2);
  }
  {
    deallog
      << "** coupling between cell1 and cell2->other_face, trapez (overlapping) **"
      << std::endl;
    fv3.reinit(cell1);
    fv5.reinit(cell2, internal_face);
    FECouplingValues<dim> fcv(fv3,
                              fv5,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::overlapping);
    inspect_fcv(fcv, v1, v2);
  }
  {
    deallog
      << "** coupling between cell1->face(0) and cell1->face(2), trapez (tensor) **"
      << std::endl;
    fv5.reinit(cell1, 0);
    fv6.reinit(cell1, 2);
    FECouplingValues<dim> fcv(fv5,
                              fv6,
                              DoFCouplingType::contiguous,
                              QuadratureCouplingType::tensor_product);
    inspect_fcv(fcv, v1, v1);
  }
}

int
main()
{
  initlog();
  test<2>();
}

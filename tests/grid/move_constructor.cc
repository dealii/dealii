// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

template <int dim>
void
print_tria_info(const Triangulation<dim>& tria)
{
  const bool manifold_0_is_flat
    = dynamic_cast<const FlatManifold<dim>*>(&tria.get_manifold(0)) != nullptr;
  deallog << (tria.n_active_cells() != 0) << ", " << (tria.n_active_hexs() != 0)
          << ", " << (tria.n_active_quads() != 0) << ", "
          << (tria.n_active_lines() != 0) << ", " << (tria.n_levels() != 0)
          << ", " << (tria.n_vertices() != 0) << ", "
          << (tria.get_periodic_face_map().size() != 0) << ", "
          << manifold_0_is_flat << std::endl;
}

template <int dim>
void
test_hyper_cube()
{
  deallog << "Dimension: " << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1.0, 1.0);
  tria.refine_global(2);

  print_tria_info(tria);

  Triangulation<dim> new_tria = std::move(tria);
  print_tria_info(new_tria);
  print_tria_info(tria);

  tria = std::move(new_tria);
  print_tria_info(new_tria);
  print_tria_info(tria);
}

template <int dim>
void
test_hyper_shell()
{
  deallog << "Dimension: " << dim << std::endl;

  const SphericalManifold<dim> boundary;

  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);

  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, boundary);

  tria.refine_global(3);
  print_tria_info(tria);

  Triangulation<dim> new_tria = std::move(tria);
  print_tria_info(new_tria);
  print_tria_info(tria);

  tria = std::move(new_tria);
  print_tria_info(new_tria);
  print_tria_info(tria);
}

template <int dim>
void
test_periodic_cube()
{
  deallog << "Dimension: " << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1.0, 1.0, true);

  typedef GridTools::PeriodicFacePair<
    typename Triangulation<dim>::cell_iterator>
                                  periodic_face_pair;
  std::vector<periodic_face_pair> periodicity_vector;
  GridTools::collect_periodic_faces(tria, 0, 1, 0, periodicity_vector);
  tria.add_periodicity(periodicity_vector);
  tria.refine_global(3);

  print_tria_info(tria);

  Triangulation<dim> new_tria = std::move(tria);
  print_tria_info(new_tria);
  print_tria_info(tria);

  tria = std::move(new_tria);
  print_tria_info(new_tria);
  print_tria_info(tria);
}

int
main()
{
  initlog();

  deallog << "Hyper-cube:" << std::endl;
  test_hyper_cube<2>();
  test_hyper_cube<3>();

  deallog << std::endl << "Hyper-shell:" << std::endl;
  test_hyper_shell<2>();
  test_hyper_shell<3>();

  deallog << std::endl << "Periodic hyper-cube:" << std::endl;
  test_periodic_cube<2>();
  test_periodic_cube<3>();

  return 0;
}

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



// adapted from mesh_3d_14, but check that quadrature points are
// correct for anisotropic refinement

#include "../tests.h"
#include "../grid/mesh_3d.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>


// declare these global in order to reduce time needed upon construction of
// these objects which is considerable
FE_DGQ<3> fe(1);
QGauss<2> quadrature(3);
MappingQ<3> mapping(2);
FEFaceValues<3> fe_face_values1 (mapping, fe, quadrature,
                                 update_quadrature_points | update_JxW_values |
                                 update_normal_vectors);
FESubfaceValues<3> fe_face_values2 (mapping, fe, quadrature,
                                    update_quadrature_points | update_JxW_values |
                                    update_normal_vectors);

void check_this (Triangulation<3> &tria)
{

  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
  for (; cell!=dof_handler.end(); ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<3>::faces_per_cell;
         ++face_no)
      if (!cell->at_boundary(face_no)
          &&
          cell->face(face_no)->has_children())
        for (unsigned int subface_no=0;
             subface_no<cell->face(face_no)->number_of_children();
             ++subface_no)
          {
            unsigned int neighbor_neighbor = cell->neighbor_face_no(face_no);

            const DoFHandler<3>::active_cell_iterator neighbor_child
              = cell->neighbor_child_on_subface(face_no, subface_no);

            fe_face_values1.reinit (neighbor_child, neighbor_neighbor);
            fe_face_values2.reinit (cell, face_no, subface_no);

            for (unsigned int q=0; q<quadrature.size(); ++q)
              {
                AssertThrow ((fe_face_values1.quadrature_point(q)-
                              fe_face_values2.quadrature_point(q)).norm_square()
                             < 1e-20,
                             ExcInternalError());

                AssertThrow (std::fabs(fe_face_values1.JxW(q)-
                                       fe_face_values2.JxW(q)) < 1e-14,
                             ExcInternalError());
                AssertThrow ((fe_face_values1.normal_vector(q) +
                              fe_face_values2.normal_vector(q)).norm_square()
                             < 1e-20,
                             ExcInternalError());
              }
          }
}


// perform the usual check, i.e. first refine a single cell (anisotropically),
// then global refinement, then global coarsening. For each step, check, that
// quadrature points both on faces and neighboring subfaces match.
void check (Triangulation<3> &tria_org)
{
  for (unsigned int c=0; c<tria_org.n_active_cells(); ++c)
    for (unsigned int i=1; i<8; ++i)
      {
        Triangulation<3> tria;
        tria.copy_triangulation(tria_org);
        Triangulation<3>::active_cell_iterator cell=tria.begin_active();
        for (unsigned int j=0; j<c; ++j)
          ++cell;
        cell->set_refine_flag (RefinementCase<3>(i));
        tria.execute_coarsening_and_refinement ();

        deallog << "Initial check, cell = "<<c<<", ref_case = " <<i<< std::endl;
        check_this (tria);

        for (unsigned int r=0; r<2; ++r)
          {
            tria.refine_global (1);
            deallog << "Check " << r << ", " << tria.n_active_cells() << " active cells"<< std::endl;
            check_this (tria);
          }

        coarsen_global (tria);
        deallog << "Check " << 0 << ", " << tria.n_active_cells() << " active cells"<< std::endl;
        check_this (tria);

        tria.refine_global (1);
        deallog << "Check " << 1 << ", " << tria.n_active_cells() << " active cells"<< std::endl;
        check_this (tria);
      }
}


// perform an additional check: simulate an isotropic refinement of a given cell
// via several anisotropic refinements. Then, perform the usual checks. This
// went wrong at some time, so check that it works now.
void check2 (Triangulation<3> &orig_tria)
{
  for (unsigned int i=0; i<orig_tria.n_active_cells(); ++i)
    {
      Triangulation<3> tria;
      tria.copy_triangulation(orig_tria);
      Triangulation<3>::cell_iterator cell = tria.begin_active(),
                                      endc = tria.end();
      for (unsigned int j=0; j<i; ++j)
        ++cell;
      cell->set_refine_flag(RefinementCase<3>::cut_z);
      tria.execute_coarsening_and_refinement ();

      cell=tria.begin();
      for (; cell!=endc; ++cell)
        if (cell->refinement_case()==RefinementCase<3>::cut_z)
          {
            cell->child(0)->set_refine_flag(RefinementCase<3>::cut_xy);
            cell->child(1)->set_refine_flag(RefinementCase<3>::cut_xy);
          }
      tria.execute_coarsening_and_refinement ();

      deallog << "2 -> Initial check, cell " << i << ", " << tria.n_active_cells() << " active cells" << std::endl;
      check_this (tria);

      for (unsigned int r=0; r<2; ++r)
        {
          tria.refine_global (1);
          deallog << "2 -> Check " << r << ", " << tria.n_active_cells() << " active cells" << std::endl;
          check_this (tria);
          deallog << "           ... done." << std::endl;
        }

      coarsen_global (tria);
      deallog << "Check " << 0 << ", " << tria.n_active_cells() << " active cells"<< std::endl;
      check_this (tria);

      tria.refine_global (1);
      deallog << "Check " << 1 << ", " << tria.n_active_cells() << " active cells"<< std::endl;
      check_this (tria);
    }
}



int main ()
{
  initlog();

  {
    Triangulation<3> coarse_grid(Triangulation<3>::allow_anisotropic_smoothing);
    GridGenerator::hyper_cube (coarse_grid);
    coarse_grid.refine_global(1);
    check2 (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_two_cubes (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    coarse_grid.reset_manifold(0);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_two_cubes_rotation (coarse_grid,1);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_two_cubes_rotation (coarse_grid,2);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_two_cubes_rotation (coarse_grid,3);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::moebius(coarse_grid, 5, 0, 10.0, 2.0);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::moebius(coarse_grid, 5, 1, 10.0, 2.0);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::moebius(coarse_grid, 5, 2, 10.0, 2.0);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::moebius(coarse_grid, 5, 3, 10.0, 2.0);
    check (coarse_grid);
  }

}




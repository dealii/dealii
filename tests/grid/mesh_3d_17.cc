// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// we used to have an abort in an assertion in
// TriaCellAccessor::neighbor_child_on_subface. it turned out that the
// assertion was wrong. make sure it works now.

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>



void check (Triangulation<3> &tria)
{
  deallog << "Coarse cell 0 vertices:" << std::endl;
  for (unsigned int i=0; i<8; ++i)
    deallog << ' ' << tria.begin_active()->vertex_index(i);
  deallog << std::endl;

  deallog << "Coarse cell 1 vertices:" << std::endl;
  for (unsigned int i=0; i<8; ++i)
    deallog << ' ' << (++tria.begin_active())->vertex_index(i);
  deallog << std::endl;


  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();

  const unsigned int dim=3;
  Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active(),
  endc = tria.end();

  // loop over all interior faces of active
  // cells with children. there should be
  // exactly one, namely the middle face of
  // the big cell (the middle faces of the
  // active cells on the refined side don't
  // have any children)
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
    for (unsigned int face_no=0;
         face_no<GeometryInfo<dim>::faces_per_cell;
         ++face_no)
      {
        const Triangulation<dim>::face_iterator
        face=cell->face(face_no);

        if (!face->at_boundary() && face->has_children())
          {
            deallog << "Interior face with children:" << std::endl;
            deallog << "  cell=" << cell << std::endl;
            deallog << "  face_no=" << face_no << std::endl;
            deallog << "  vertices="
                    << face->vertex_index(0) << ' '
                    << face->vertex_index(1) << ' '
                    << face->vertex_index(2) << ' '
                    << face->vertex_index(3) << std::endl;
            deallog << "  face orientation="
                    << (cell->face_orientation(face_no) ? "true" : "false")
                    << std::endl;

            const Triangulation<dim>::cell_iterator neighbor
              = cell->neighbor(face_no);

            const unsigned int nb_of_nb
              = cell->neighbor_of_neighbor(face_no);

            // check that
            // neighbor_of_neighbor works:
            Assert (neighbor->neighbor(nb_of_nb) == cell,
                    ExcInternalError());

            const Triangulation<dim>::face_iterator neighbor_face
              = neighbor->face(nb_of_nb);

            deallog << "  neighbor face vertices="
                    << neighbor_face->vertex_index(0) << ' '
                    << neighbor_face->vertex_index(1) << ' '
                    << neighbor_face->vertex_index(2) << ' '
                    << neighbor_face->vertex_index(3)
                    << std::endl;
            deallog << "  neighbor->face_orientation(nb_of_nb)="
                    << (neighbor->face_orientation(nb_of_nb) ?
                        "true" : "false")
                    << std::endl;


            deallog << "  checking subfaces:" << std::endl;
            for (unsigned int subface_no=0;
                 subface_no<GeometryInfo<dim>::max_children_per_face;
                 ++subface_no)
              {
                deallog << "    subface_no=" << subface_no << std::endl;
                deallog << "      center=" << face->child(subface_no)->center()
                        << std::endl;

                // we used to abort in the
                // following call, due to a
                // wrong Assert condition:
                cell->neighbor_child_on_subface(face_no, subface_no);
              }
          }
      }
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<3> coarse_grid;
  GridIn<3> in;
  in.attach_triangulation(coarse_grid);
  std::ifstream ucd_file(SOURCE_DIR "/two_cubes.inp");
  in.read_ucd(ucd_file);
  ucd_file.close();

  check (coarse_grid);
}


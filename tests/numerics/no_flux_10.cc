// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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



// test by Jennifer Worthen -- we get an error of the following kind:
// --------------------------------------------------------
// An error occurred in line <304> of file <.../constraint_matrix.cc> in function
//     void dealii::ConstraintMatrix::close()
// The violated condition was:
//     dof_index != line->line
// The name and call sequence of the exception was:
//     ExcMessage ("Cycle in constraints detected!")
// Additional Information:
// Cycle in constraints detected!
// --------------------------------------------------------



#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>


void
colorize_sixty_deg_hyper_shell(Triangulation<3> &tria,
                               const Point<3> &center,
                               const double inner_radius,
                               const double outer_radius)
{

  //    if (tria.n_cells() != 4)
  //      AssertThrow (false, ExcNotImplemented());

  double middle = (outer_radius-inner_radius)/2e0 + inner_radius;
  double eps = 1e-3*middle;
  Triangulation<3>::cell_iterator cell = tria.begin();

  for (; cell!=tria.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      {
        if (!cell->face(f)->at_boundary())
          continue;

        double radius = cell->face(f)->center().norm() - center.norm();
        if (std::fabs(cell->face(f)->center()(2) - sqrt(3.)*cell->face(f)->center()(0)) < eps ) // z = sqrt(3)x set boundary 2
          {
            cell->face(f)->set_boundary_id(2);
            for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
              if (cell->face(f)->line(j)->at_boundary())
                if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                  cell->face(f)->line(j)->set_boundary_id(2);
          }
        else if (std::fabs(cell->face(f)->center()(2) + sqrt(3.)*cell->face(f)->center()(0)) < eps) // z = -sqrt(3)x set boundary 3
          {
            cell->face(f)->set_boundary_id(3);
            for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
              if (cell->face(f)->line(j)->at_boundary())
                if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                  cell->face(f)->line(j)->set_boundary_id(3);
          }
        else if (std::fabs(cell->face(f)->center()(2) - sqrt(3.)*cell->face(f)->center()(1)) < eps ) // z = sqrt(3)y set boundary 4
          {
            cell->face(f)->set_boundary_id(4);
            for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
              if (cell->face(f)->line(j)->at_boundary())
                if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                  cell->face(f)->line(j)->set_boundary_id(4);
          }
        else if (std::fabs(cell->face(f)->center()(2) + sqrt(3.)*cell->face(f)->center()(1)) < eps ) // z = -sqrt(3)y set boundary 5
          {
            cell->face(f)->set_boundary_id(5);
            for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
              if (cell->face(f)->line(j)->at_boundary())
                if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                  cell->face(f)->line(j)->set_boundary_id(5);
          }
        else if (radius < middle) // inner radius set boundary 0
          {
            cell->face(f)->set_boundary_id(0);
            for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
              if (cell->face(f)->line(j)->at_boundary())
                if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) < eps)
                  cell->face(f)->line(j)->set_boundary_id(0);
          }
        else if (radius > middle) // outer radius set boundary 1
          {
            cell->face(f)->set_boundary_id(1);
            for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
              if (cell->face(f)->line(j)->at_boundary())
                if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) < eps)
                  cell->face(f)->line(j)->set_boundary_id(1);
          }
        else
          AssertThrow (false, ExcInternalError());
      }
}

void sixty_deg_hyper_shell (Triangulation<3> &tria,
                            const Point<3> &center,
                            const double inner_radius,
                            const double outer_radius)
{
  const double r0 = inner_radius;
  const double r1 = outer_radius;

  std::vector<Point<3> > vertices;

  vertices.push_back (center+Point<3>( 1.0/sqrt(5.)*r0, 1.0/sqrt(5.)*r0, sqrt(3./5.)*r0));   //8 -> 0
  vertices.push_back (center+Point<3>( 1.0/sqrt(5.)*r1, 1.0/sqrt(5.)*r1, sqrt(3./5.)*r1));   //9 -> 1
  vertices.push_back (center+Point<3>( 1.0/sqrt(5.)*r0, -1.0/sqrt(5.)*r0, sqrt(3./5.)*r0));  //10 -> 2
  vertices.push_back (center+Point<3>( 1.0/sqrt(5.)*r1, -1.0/sqrt(5.)*r1, sqrt(3./5.)*r1));  //11 -> 3
  vertices.push_back (center+Point<3>( -1.0/sqrt(5.)*r0, 1.0/sqrt(5.)*r0, sqrt(3./5.)*r0));  //14 -> 4
  vertices.push_back (center+Point<3>( -1.0/sqrt(5.)*r1, 1.0/sqrt(5.)*r1, sqrt(3./5.)*r1));  //15 -> 5
  vertices.push_back (center+Point<3>( -1.0/sqrt(5.)*r0, -1.0/sqrt(5.)*r0, sqrt(3./5.)*r0)); //16 -> 6
  vertices.push_back (center+Point<3>( -1.0/sqrt(5.)*r1, -1.0/sqrt(5.)*r1, sqrt(3./5.)*r1)); //17 -> 7

  const int cell_vertices[1][8] =
  {
    {6, 2, 4, 0, 7, 3, 5, 1},
  };

  std::vector<CellData<3> > cells(1);

  for (unsigned int i=0; i<1; ++i)
    {
      for (unsigned int j=0; j<8; ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  tria.create_triangulation ( vertices, cells, SubCellData());      // no boundary information

  colorize_sixty_deg_hyper_shell(tria, center, inner_radius, outer_radius);
}


template <int dim>
void run()
{
  Triangulation<dim> triangulation;
  FESystem<dim> fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler (triangulation);
  ConstraintMatrix constraints;

  sixty_deg_hyper_shell (triangulation,
                         Point<dim>(),
                         0.5,
                         1.0);
  GridTools::copy_boundary_to_manifold_id(triangulation);

  static SphericalManifold<dim> boundary((Point<dim>()));
  triangulation.set_manifold (0, boundary);

  // write out the mesh. may not be
  // strictly necessary here, but is
  // a good test for functionality
  // that may not actually be tested
  // anywhere else
  MappingQ<3> m(4);
  GridOut go;
  GridOutFlags::Gnuplot gof;
  gof.n_boundary_face_points = 6;
  go.set_flags (gof);
  go.write_gnuplot (triangulation, deallog.get_file_stream(), &m);

  dof_handler.distribute_dofs (fe);

  for (unsigned int f=0; f<6; ++f)
    deallog << "Face=" << f << ", boundary_id="
            << (int)triangulation.begin_active()->face(f)->boundary_id() << std::endl;

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  no_normal_flux_boundaries.insert (2);
  VectorTools::compute_no_normal_flux_constraints
  (dof_handler, 0,
   no_normal_flux_boundaries,
   constraints);

  constraints.print(deallog.get_file_stream());

  deallog.get_file_stream() << std::flush;
  constraints.close();

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  logfile.precision (4);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);

  run<3> ();
}

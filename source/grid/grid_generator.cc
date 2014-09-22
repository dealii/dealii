// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/distributed/tria.h>

#include <iostream>
#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN


namespace GridGenerator
{
  namespace
  {
    // Corner points of the cube [-1,1]^3
    const Point<3> hexahedron[8] =
    {
      Point<3>(-1,-1,-1),
      Point<3>(+1,-1,-1),
      Point<3>(-1,+1,-1),
      Point<3>(+1,+1,-1),
      Point<3>(-1,-1,+1),
      Point<3>(+1,-1,+1),
      Point<3>(-1,+1,+1),
      Point<3>(+1,+1,+1)
    };

    // Octahedron inscribed in the cube
    // [-1,1]^3
    const Point<3> octahedron[6] =
    {
      Point<3>(-1, 0, 0),
      Point<3>( 1, 0, 0),
      Point<3>( 0,-1, 0),
      Point<3>( 0, 1, 0),
      Point<3>( 0, 0,-1),
      Point<3>( 0, 0, 1)
    };


    /**
     * Perform the action specified by the @p colorize flag of the
     * hyper_rectangle() function of this class.
     */
    template <int dim, int spacedim>
    void
    colorize_hyper_rectangle (Triangulation<dim,spacedim> &tria)
    {
      // there is nothing to do in 1d
      if (dim > 1)
        {
          // there is only one cell, so
          // simple task
          const typename Triangulation<dim,spacedim>::cell_iterator
          cell = tria.begin();
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            cell->face(f)->set_boundary_indicator (f);
        }
    }



    void
    colorize_subdivided_hyper_rectangle (Triangulation<1> &tria,
                                         const Point<1> &,
                                         const Point<1> &,
                                         const double)
    {
      for (Triangulation<1>::cell_iterator cell = tria.begin();
           cell != tria.end(); ++cell)
        if (cell->center()(0) > 0)
          cell->set_material_id(1);
      // boundary indicators are set to
      // 0 (left) and 1 (right) by default.
    }



    template <int dim>
    void
    colorize_subdivided_hyper_rectangle (Triangulation<dim> &tria,
                                         const Point<dim>   &p1,
                                         const Point<dim>   &p2,
                                         const double        epsilon)
    {

      // run through all faces and check
      // if one of their center coordinates matches
      // one of the corner points. Comparisons
      // are made using an epsilon which
      // should be smaller than the smallest cell
      // diameter.

      typename Triangulation<dim>::face_iterator face = tria.begin_face(),
                                                 endface = tria.end_face();
      for (; face!=endface; ++face)
        {
          if (face->boundary_indicator() == 0)
            {
              const Point<dim> center (face->center());
              if (std::abs(center(0)-p1[0]) < epsilon)
                face->set_boundary_indicator(0);
              else if (std::abs(center(0) - p2[0]) < epsilon)
                face->set_boundary_indicator(1);
              else if (dim > 1 && std::abs(center(1) - p1[1]) < epsilon)
                face->set_boundary_indicator(2);
              else if (dim > 1 && std::abs(center(1) - p2[1]) < epsilon)
                face->set_boundary_indicator(3);
              else if (dim > 2 && std::abs(center(2) - p1[2]) < epsilon)
                face->set_boundary_indicator(4);
              else if (dim > 2 && std::abs(center(2) - p2[2]) < epsilon)
                face->set_boundary_indicator(5);
              else
                // triangulation says it
                // is on the boundary,
                // but we could not find
                // on which boundary.
                Assert (false, ExcInternalError());

            }
        }
      for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
           cell != tria.end(); ++cell)
        {
          char id = 0;
          for (unsigned int d=0; d<dim; ++d)
            if (cell->center()(d) > 0) id += 1 << d;
          cell->set_material_id(id);
        }
    }


    /**
     * Assign boundary number zero to the inner shell boundary and 1 to the
     * outer.
     */
    void
    colorize_hyper_shell (Triangulation<2> &tria,
                          const Point<2> &,
                          const double,
                          const double)
    {
      // In spite of receiving geometrical
      // data, we do this only based on
      // topology.

      // For the mesh based on  cube,
      // this is highly irregular
      for (Triangulation<2>::cell_iterator cell = tria.begin ();
           cell != tria.end (); ++cell)
        {
          Assert(cell->face(2)->at_boundary(), ExcInternalError());
          cell->face (2)->set_boundary_indicator (1);
        }
    }


    /**
     * Assign boundary number zero to the inner shell boundary and 1 to the
     * outer.
     */
    void
    colorize_hyper_shell (Triangulation<3> &tria,
                          const Point<3> &,
                          const double,
                          const double)
    {
      // the following uses a good amount
      // of knowledge about the
      // orientation of cells. this is
      // probably not good style...
      if (tria.n_cells() == 6)
        {
          Triangulation<3>::cell_iterator cell = tria.begin();

          cell->face(4)->set_boundary_indicator(1);
          Assert (cell->face(4)->at_boundary(), ExcInternalError());

          (++cell)->face(2)->set_boundary_indicator(1);
          Assert (cell->face(2)->at_boundary(), ExcInternalError());

          (++cell)->face(2)->set_boundary_indicator(1);
          Assert (cell->face(2)->at_boundary(), ExcInternalError());

          (++cell)->face(0)->set_boundary_indicator(1);
          Assert (cell->face(0)->at_boundary(), ExcInternalError());

          (++cell)->face(2)->set_boundary_indicator(1);
          Assert (cell->face(2)->at_boundary(), ExcInternalError());

          (++cell)->face(0)->set_boundary_indicator(1);
          Assert (cell->face(0)->at_boundary(), ExcInternalError());
        }
      else if (tria.n_cells() == 12)
        {
          // again use some internal
          // knowledge
          for (Triangulation<3>::cell_iterator cell = tria.begin();
               cell != tria.end(); ++cell)
            {
              Assert (cell->face(5)->at_boundary(), ExcInternalError());
              cell->face(5)->set_boundary_indicator(1);
            }
        }
      else if (tria.n_cells() == 96)
        {
          // the 96-cell hypershell is
          // based on a once refined
          // 12-cell mesh. consequently,
          // since the outer faces all
          // are face_no==5 above, so
          // they are here (unless they
          // are in the interior). Use
          // this to assign boundary
          // indicators, but also make
          // sure that we encounter
          // exactly 48 such faces
          unsigned int count = 0;
          for (Triangulation<3>::cell_iterator cell = tria.begin();
               cell != tria.end(); ++cell)
            if (cell->face(5)->at_boundary())
              {
                cell->face(5)->set_boundary_indicator(1);
                ++count;
              }
          Assert (count == 48, ExcInternalError());
        }
      else
        Assert (false, ExcNotImplemented());
    }



    /**
     * Assign boundary number zero the inner shell boundary, one to the outer
     * shell boundary, two to the face with x=0, three to the face with y=0,
     * four to the face with z=0.
     */
    void
    colorize_quarter_hyper_shell(Triangulation<3> &tria,
                                 const Point<3> &center,
                                 const double inner_radius,
                                 const double outer_radius)
    {
      if (tria.n_cells() != 3)
        AssertThrow (false, ExcNotImplemented());

      double middle = (outer_radius-inner_radius)/2e0 + inner_radius;
      double eps = 1e-3*middle;
      Triangulation<3>::cell_iterator cell = tria.begin();

      for (; cell!=tria.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
          {
            if (!cell->face(f)->at_boundary())
              continue;

            double radius = cell->face(f)->center().norm() - center.norm();
            if (std::fabs(cell->face(f)->center()(0)) < eps ) // x = 0 set boundary 2
              {
                cell->face(f)->set_boundary_indicator(2);
                for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
                  if (cell->face(f)->line(j)->at_boundary())
                    if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                      cell->face(f)->line(j)->set_boundary_indicator(2);
              }
            else if (std::fabs(cell->face(f)->center()(1)) < eps) // y = 0 set boundary 3
              {
                cell->face(f)->set_boundary_indicator(3);
                for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
                  if (cell->face(f)->line(j)->at_boundary())
                    if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                      cell->face(f)->line(j)->set_boundary_indicator(3);
              }
            else if (std::fabs(cell->face(f)->center()(2)) < eps ) // z = 0 set boundary 4
              {
                cell->face(f)->set_boundary_indicator(4);
                for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
                  if (cell->face(f)->line(j)->at_boundary())
                    if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) > eps)
                      cell->face(f)->line(j)->set_boundary_indicator(4);
              }
            else if (radius < middle) // inner radius set boundary 0
              {
                cell->face(f)->set_boundary_indicator(0);
                for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
                  if (cell->face(f)->line(j)->at_boundary())
                    if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) < eps)
                      cell->face(f)->line(j)->set_boundary_indicator(0);
              }
            else if (radius > middle) // outer radius set boundary 1
              {
                cell->face(f)->set_boundary_indicator(1);
                for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
                  if (cell->face(f)->line(j)->at_boundary())
                    if (std::fabs(cell->face(f)->line(j)->vertex(0).norm() - cell->face(f)->line(j)->vertex(1).norm()) < eps)
                      cell->face(f)->line(j)->set_boundary_indicator(1);
              }
            else
              AssertThrow (false, ExcInternalError());
          }
    }

  }


  template <int dim, int spacedim>
  void
  hyper_rectangle (Triangulation<dim,spacedim> &tria,
                   const Point<spacedim>   &p_1,
                   const Point<spacedim>   &p_2,
                   const bool          colorize)
  {
    // First, normalize input such that
    // p1 is lower in all coordinate directions.
    Point<spacedim> p1(p_1);
    Point<spacedim> p2(p_2);

    for (unsigned int i=0; i<spacedim; ++i)
      if (p1(i) > p2(i))
        std::swap (p1(i), p2(i));

    std::vector<Point<spacedim> > vertices (GeometryInfo<dim>::vertices_per_cell);
    switch (dim)
      {
      case 1:
        vertices[0] = p1;
        vertices[1] = p2;
        break;
      case 2:
        vertices[0] = vertices[1] = p1;
        vertices[2] = vertices[3] = p2;

        vertices[1](0) = p2(0);
        vertices[2](0) = p1(0);
        break;
      case 3:
        vertices[0] = vertices[1] = vertices[2] = vertices[3] = p1;
        vertices[4] = vertices[5] = vertices[6] = vertices[7] = p2;

        vertices[1](0) = p2(0);
        vertices[2](1) = p2(1);
        vertices[3](0) = p2(0);
        vertices[3](1) = p2(1);

        vertices[4](0) = p1(0);
        vertices[4](1) = p1(1);
        vertices[5](1) = p1(1);
        vertices[6](0) = p1(0);

        break;
      default:
        Assert (false, ExcNotImplemented ());
      }

    // Prepare cell data
    std::vector<CellData<dim> > cells (1);
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      cells[0].vertices[i] = i;
    cells[0].material_id = 0;

    tria.create_triangulation (vertices, cells, SubCellData());

    // Assign boundary indicators
    if (colorize)
      colorize_hyper_rectangle (tria);
  }


  template <int dim, int spacedim>
  void hyper_cube (Triangulation<dim,spacedim> &tria,
                   const double                 left,
                   const double                 right,
                   const bool                   colorize)
  {
    Assert (left < right,
            ExcMessage ("Invalid left-to-right bounds of hypercube"));

    Point<spacedim> p1;
    Point<spacedim> p2;

    p1(spacedim-1) = 0;
    p2(spacedim-1) = 0;

    for (unsigned int i=0; i<dim; ++i)
      {
        p1(i) = left;
        p2(i) = right;
      }
    hyper_rectangle (tria, p1, p2, colorize);
  }



  void
  moebius (Triangulation<3>  &tria,
           const unsigned int      n_cells,
           const unsigned int   n_rotations,
           const double         R,
           const double         r)
  {
    const unsigned int dim=3;
    Assert (n_cells>4, ExcMessage("More than 4 cells are needed to create a moebius grid."));
    Assert (r>0 && R>0, ExcMessage("Outer and inner radius must be positive."));
    Assert (R>r, ExcMessage("Outer radius must be greater than inner radius."));


    std::vector<Point<dim> > vertices (4*n_cells);
    double beta_step=n_rotations*numbers::PI/2.0/n_cells;
    double alpha_step=2.0*numbers::PI/n_cells;

    for (unsigned int i=0; i<n_cells; ++i)
      for (unsigned int j=0; j<4; ++j)
        {
          vertices[4*i+j][0]=R*std::cos(i*alpha_step)+r*std::cos(i*beta_step+j*numbers::PI/2.0)*std::cos(i*alpha_step);
          vertices[4*i+j][1]=R*std::sin(i*alpha_step)+r*std::cos(i*beta_step+j*numbers::PI/2.0)*std::sin(i*alpha_step);
          vertices[4*i+j][2]=r*std::sin(i*beta_step+j*numbers::PI/2.0);
        }

    unsigned int offset=0;

    std::vector<CellData<dim> > cells (n_cells);
    for (unsigned int i=0; i<n_cells; ++i)
      {
        for (unsigned int j=0; j<2; ++j)
          {
            cells[i].vertices[0+4*j]=offset+0+4*j;
            cells[i].vertices[1+4*j]=offset+3+4*j;
            cells[i].vertices[2+4*j]=offset+2+4*j;
            cells[i].vertices[3+4*j]=offset+1+4*j;
          }
        offset+=4;
        cells[i].material_id=0;
      }

    // now correct the last four vertices
    cells[n_cells-1].vertices[4]=(0+n_rotations)%4;
    cells[n_cells-1].vertices[5]=(3+n_rotations)%4;
    cells[n_cells-1].vertices[6]=(2+n_rotations)%4;
    cells[n_cells-1].vertices[7]=(1+n_rotations)%4;

    GridReordering<dim>::invert_all_cells_of_negative_grid(vertices,cells);
    tria.create_triangulation_compatibility (vertices, cells, SubCellData());
  }



  void
  torus (Triangulation<2,3>  &tria,
         const double         R,
         const double         r)
  {
    Assert (R>r, ExcMessage("Outer radius must be greater than inner radius."));

    const unsigned int dim=2;
    const unsigned int spacedim=3;
    std::vector<Point<spacedim> > vertices (16);

    vertices[0]=Point<spacedim>(R-r,0,0);
    vertices[1]=Point<spacedim>(R,-r,0);
    vertices[2]=Point<spacedim>(R+r,0,0);
    vertices[3]=Point<spacedim>(R, r,0);
    vertices[4]=Point<spacedim>(0,0,R-r);
    vertices[5]=Point<spacedim>(0,-r,R);
    vertices[6]=Point<spacedim>(0,0,R+r);
    vertices[7]=Point<spacedim>(0,r,R);
    vertices[8]=Point<spacedim>(-(R-r),0,0);
    vertices[9]=Point<spacedim>(-R,-r,0);
    vertices[10]=Point<spacedim>(-(R+r),0,0);
    vertices[11]=Point<spacedim>(-R, r,0);
    vertices[12]=Point<spacedim>(0,0,-(R-r));
    vertices[13]=Point<spacedim>(0,-r,-R);
    vertices[14]=Point<spacedim>(0,0,-(R+r));
    vertices[15]=Point<spacedim>(0,r,-R);

    std::vector<CellData<dim> > cells (16);
    //Right Hand Orientation
    cells[0].vertices[0] =  0;
    cells[0].vertices[1] =  4;
    cells[0].vertices[2] =  7;
    cells[0].vertices[3] =  3;
    cells[0].material_id = 0;

    cells[1].vertices[0] =  1;
    cells[1].vertices[1] =  5;
    cells[1].vertices[2] =  4;
    cells[1].vertices[3] =  0;
    cells[1].material_id = 0;

    cells[2].vertices[0] =  2;
    cells[2].vertices[1] =  6;
    cells[2].vertices[2] =  5;
    cells[2].vertices[3] =  1;
    cells[2].material_id = 0;

    cells[3].vertices[0] =  3;
    cells[3].vertices[1] =  7;
    cells[3].vertices[2] =  6;
    cells[3].vertices[3] =  2;
    cells[3].material_id = 0;

    cells[4].vertices[0] =  4;
    cells[4].vertices[1] =  8;
    cells[4].vertices[2] =  11;
    cells[4].vertices[3] =  7;
    cells[4].material_id = 0;

    cells[5].vertices[0] =  5;
    cells[5].vertices[1] =  9;
    cells[5].vertices[2] =  8;
    cells[5].vertices[3] =  4;
    cells[5].material_id = 0;

    cells[6].vertices[0] =  6;
    cells[6].vertices[1] =  10;
    cells[6].vertices[2] =  9;
    cells[6].vertices[3] =  5;
    cells[6].material_id = 0;

    cells[7].vertices[0] =  7;
    cells[7].vertices[1] =  11;
    cells[7].vertices[2] =  10;
    cells[7].vertices[3] =  6;
    cells[7].material_id = 0;

    cells[8].vertices[0] =  8;
    cells[8].vertices[1] =  12;
    cells[8].vertices[2] =  15;
    cells[8].vertices[3] =  11;
    cells[8].material_id = 0;

    cells[9].vertices[0] =  9;
    cells[9].vertices[1] =  13;
    cells[9].vertices[2] =  12;
    cells[9].vertices[3] =  8;
    cells[9].material_id = 0;

    cells[10].vertices[0] =  10;
    cells[10].vertices[1] =  14;
    cells[10].vertices[2] =  13;
    cells[10].vertices[3] =  9;
    cells[10].material_id = 0;

    cells[11].vertices[0] =  11;
    cells[11].vertices[1] =  15;
    cells[11].vertices[2] =  14;
    cells[11].vertices[3] =  10;
    cells[11].material_id = 0;

    cells[12].vertices[0] =  12;
    cells[12].vertices[1] =  0;
    cells[12].vertices[2] =  3;
    cells[12].vertices[3] =  15;
    cells[12].material_id = 0;

    cells[13].vertices[0] =  13;
    cells[13].vertices[1] =  1;
    cells[13].vertices[2] =  0;
    cells[13].vertices[3] =  12;
    cells[13].material_id = 0;

    cells[14].vertices[0] =  14;
    cells[14].vertices[1] =  2;
    cells[14].vertices[2] =  1;
    cells[14].vertices[3] =  13;
    cells[14].material_id = 0;

    cells[15].vertices[0] =  15;
    cells[15].vertices[1] =  3;
    cells[15].vertices[2] =  2;
    cells[15].vertices[3] =  14;
    cells[15].material_id = 0;

    // Must call this to be able to create a
    // correct triangulation in dealii, read
    // GridReordering<> doc
    GridReordering<dim,spacedim>::reorder_cells (cells);
    tria.create_triangulation_compatibility (vertices, cells, SubCellData());
  }



  template<>
  void
  parallelogram (Triangulation<3>  &tria,
                 const Point<3> (&corners)[3],
                 const bool         colorize)
  {
    Assert (false, ExcNotImplemented());
  }

  template<>
  void
  parallelogram (Triangulation<1>  &tria,
                 const Point<1> (&corners)[1],
                 const bool         colorize)
  {
    Assert (false, ExcNotImplemented());
  }

// Implementation for 2D only
  template<>
  void
  parallelogram (Triangulation<2>  &tria,
                 const Point<2> (&corners)[2],
                 const bool         colorize)
  {
    std::vector<Point<2> > vertices (GeometryInfo<2>::vertices_per_cell);

    vertices[1] = corners[0];
    vertices[2] = corners[1];
    vertices[3] = vertices[1] + vertices[2];
    // Prepare cell data
    std::vector<CellData<2> > cells (1);
    for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
      cells[0].vertices[i] = i;
    cells[0].material_id = 0;

    tria.create_triangulation (vertices, cells, SubCellData());

    // Assign boundary indicators
    if (colorize)
      colorize_hyper_rectangle (tria);
  }


  template<>
  void
  parallelogram (Triangulation<1>  &tria,
                 const Tensor<2,1> &corners,
                 const bool         colorize)
  {
    Assert (false, ExcNotImplemented());
  }

  template<>
  void
  parallelogram (Triangulation<3>  &tria,
                 const Tensor<2,3> &corners,
                 const bool         colorize)
  {
    Assert (false, ExcNotImplemented());
  }


  template<>
  void
  parallelogram (Triangulation<2>  &tria,
                 const Tensor<2,2> &corners,
                 const bool         colorize)
  {
    // simply pass everything to the other function of same name
    const Point<2> x[2] = { corners[0], corners[1] };
    parallelogram (tria, x, colorize);
  }



// Parallelepiped implementation in 1d, 2d, and 3d. @note The
// implementation in 1d is similar to hyper_rectangle(), and in 2d is
// similar to parallelogram().
//
// The GridReordering::reorder_grid is made use of towards the end of
// this function. Thus the triangulation is explicitly constructed for
// all dim here since it is slightly different in that respect
// (cf. hyper_rectangle(), parallelogram()).
  template<int dim>
  void
  parallelepiped (Triangulation<dim>  &tria,
                  const Point<dim>   (&corners) [dim],
                  const bool           colorize)
  {
    // Check that none of the user defined vertices overlap
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=i+1; j<dim; ++j)
        Assert ((corners[i]!=corners[j]),
                ExcMessage ("Invalid distance between corner points of parallelepiped."));

    // Note: vertex[0] is the origin and is initialised as so here:
    std::vector<Point<dim> > vertices (GeometryInfo<dim>::vertices_per_cell);

    switch (dim)
      {
      // A line (1d parallelepiped)
      case 1:
        vertices[1] = corners[0];
        break;

      // A parallelogram (2d parallelepiped)
      case 2:
        // assign corners to vertices:
        vertices[1] = corners[0];
        vertices[2] = corners[1];

        // compose the remaining vertex:
        vertices[3] = vertices[1] + vertices[2];
        break;

      // A parallelepiped (3d parallelepiped)
      case 3:
        // assign corners to vertices:
        vertices[1] = corners[0];
        vertices[2] = corners[1];
        vertices[4] = corners[2];

        // compose the remaining vertices:
        vertices[3] = vertices[1] + vertices[2];
        vertices[5] = vertices[1] + vertices[4];
        vertices[6] = vertices[2] + vertices[4];
        vertices[7] = vertices[1] + vertices[2] + vertices[4];
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    // Prepare cell data and wipe material identity
    std::vector<CellData<dim> > cells (1);
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      cells[0].vertices[i] = i;
    cells[0].material_id = 0;

    // Check ordering of vertices and create triangulation
    GridReordering<dim>::reorder_cells (cells);
    tria.create_triangulation (vertices, cells, SubCellData());

    // Finally assign boundary indicators according to hyper_rectangle
    if (colorize)
      colorize_hyper_rectangle (tria);
  }

  template<int dim>
  void
  subdivided_parallelepiped (Triangulation<dim>  &tria,
                             const unsigned int      n_subdivisions,
                             const Point<dim>   (&corners) [dim],
                             const bool           colorize)
  {
    // Equalise number of subdivisions in each dim-direction, heir
    // validity will be checked later
    unsigned int (n_subdivisions_) [dim];
    for (unsigned int i=0; i<dim; ++i)
      n_subdivisions_[i] = n_subdivisions;

    // and call the function below
    subdivided_parallelepiped (tria, n_subdivisions_,
                               corners,
                               colorize);
  }

  template<int dim>
  void
  subdivided_parallelepiped (Triangulation<dim>  &tria,
                             const unsigned int    ( n_subdivisions) [dim],
                             const Point<dim>   (&corners) [dim],
                             const bool           colorize)
  {
    // Zero n_subdivisions is the origin only, which makes no sense
    for (unsigned int i=0; i<dim; ++i)
      Assert (n_subdivisions[i]>0, ExcInvalidRepetitions(n_subdivisions[i]));

    // Check corners do not overlap (unique)
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=i+1; j<dim; ++j)
        Assert ((corners[i]!=corners[j]),
                ExcMessage ("Invalid distance between corner points of parallelepiped."));

    // Create a list of points
    Point<dim> (delta) [dim];

    for (unsigned int i=0; i<dim; ++i)
      delta[i] = corners[i]/n_subdivisions[i];
    std::vector<Point<dim> > points;

    switch (dim)
      {
      case 1:
        for (unsigned int x=0; x<=n_subdivisions[0]; ++x)
          points.push_back (Point<dim> (x*delta[0]));
        break;

      case 2:
        for (unsigned int y=0; y<=n_subdivisions[1]; ++y)
          for (unsigned int x=0; x<=n_subdivisions[0]; ++x)
            points.push_back (Point<dim> (x*delta[0] + y*delta[1]));
        break;

      case 3:
        for (unsigned int z=0; z<=n_subdivisions[2]; ++z)
          for (unsigned int y=0; y<=n_subdivisions[1]; ++y)
            for (unsigned int x=0; x<=n_subdivisions[0]; ++x)
              points.push_back (Point<dim> (x*delta[0] + y*delta[1] + z*delta[2]));
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    // Prepare cell data
    unsigned int n_cells = 1;
    for (unsigned int i=0; i<dim; ++i)
      n_cells *= n_subdivisions[i];
    std::vector<CellData<dim> > cells (n_cells);

    // Create fixed ordering of
    switch (dim)
      {
      case 1:
        for (unsigned int x=0; x<n_subdivisions[0]; ++x)
          {
            cells[x].vertices[0] = x;
            cells[x].vertices[1] = x+1;

            // wipe material id
            cells[x].material_id = 0;
          }
        break;

      case 2:
      {
        // Shorthand
        const unsigned int n_dy = n_subdivisions[1];
        const unsigned int n_dx = n_subdivisions[0];

        for (unsigned int y=0; y<n_dy; ++y)
          for (unsigned int x=0; x<n_dx; ++x)
            {
              const unsigned int    c = y*n_dx         + x;
              cells[c].vertices[0] = y*(n_dx+1)     + x;
              cells[c].vertices[1] = y*(n_dx+1)     + x+1;
              cells[c].vertices[2] = (y+1)*(n_dx+1) + x;
              cells[c].vertices[3] = (y+1)*(n_dx+1) + x+1;

              // wipe material id
              cells[c].material_id = 0;
            }
      }
      break;

      case 3:
      {
        // Shorthand
        const unsigned int n_dz = n_subdivisions[2];
        const unsigned int n_dy = n_subdivisions[1];
        const unsigned int n_dx = n_subdivisions[0];

        for (unsigned int z=0; z<n_dz; ++z)
          for (unsigned int y=0; y<n_dy; ++y)
            for (unsigned int x=0; x<n_dx; ++x)
              {
                const unsigned int    c = z*n_dy*n_dx             + y*n_dx         + x;

                cells[c].vertices[0] = z*(n_dy+1)*(n_dx+1)     + y*(n_dx+1)     + x;
                cells[c].vertices[1] = z*(n_dy+1)*(n_dx+1)     + y*(n_dx+1)     + x+1;
                cells[c].vertices[2] = z*(n_dy+1)*(n_dx+1)     + (y+1)*(n_dx+1) + x;
                cells[c].vertices[3] = z*(n_dy+1)*(n_dx+1)     + (y+1)*(n_dx+1) + x+1;
                cells[c].vertices[4] = (z+1)*(n_dy+1)*(n_dx+1) + y*(n_dx+1)     + x;
                cells[c].vertices[5] = (z+1)*(n_dy+1)*(n_dx+1) + y*(n_dx+1)     + x+1;
                cells[c].vertices[6] = (z+1)*(n_dy+1)*(n_dx+1) + (y+1)*(n_dx+1) + x;
                cells[c].vertices[7] = (z+1)*(n_dy+1)*(n_dx+1) + (y+1)*(n_dx+1) + x+1;

                // wipe material id
                cells[c].material_id = 0;
              }
        break;
      }

      default:
        Assert (false, ExcNotImplemented());
      }

    // Create triangulation
    tria.create_triangulation (points, cells, SubCellData());

    // Finally assign boundary indicators according to hyper_rectangle
    if (colorize)
      colorize_hyper_rectangle (tria);
  }


  template <int dim>
  void
  subdivided_hyper_cube (Triangulation<dim> &tria,
                         const unsigned int  repetitions,
                         const double        left,
                         const double        right)
  {
    Assert (repetitions >= 1, ExcInvalidRepetitions(repetitions));
    Assert (left < right,
            ExcMessage ("Invalid left-to-right bounds of hypercube"));

    // first generate the necessary
    // points
    const double delta = (right-left)/repetitions;
    std::vector<Point<dim> > points;
    switch (dim)
      {
      case 1:
        for (unsigned int x=0; x<=repetitions; ++x)
          points.push_back (Point<dim> (left+x*delta));
        break;

      case 2:
        for (unsigned int y=0; y<=repetitions; ++y)
          for (unsigned int x=0; x<=repetitions; ++x)
            points.push_back (Point<dim> (left+x*delta,
                                          left+y*delta));
        break;

      case 3:
        for (unsigned int z=0; z<=repetitions; ++z)
          for (unsigned int y=0; y<=repetitions; ++y)
            for (unsigned int x=0; x<=repetitions; ++x)
              points.push_back (Point<dim> (left+x*delta,
                                            left+y*delta,
                                            left+z*delta));
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    // next create the cells
    // Prepare cell data
    std::vector<CellData<dim> > cells;
    // Define these as abbreviations
    // for the step sizes below. The
    // number of points in a single
    // direction is repetitions+1
    const unsigned int dy = repetitions+1;
    const unsigned int dz = dy*dy;
    switch (dim)
      {
      case 1:
        cells.resize (repetitions);
        for (unsigned int x=0; x<repetitions; ++x)
          {
            cells[x].vertices[0] = x;
            cells[x].vertices[1] = x+1;
            cells[x].material_id = 0;
          }
        break;

      case 2:
        cells.resize (repetitions*repetitions);
        for (unsigned int y=0; y<repetitions; ++y)
          for (unsigned int x=0; x<repetitions; ++x)
            {
              const unsigned int c = x  +y*repetitions;
              cells[c].vertices[0] = x  +y*dy;
              cells[c].vertices[1] = x+1+y*dy;
              cells[c].vertices[2] = x  +(y+1)*dy;
              cells[c].vertices[3] = x+1+(y+1)*dy;
              cells[c].material_id = 0;
            }
        break;

      case 3:
        cells.resize (repetitions*repetitions*repetitions);
        for (unsigned int z=0; z<repetitions; ++z)
          for (unsigned int y=0; y<repetitions; ++y)
            for (unsigned int x=0; x<repetitions; ++x)
              {
                const unsigned int c = x+y*repetitions
                                       +z*repetitions*repetitions;
                cells[c].vertices[0] = x  +y*dy    +z*dz;
                cells[c].vertices[1] = x+1+y*dy    +z*dz;
                cells[c].vertices[2] = x  +(y+1)*dy+z*dz;
                cells[c].vertices[3] = x+1+(y+1)*dy+z*dz;
                cells[c].vertices[4] = x  +y*dy    +(z+1)*dz;
                cells[c].vertices[5] = x+1+y*dy    +(z+1)*dz;
                cells[c].vertices[6] = x  +(y+1)*dy+(z+1)*dz;
                cells[c].vertices[7] = x+1+(y+1)*dy+(z+1)*dz;
                cells[c].material_id = 0;
              }
        break;

      default:
        // should be trivial to
        // do for 3d as well, but
        // am too tired at this
        // point of the night to
        // do that...
        //
        // contributions are welcome!
        Assert (false, ExcNotImplemented());
      }

    tria.create_triangulation (points, cells, SubCellData());
  }



  template <int dim>
  void
  subdivided_hyper_rectangle (
    Triangulation<dim>              &tria,
    const std::vector<unsigned int> &repetitions,
    const Point<dim>                &p_1,
    const Point<dim>                &p_2,
    const bool                       colorize)
  {
    // contributed by Joerg R. Weimar
    // (j.weimar@jweimar.de) 2003
    Assert(repetitions.size() == dim,
           ExcInvalidRepetitionsDimension(dim));
    // First, normalize input such that
    // p1 is lower in all coordinate
    // directions.
    Point<dim> p1(p_1);
    Point<dim> p2(p_2);

    for (unsigned int i=0; i<dim; ++i)
      if (p1(i) > p2(i))
        std::swap (p1(i), p2(i));

    // then check that all repetitions
    // are >= 1, and calculate deltas
    // convert repetitions from double
    // to int by taking the ceiling.
    Point<dim> delta;

    for (unsigned int i=0; i<dim; ++i)
      {
        Assert (repetitions[i] >= 1, ExcInvalidRepetitions(repetitions[i]));

        delta[i] = (p2[i]-p1[i])/repetitions[i];
      }

    // then generate the necessary
    // points
    std::vector<Point<dim> > points;
    switch (dim)
      {
      case 1:
        for (unsigned int x=0; x<=repetitions[0]; ++x)
          points.push_back (Point<dim> (p1[0]+x*delta[0]));
        break;

      case 2:
        for (unsigned int y=0; y<=repetitions[1]; ++y)
          for (unsigned int x=0; x<=repetitions[0]; ++x)
            points.push_back (Point<dim> (p1[0]+x*delta[0],
                                          p1[1]+y*delta[1]));
        break;

      case 3:
        for (unsigned int z=0; z<=repetitions[2]; ++z)
          for (unsigned int y=0; y<=repetitions[1]; ++y)
            for (unsigned int x=0; x<=repetitions[0]; ++x)
              points.push_back (Point<dim> (p1[0]+x*delta[0],
                                            p1[1]+y*delta[1],
                                            p1[2]+z*delta[2]));
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    // next create the cells
    // Prepare cell data
    std::vector<CellData<dim> > cells;
    switch (dim)
      {
      case 1:
      {
        cells.resize (repetitions[0]);
        for (unsigned int x=0; x<repetitions[0]; ++x)
          {
            cells[x].vertices[0] = x;
            cells[x].vertices[1] = x+1;
            cells[x].material_id = 0;
          }
        break;
      }

      case 2:
      {
        cells.resize (repetitions[1]*repetitions[0]);
        for (unsigned int y=0; y<repetitions[1]; ++y)
          for (unsigned int x=0; x<repetitions[0]; ++x)
            {
              const unsigned int c = x+y*repetitions[0];
              cells[c].vertices[0] = y*(repetitions[0]+1)+x;
              cells[c].vertices[1] = y*(repetitions[0]+1)+x+1;
              cells[c].vertices[2] = (y+1)*(repetitions[0]+1)+x;
              cells[c].vertices[3] = (y+1)*(repetitions[0]+1)+x+1;
              cells[c].material_id = 0;
            }
        break;
      }

      case 3:
      {
        const unsigned int n_x  = (repetitions[0]+1);
        const unsigned int n_xy = (repetitions[0]+1)*(repetitions[1]+1);

        cells.resize (repetitions[2]*repetitions[1]*repetitions[0]);
        for (unsigned int z=0; z<repetitions[2]; ++z)
          for (unsigned int y=0; y<repetitions[1]; ++y)
            for (unsigned int x=0; x<repetitions[0]; ++x)
              {
                const unsigned int c = x+y*repetitions[0] +
                                       z*repetitions[0]*repetitions[1];
                cells[c].vertices[0] = z*n_xy + y*n_x + x;
                cells[c].vertices[1] = z*n_xy + y*n_x + x+1;
                cells[c].vertices[2] = z*n_xy + (y+1)*n_x + x;
                cells[c].vertices[3] = z*n_xy + (y+1)*n_x + x+1;
                cells[c].vertices[4] = (z+1)*n_xy + y*n_x + x;
                cells[c].vertices[5] = (z+1)*n_xy + y*n_x + x+1;
                cells[c].vertices[6] = (z+1)*n_xy + (y+1)*n_x + x;
                cells[c].vertices[7] = (z+1)*n_xy + (y+1)*n_x + x+1;
                cells[c].material_id = 0;
              }
        break;

      }

      default:
        Assert (false, ExcNotImplemented());
      }

    tria.create_triangulation (points, cells, SubCellData());

    if (colorize)
      {
        // to colorize, run through all
        // faces of all cells and set
        // boundary indicator to the
        // correct value if it was 0.

        // use a large epsilon to
        // compare numbers to avoid
        // roundoff problems.
        const double epsilon
          = 0.01 * *std::min_element (&delta[0], &delta[0]+dim);
        Assert (epsilon > 0,
                ExcMessage ("The distance between corner points must be positive."))

        // actual code is external since
        // 1-D is different from 2/3D.
        colorize_subdivided_hyper_rectangle (tria, p1, p2, epsilon);
      }
  }



  template <int dim>
  void
  subdivided_hyper_rectangle(
    Triangulation<dim>              &tria,
    const std::vector<std::vector<double> > &step_sz,
    const Point<dim>                &p_1,
    const Point<dim>                &p_2,
    const bool                       colorize)
  {
    // contributed by Joerg R. Weimar
    // (j.weimar@jweimar.de) 2003
    // modified by Yaqi Wang 2006
    Assert(step_sz.size() == dim,
           ExcInvalidRepetitionsDimension(dim));


    // First, normalize input such that
    // p1 is lower in all coordinate
    // directions.

    // and check the consistency of
    // step sizes, i.e. that they all
    // add up to the sizes specified by
    // p_1 and p_2
    Point<dim> p1(p_1);
    Point<dim> p2(p_2);
    std::vector< std::vector<double> > step_sizes(step_sz);

    for (unsigned int i=0; i<dim; ++i)
      {
        if (p1(i) > p2(i))
          {
            std::swap (p1(i), p2(i));
            std::reverse (step_sizes[i].begin(), step_sizes[i].end());
          }

        double x = 0;
        for (unsigned int j=0; j<step_sizes.at(i).size(); j++)
          x += step_sizes[i][j];
        Assert(std::fabs(x - (p2(i)-p1(i))) <= 1e-12*std::fabs(x),
               ExcInvalidRepetitions (i) );
      }


    // then generate the necessary
    // points
    std::vector<Point<dim> > points;
    switch (dim)
      {
      case 1:
      {
        double x=0;
        for (unsigned int i=0; ; ++i)
          {
            points.push_back (Point<dim> (p1[0]+x));

            // form partial sums. in
            // the last run through
            // avoid accessing
            // non-existent values
            // and exit early instead
            if (i == step_sizes[0].size())
              break;

            x += step_sizes[0][i];
          }
        break;
      }

      case 2:
      {
        double y=0;
        for (unsigned int j=0; ; ++j)
          {
            double x=0;
            for (unsigned int i=0; ; ++i)
              {
                points.push_back (Point<dim> (p1[0]+x,
                                              p1[1]+y));
                if (i == step_sizes[0].size())
                  break;

                x += step_sizes[0][i];
              }

            if (j == step_sizes[1].size())
              break;

            y += step_sizes[1][j];
          }
        break;

      }
      case 3:
      {
        double z=0;
        for (unsigned int k=0; ; ++k)
          {
            double y=0;
            for (unsigned int j=0; ; ++j)
              {
                double x=0;
                for (unsigned int i=0; ; ++i)
                  {
                    points.push_back (Point<dim> (p1[0]+x,
                                                  p1[1]+y,
                                                  p1[2]+z));
                    if (i == step_sizes[0].size())
                      break;

                    x += step_sizes[0][i];
                  }

                if (j == step_sizes[1].size())
                  break;

                y += step_sizes[1][j];
              }

            if (k == step_sizes[2].size())
              break;

            z += step_sizes[2][k];
          }
        break;
      }

      default:
        Assert (false, ExcNotImplemented());
      }

    // next create the cells
    // Prepare cell data
    std::vector<CellData<dim> > cells;
    switch (dim)
      {
      case 1:
      {
        cells.resize (step_sizes[0].size());
        for (unsigned int x=0; x<step_sizes[0].size(); ++x)
          {
            cells[x].vertices[0] = x;
            cells[x].vertices[1] = x+1;
            cells[x].material_id = 0;
          }
        break;
      }

      case 2:
      {
        cells.resize (step_sizes[1].size()*step_sizes[0].size());
        for (unsigned int y=0; y<step_sizes[1].size(); ++y)
          for (unsigned int x=0; x<step_sizes[0].size(); ++x)
            {
              const unsigned int c = x+y*step_sizes[0].size();
              cells[c].vertices[0] = y*(step_sizes[0].size()+1)+x;
              cells[c].vertices[1] = y*(step_sizes[0].size()+1)+x+1;
              cells[c].vertices[2] = (y+1)*(step_sizes[0].size()+1)+x;
              cells[c].vertices[3] = (y+1)*(step_sizes[0].size()+1)+x+1;
              cells[c].material_id = 0;
            }
        break;
      }

      case 3:
      {
        const unsigned int n_x  = (step_sizes[0].size()+1);
        const unsigned int n_xy = (step_sizes[0].size()+1)*(step_sizes[1].size()+1);

        cells.resize (step_sizes[2].size()*step_sizes[1].size()*step_sizes[0].size());
        for (unsigned int z=0; z<step_sizes[2].size(); ++z)
          for (unsigned int y=0; y<step_sizes[1].size(); ++y)
            for (unsigned int x=0; x<step_sizes[0].size(); ++x)
              {
                const unsigned int    c = x+y*step_sizes[0].size() +
                                          z*step_sizes[0].size()*step_sizes[1].size();
                cells[c].vertices[0] = z*n_xy + y*n_x + x;
                cells[c].vertices[1] = z*n_xy + y*n_x + x+1;
                cells[c].vertices[2] = z*n_xy + (y+1)*n_x + x;
                cells[c].vertices[3] = z*n_xy + (y+1)*n_x + x+1;
                cells[c].vertices[4] = (z+1)*n_xy + y*n_x + x;
                cells[c].vertices[5] = (z+1)*n_xy + y*n_x + x+1;
                cells[c].vertices[6] = (z+1)*n_xy + (y+1)*n_x + x;
                cells[c].vertices[7] = (z+1)*n_xy + (y+1)*n_x + x+1;
                cells[c].material_id = 0;
              }
        break;

      }

      default:
        Assert (false, ExcNotImplemented());
      }

    tria.create_triangulation (points, cells, SubCellData());

    if (colorize)
      {
        // to colorize, run through all
        // faces of all cells and set
        // boundary indicator to the
        // correct value if it was 0.

        // use a large epsilon to
        // compare numbers to avoid
        // roundoff problems.
        double min_size = *std::min_element (step_sizes[0].begin(),
                                             step_sizes[0].end());
        for (unsigned int i=1; i<dim; ++i)
          min_size = std::min (min_size,
                               *std::min_element (step_sizes[i].begin(),
                                                  step_sizes[i].end()));
        const double epsilon = 0.01 * min_size;

        // actual code is external since
        // 1-D is different from 2/3D.
        colorize_subdivided_hyper_rectangle (tria, p1, p2, epsilon);
      }
  }



  template <>
  void
  subdivided_hyper_rectangle (
    Triangulation<1>                             &tria,
    const std::vector< std::vector<double> >     &spacing,
    const Point<1>                               &p,
    const Table<1,types::material_id>                 &material_id,
    const bool                                    colorize)
  {
    // contributed by Yaqi Wang 2006
    Assert(spacing.size() == 1,
           ExcInvalidRepetitionsDimension(1));

    const unsigned int n_cells = material_id.size(0);

    Assert(spacing[0].size() == n_cells,
           ExcInvalidRepetitionsDimension(1));

    double delta = std::numeric_limits<double>::max();
    for (unsigned int i=0; i<n_cells; i++)
      {
        Assert (spacing[0][i] >= 0, ExcInvalidRepetitions(-1));
        delta = std::min (delta, spacing[0][i]);
      }

    // generate the necessary points
    std::vector<Point<1> > points;
    double ax = p[0];
    for (unsigned int x=0; x<=n_cells; ++x)
      {
        points.push_back (Point<1> (ax));
        if (x<n_cells)
          ax += spacing[0][x];
      }
    // create the cells
    unsigned int n_val_cells = 0;
    for (unsigned int i=0; i<n_cells; i++)
      if (material_id[i]!=numbers::invalid_material_id) n_val_cells++;

    std::vector<CellData<1> > cells(n_val_cells);
    unsigned int id = 0;
    for (unsigned int x=0; x<n_cells; ++x)
      if (material_id[x] != numbers::invalid_material_id)
        {
          cells[id].vertices[0] = x;
          cells[id].vertices[1] = x+1;
          cells[id].material_id = material_id[x];
          id++;
        }
    // create triangulation
    SubCellData t;
    GridTools::delete_unused_vertices (points, cells, t);

    tria.create_triangulation (points, cells, t);

    // set boundary indicator
    if (colorize)
      Assert (false, ExcNotImplemented());
  }


  template <>
  void
  subdivided_hyper_rectangle (
    Triangulation<2>                         &tria,
    const std::vector< std::vector<double> >     &spacing,
    const Point<2>                               &p,
    const Table<2,types::material_id>          &material_id,
    const bool                                    colorize)
  {
    // contributed by Yaqi Wang 2006
    Assert(spacing.size() == 2,
           ExcInvalidRepetitionsDimension(2));

    std::vector<unsigned int> repetitions(2);
    unsigned int n_cells = 1;
    double delta = std::numeric_limits<double>::max();
    for (unsigned int i=0; i<2; i++)
      {
        repetitions[i] = spacing[i].size();
        n_cells *= repetitions[i];
        for (unsigned int j=0; j<repetitions[i]; j++)
          {
            Assert (spacing[i][j] >= 0, ExcInvalidRepetitions(-1));
            delta = std::min (delta, spacing[i][j]);
          }
        Assert(material_id.size(i) == repetitions[i],
               ExcInvalidRepetitionsDimension(i));
      }

    // generate the necessary points
    std::vector<Point<2> > points;
    double ay = p[1];
    for (unsigned int y=0; y<=repetitions[1]; ++y)
      {
        double ax = p[0];
        for (unsigned int x=0; x<=repetitions[0]; ++x)
          {
            points.push_back (Point<2> (ax,ay));
            if (x<repetitions[0])
              ax += spacing[0][x];
          }
        if (y<repetitions[1])
          ay += spacing[1][y];
      }

    // create the cells
    unsigned int n_val_cells = 0;
    for (unsigned int i=0; i<material_id.size(0); i++)
      for (unsigned int j=0; j<material_id.size(1); j++)
        if (material_id[i][j] != numbers::invalid_material_id)
          n_val_cells++;

    std::vector<CellData<2> > cells(n_val_cells);
    unsigned int id = 0;
    for (unsigned int y=0; y<repetitions[1]; ++y)
      for (unsigned int x=0; x<repetitions[0]; ++x)
        if (material_id[x][y]!=numbers::invalid_material_id)
          {
            cells[id].vertices[0] = y*(repetitions[0]+1)+x;
            cells[id].vertices[1] = y*(repetitions[0]+1)+x+1;
            cells[id].vertices[2] = (y+1)*(repetitions[0]+1)+x;
            cells[id].vertices[3] = (y+1)*(repetitions[0]+1)+x+1;
            cells[id].material_id = material_id[x][y];
            id++;
          }

    // create triangulation
    SubCellData t;
    GridTools::delete_unused_vertices (points, cells, t);

    tria.create_triangulation (points, cells, t);

    // set boundary indicator
    if (colorize)
      {
        double eps = 0.01 * delta;
        Triangulation<2>::cell_iterator cell = tria.begin(),
                                        endc = tria.end();
        for (; cell !=endc; ++cell)
          {
            Point<2> cell_center = cell->center();
            for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
              if (cell->face(f)->boundary_indicator() == 0)
                {
                  Point<2> face_center = cell->face(f)->center();
                  for (unsigned int i=0; i<2; ++i)
                    {
                      if (face_center[i]<cell_center[i]-eps)
                        cell->face(f)->set_boundary_indicator(i*2);
                      if (face_center[i]>cell_center[i]+eps)
                        cell->face(f)->set_boundary_indicator(i*2+1);
                    }
                }
          }
      }
  }


  template <>
  void
  subdivided_hyper_rectangle (
    Triangulation<3>                           &tria,
    const std::vector< std::vector<double> >     &spacing,
    const Point<3>                             &p,
    const Table<3,types::material_id>               &material_id,
    const bool                                    colorize)
  {
    // contributed by Yaqi Wang 2006
    const unsigned int dim = 3;

    Assert(spacing.size() == dim,
           ExcInvalidRepetitionsDimension(dim));

    std::vector<unsigned int > repetitions(dim);
    unsigned int n_cells = 1;
    double delta = std::numeric_limits<double>::max();
    for (unsigned int i=0; i<dim; i++)
      {
        repetitions[i] = spacing[i].size();
        n_cells *= repetitions[i];
        for (unsigned int j=0; j<repetitions[i]; j++)
          {
            Assert (spacing[i][j] >= 0, ExcInvalidRepetitions(-1));
            delta = std::min (delta, spacing[i][j]);
          }
        Assert(material_id.size(i) == repetitions[i],
               ExcInvalidRepetitionsDimension(i));
      }

    // generate the necessary points
    std::vector<Point<dim> > points;
    double az = p[2];
    for (unsigned int z=0; z<=repetitions[2]; ++z)
      {
        double ay = p[1];
        for (unsigned int y=0; y<=repetitions[1]; ++y)
          {
            double ax = p[0];
            for (unsigned int x=0; x<=repetitions[0]; ++x)
              {
                points.push_back (Point<dim> (ax,ay,az));
                if (x<repetitions[0])
                  ax += spacing[0][x];
              }
            if (y<repetitions[1])
              ay += spacing[1][y];
          }
        if (z<repetitions[2])
          az += spacing[2][z];
      }

    // create the cells
    unsigned int n_val_cells = 0;
    for (unsigned int i=0; i<material_id.size(0); i++)
      for (unsigned int j=0; j<material_id.size(1); j++)
        for (unsigned int k=0; k<material_id.size(2); k++)
          if (material_id[i][j][k]!=numbers::invalid_material_id)
            n_val_cells++;

    std::vector<CellData<dim> > cells(n_val_cells);
    unsigned int id = 0;
    const unsigned int n_x  = (repetitions[0]+1);
    const unsigned int n_xy = (repetitions[0]+1)*(repetitions[1]+1);
    for (unsigned int z=0; z<repetitions[2]; ++z)
      for (unsigned int y=0; y<repetitions[1]; ++y)
        for (unsigned int x=0; x<repetitions[0]; ++x)
          if (material_id[x][y][z]!=numbers::invalid_material_id)
            {
              cells[id].vertices[0] = z*n_xy + y*n_x + x;
              cells[id].vertices[1] = z*n_xy + y*n_x + x+1;
              cells[id].vertices[2] = z*n_xy + (y+1)*n_x + x;
              cells[id].vertices[3] = z*n_xy + (y+1)*n_x + x+1;
              cells[id].vertices[4] = (z+1)*n_xy + y*n_x + x;
              cells[id].vertices[5] = (z+1)*n_xy + y*n_x + x+1;
              cells[id].vertices[6] = (z+1)*n_xy + (y+1)*n_x + x;
              cells[id].vertices[7] = (z+1)*n_xy + (y+1)*n_x + x+1;
              cells[id].material_id = material_id[x][y][z];
              id++;
            }

    // create triangulation
    SubCellData t;
    GridTools::delete_unused_vertices (points, cells, t);

    tria.create_triangulation (points, cells, t);

    // set boundary indicator
    if (colorize && dim>1)
      {
        double eps = 0.01 * delta;
        Triangulation<dim>::cell_iterator cell = tria.begin(),
                                          endc = tria.end();
        for (; cell !=endc; ++cell)
          {
            Point<dim> cell_center = cell->center();
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->face(f)->boundary_indicator() == 0)
                {
                  Point<dim> face_center = cell->face(f)->center();
                  for (unsigned int i=0; i<dim; ++i)
                    {
                      if (face_center[i]<cell_center[i]-eps)
                        cell->face(f)->set_boundary_indicator(i*2);
                      if (face_center[i]>cell_center[i]+eps)
                        cell->face(f)->set_boundary_indicator(i*2+1);
                    }
                }
          }
      }
  }



  template <>
  void hyper_cube_slit (Triangulation<1> &,
                        const double,
                        const double,
                        const bool)
  {
    Assert (false, ExcNotImplemented());
  }



  template <>
  void enclosed_hyper_cube (Triangulation<1> &,
                            const double,
                            const double,
                            const double,
                            const bool)
  {
    Assert (false, ExcNotImplemented());
  }



  template <>
  void hyper_L (Triangulation<1> &,
                const double,
                const double)
  {
    Assert (false, ExcNotImplemented());
  }



  template <>
  void hyper_ball (Triangulation<1> &,
                   const Point<1> &,
                   const double)
  {
    Assert (false, ExcNotImplemented());
  }



  template <>
  void cylinder (Triangulation<1> &,
                 const double,
                 const double)
  {
    Assert (false, ExcNotImplemented());
  }



  template <>
  void truncated_cone (Triangulation<1> &,
                       const double,
                       const double,
                       const double)
  {
    Assert (false, ExcNotImplemented());
  }



  template <>
  void hyper_shell (Triangulation<1> &,
                    const Point<1> &,
                    const double,
                    const double,
                    const unsigned int ,
                    const bool)
  {
    Assert (false, ExcNotImplemented());
  }


  template <>
  void cylinder_shell (Triangulation<1> &,
                       const double,
                       const double,
                       const double,
                       const unsigned int ,
                       const unsigned int )
  {
    Assert (false, ExcNotImplemented());
  }


  template <>
  void
  half_hyper_ball (Triangulation<1> &,
                   const Point<1> &,
                   const double)
  {
    Assert (false, ExcNotImplemented());
  }


  template <>
  void
  half_hyper_shell (Triangulation<1> &,
                    const Point<1> &,
                    const double,
                    const double,
                    const unsigned int ,
                    const bool)
  {
    Assert (false, ExcNotImplemented());
  }

  template <>
  void quarter_hyper_shell (Triangulation<1> &,
                            const Point<1> &,
                            const double,
                            const double,
                            const unsigned int ,
                            const bool)
  {
    Assert (false, ExcNotImplemented());
  }

  template <>
  void enclosed_hyper_cube (Triangulation<2> &tria,
                            const double        left,
                            const double        right,
                            const double        thickness,
                            const bool          colorize)
  {
    Assert(left<right,
           ExcMessage ("Invalid left-to-right bounds of enclosed hypercube"));

    std::vector<Point<2> > vertices(16);
    double coords[4];
    coords[0] = left-thickness;
    coords[1] = left;
    coords[2] = right;
    coords[3] = right+thickness;

    unsigned int k=0;
    for (unsigned int i0=0; i0<4; ++i0)
      for (unsigned int i1=0; i1<4; ++i1)
        vertices[k++] = Point<2>(coords[i1], coords[i0]);

    const types::material_id materials[9] = { 5, 4, 6,
                                              1, 0, 2,
                                              9, 8,10
                                            };

    std::vector<CellData<2> > cells(9);
    k = 0;
    for (unsigned int i0=0; i0<3; ++i0)
      for (unsigned int i1=0; i1<3; ++i1)
        {
          cells[k].vertices[0] = i1+4*i0;
          cells[k].vertices[1] = i1+4*i0+1;
          cells[k].vertices[2] = i1+4*i0+4;
          cells[k].vertices[3] = i1+4*i0+5;
          if (colorize)
            cells[k].material_id = materials[k];
          ++k;
        }
    tria.create_triangulation (vertices,
                               cells,
                               SubCellData());       // no boundary information
  }



// Implementation for 2D only
  template <>
  void
  hyper_cube_slit (Triangulation<2> &tria,
                   const double left,
                   const double right,
                   const bool colorize)
  {
    const double rl2=(right+left)/2;
    const Point<2> vertices[10] = { Point<2>(left, left ),
                                    Point<2>(rl2,  left ),
                                    Point<2>(rl2,  rl2  ),
                                    Point<2>(left, rl2  ),
                                    Point<2>(right,left ),
                                    Point<2>(right,rl2  ),
                                    Point<2>(rl2,  right),
                                    Point<2>(left, right),
                                    Point<2>(right,right),
                                    Point<2>(rl2,  left )
                                  };
    const int cell_vertices[4][4] = { { 0,1,3,2 },
      { 9,4,2,5 },
      { 3,2,7,6 },
      { 2,5,6,8 }
    };
    std::vector<CellData<2> > cells (4, CellData<2>());
    for (unsigned int i=0; i<4; ++i)
      {
        for (unsigned int j=0; j<4; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };
    tria.create_triangulation (
      std::vector<Point<2> >(&vertices[0], &vertices[10]),
      cells,
      SubCellData());       // no boundary information

    if (colorize)
      {
        Triangulation<2>::cell_iterator cell = tria.begin();
        cell->face(1)->set_boundary_indicator(1);
        ++cell;
        cell->face(3)->set_boundary_indicator(2);
      }
  }



  template <>
  void truncated_cone (Triangulation<2> &triangulation,
                       const double radius_0,
                       const double radius_1,
                       const double half_length)
  {
    Point<2> vertices_tmp[4];

    vertices_tmp[0] = Point<2> (-half_length, -radius_0);
    vertices_tmp[1] = Point<2> (half_length, -radius_1);
    vertices_tmp[2] = Point<2> (-half_length, radius_0);
    vertices_tmp[3] = Point<2> (half_length, radius_1);

    const std::vector<Point<2> > vertices (&vertices_tmp[0], &vertices_tmp[4]);
    unsigned int cell_vertices[1][GeometryInfo<2>::vertices_per_cell];

    for (unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i)
      cell_vertices[0][i] = i;

    std::vector<CellData<2> > cells (1, CellData<2> ());

    for (unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i)
      cells[0].vertices[i] = cell_vertices[0][i];

    cells[0].material_id = 0;
    triangulation.create_triangulation (vertices, cells, SubCellData ());

    Triangulation<2>::cell_iterator cell = triangulation.begin ();

    cell->face (0)->set_boundary_indicator (1);
    cell->face (1)->set_boundary_indicator (2);

    for (unsigned int i = 2; i < 4; ++i)
      cell->face (i)->set_boundary_indicator (0);
  }



//TODO: Colorize edges as circumference, left and right radius
// Implementation for 2D only
  template <>
  void
  hyper_L (Triangulation<2> &tria,
           const double a,
           const double b)
  {
    const Point<2> vertices[8] = { Point<2> (a,a),
                                   Point<2> ((a+b)/2,a),
                                   Point<2> (b,a),
                                   Point<2> (a,(a+b)/2),
                                   Point<2> ((a+b)/2,(a+b)/2),
                                   Point<2> (b,(a+b)/2),
                                   Point<2> (a,b),
                                   Point<2> ((a+b)/2,b)
                                 };
    const int cell_vertices[3][4] = {{0, 1, 3, 4},
      {1, 2, 4, 5},
      {3, 4, 6, 7}
    };

    std::vector<CellData<2> > cells (3, CellData<2>());

    for (unsigned int i=0; i<3; ++i)
      {
        for (unsigned int j=0; j<4; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<2> >(&vertices[0], &vertices[8]),
      cells,
      SubCellData());       // no boundary information
  }



// Implementation for 2D only
  template <>
  void
  hyper_ball (Triangulation<2> &tria,
              const Point<2>   &p,
              const double      radius)
  {
    // equilibrate cell sizes at
    // transition from the inner part
    // to the radial cells
    const double a = 1./(1+std::sqrt(2.0));
    const Point<2> vertices[8] = { p+Point<2>(-1,-1) *(radius/std::sqrt(2.0)),
                                   p+Point<2>(+1,-1) *(radius/std::sqrt(2.0)),
                                   p+Point<2>(-1,-1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(+1,-1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(-1,+1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(+1,+1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(-1,+1) *(radius/std::sqrt(2.0)),
                                   p+Point<2>(+1,+1) *(radius/std::sqrt(2.0))
                                 };

    const int cell_vertices[5][4] = {{0, 1, 2, 3},
      {0, 2, 6, 4},
      {2, 3, 4, 5},
      {1, 7, 3, 5},
      {6, 4, 7, 5}
    };

    std::vector<CellData<2> > cells (5, CellData<2>());

    for (unsigned int i=0; i<5; ++i)
      {
        for (unsigned int j=0; j<4; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<2> >(&vertices[0], &vertices[8]),
      cells,
      SubCellData());       // no boundary information
  }



  template <>
  void hyper_shell (Triangulation<2> &tria,
                    const Point<2>   &center,
                    const double      inner_radius,
                    const double      outer_radius,
                    const unsigned int   n_cells,
                    const bool colorize)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());

    const double pi = numbers::PI;

    // determine the number of cells
    // for the grid. if not provided by
    // the user determine it such that
    // the length of each cell on the
    // median (in the middle between
    // the two circles) is equal to its
    // radial extent (which is the
    // difference between the two
    // radii)
    const unsigned int N = (n_cells == 0 ?
                            static_cast<unsigned int>
                            (std::ceil((2*pi* (outer_radius + inner_radius)/2) /
                                       (outer_radius - inner_radius))) :
                            n_cells);

    // set up N vertices on the
    // outer and N vertices on
    // the inner circle. the
    // first N ones are on the
    // outer one, and all are
    // numbered counter-clockwise
    std::vector<Point<2> > vertices(2*N);
    for (unsigned int i=0; i<N; ++i)
      {
        vertices[i] = Point<2>( std::cos(2*pi * i/N),
                                std::sin(2*pi * i/N)) * outer_radius;
        vertices[i+N] = vertices[i] * (inner_radius/outer_radius);

        vertices[i]   += center;
        vertices[i+N] += center;
      };

    std::vector<CellData<2> > cells (N, CellData<2>());

    for (unsigned int i=0; i<N; ++i)
      {
        cells[i].vertices[0] = i;
        cells[i].vertices[1] = (i+1)%N;
        cells[i].vertices[2] = N+i;
        cells[i].vertices[3] = N+((i+1)%N);

        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      vertices, cells, SubCellData());

    if (colorize)
      colorize_hyper_shell(tria, center, inner_radius, outer_radius);
  }


// Implementation for 2D only
  template <>
  void
  cylinder (Triangulation<2> &tria,
            const double radius,
            const double half_length)
  {
    Point<2> p1 (-half_length, -radius);
    Point<2> p2 (half_length, radius);

    hyper_rectangle(tria, p1, p2, true);

    Triangulation<2>::face_iterator f = tria.begin_face();
    Triangulation<2>::face_iterator end = tria.end_face();
    while (f != end)
      {
        switch (f->boundary_indicator())
          {
          case 0:
            f->set_boundary_indicator(1);
            break;
          case 1:
            f->set_boundary_indicator(2);
            break;
          default:
            f->set_boundary_indicator(0);
            break;
          }
        ++f;
      }
  }



// Implementation for 2D only
  template <>
  void cylinder_shell (Triangulation<2> &,
                       const double,
                       const double,
                       const double,
                       const unsigned int,
                       const unsigned int)
  {
    Assert (false, ExcNotImplemented());
  }


  template <>
  void
  half_hyper_ball (Triangulation<2> &tria,
                   const Point<2>   &p,
                   const double      radius)
  {
    // equilibrate cell sizes at
    // transition from the inner part
    // to the radial cells
    const double a = 1./(1+std::sqrt(2.0));
    const Point<2> vertices[8] = { p+Point<2>(0,-1) *radius,
                                   p+Point<2>(+1,-1) *(radius/std::sqrt(2.0)),
                                   p+Point<2>(0,-1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(+1,-1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(0,+1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(+1,+1) *(radius/std::sqrt(2.0)*a),
                                   p+Point<2>(0,+1) *radius,
                                   p+Point<2>(+1,+1) *(radius/std::sqrt(2.0))
                                 };

    const int cell_vertices[5][4] = {{0, 1, 2, 3},
      {2, 3, 4, 5},
      {1, 7, 3, 5},
      {6, 4, 7, 5}
    };

    std::vector<CellData<2> > cells (4, CellData<2>());

    for (unsigned int i=0; i<4; ++i)
      {
        for (unsigned int j=0; j<4; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<2> >(&vertices[0], &vertices[8]),
      cells,
      SubCellData());       // no boundary information

    Triangulation<2>::cell_iterator cell = tria.begin();
    Triangulation<2>::cell_iterator end = tria.end();


    while (cell != end)
      {
        for (unsigned int i=0; i<GeometryInfo<2>::faces_per_cell; ++i)
          {
            if (cell->face(i)->boundary_indicator() == numbers::internal_face_boundary_id)
              continue;

            // If x is zero, then this is part of the plane
            if (cell->face(i)->center()(0) < p(0)+1.e-5)
              cell->face(i)->set_boundary_indicator(1);
          }
        ++cell;
      }
  }



// Implementation for 2D only
  template <>
  void
  half_hyper_shell (Triangulation<2> &tria,
                    const Point<2>   &center,
                    const double      inner_radius,
                    const double      outer_radius,
                    const unsigned int   n_cells,
                    const bool colorize)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());

    const double pi     = numbers::PI;
    // determine the number of cells
    // for the grid. if not provided by
    // the user determine it such that
    // the length of each cell on the
    // median (in the middle between
    // the two circles) is equal to its
    // radial extent (which is the
    // difference between the two
    // radii)
    const unsigned int N = (n_cells == 0 ?
                            static_cast<unsigned int>
                            (std::ceil((pi* (outer_radius + inner_radius)/2) /
                                       (outer_radius - inner_radius))) :
                            n_cells);

    // set up N+1 vertices on the
    // outer and N+1 vertices on
    // the inner circle. the
    // first N+1 ones are on the
    // outer one, and all are
    // numbered counter-clockwise
    std::vector<Point<2> > vertices(2*(N+1));
    for (unsigned int i=0; i<=N; ++i)
      {
        // enforce that the x-coordinates
        // of the first and last point of
        // each half-circle are exactly
        // zero (contrary to what we may
        // compute using the imprecise
        // value of pi)
        vertices[i] =  Point<2>( ( (i==0) || (i==N) ?
                                   0 :
                                   std::cos(pi * i/N - pi/2) ),
                                 std::sin(pi * i/N - pi/2)) * outer_radius;
        vertices[i+N+1] = vertices[i] * (inner_radius/outer_radius);

        vertices[i]     += center;
        vertices[i+N+1] += center;
      };


    std::vector<CellData<2> > cells (N, CellData<2>());

    for (unsigned int i=0; i<N; ++i)
      {
        cells[i].vertices[0] = i;
        cells[i].vertices[1] = (i+1)%(N+1);
        cells[i].vertices[2] = N+1+i;
        cells[i].vertices[3] = N+1+((i+1)%(N+1));

        cells[i].material_id = 0;
      };

    tria.create_triangulation (vertices, cells, SubCellData());

    if (colorize)
      {
        Triangulation<2>::cell_iterator cell = tria.begin();
        for (; cell!=tria.end(); ++cell)
          {
            cell->face(2)->set_boundary_indicator(1);
          }
        tria.begin()->face(0)->set_boundary_indicator(3);

        tria.last()->face(1)->set_boundary_indicator(2);
      }
  }


  template <>
  void quarter_hyper_shell (Triangulation<2> &tria,
                            const Point<2>   &center,
                            const double      inner_radius,
                            const double      outer_radius,
                            const unsigned int   n_cells,
                            const bool colorize)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());

    const double pi     = numbers::PI;
    // determine the number of cells
    // for the grid. if not provided by
    // the user determine it such that
    // the length of each cell on the
    // median (in the middle between
    // the two circles) is equal to its
    // radial extent (which is the
    // difference between the two
    // radii)
    const unsigned int N = (n_cells == 0 ?
                            static_cast<unsigned int>
                            (std::ceil((pi* (outer_radius + inner_radius)/4) /
                                       (outer_radius - inner_radius))) :
                            n_cells);

    // set up N+1 vertices on the
    // outer and N+1 vertices on
    // the inner circle. the
    // first N+1 ones are on the
    // outer one, and all are
    // numbered counter-clockwise
    std::vector<Point<2> > vertices(2*(N+1));
    for (unsigned int i=0; i<=N; ++i)
      {
        // enforce that the x-coordinates
        // of the last point is exactly
        // zero (contrary to what we may
        // compute using the imprecise
        // value of pi)
        vertices[i] =  Point<2>( ( (i==N) ?
                                   0 :
                                   std::cos(pi * i/N/2) ),
                                 std::sin(pi * i/N/2)) * outer_radius;
        vertices[i+N+1] = vertices[i] * (inner_radius/outer_radius);

        vertices[i]     += center;
        vertices[i+N+1] += center;
      };


    std::vector<CellData<2> > cells (N, CellData<2>());

    for (unsigned int i=0; i<N; ++i)
      {
        cells[i].vertices[0] = i;
        cells[i].vertices[1] = (i+1)%(N+1);
        cells[i].vertices[2] = N+1+i;
        cells[i].vertices[3] = N+1+((i+1)%(N+1));

        cells[i].material_id = 0;
      };

    tria.create_triangulation (vertices, cells, SubCellData());

    if (colorize)
      {
        Triangulation<2>::cell_iterator cell = tria.begin();
        for (; cell!=tria.end(); ++cell)
          {
            cell->face(2)->set_boundary_indicator(1);
          }
        tria.begin()->face(0)->set_boundary_indicator(3);

        tria.last()->face(1)->set_boundary_indicator(2);
      }
  }



// Implementation for 3D only
  template <>
  void hyper_cube_slit (Triangulation<3> &tria,
                        const double left,
                        const double right,
                        const bool colorize)
  {
    const double rl2=(right+left)/2;
    const double len = (right-left)/2.;

    const Point<3> vertices[20] =
    {
      Point<3>(left, left , -len/2.),
      Point<3>(rl2,  left , -len/2.),
      Point<3>(rl2,  rl2  , -len/2.),
      Point<3>(left, rl2  , -len/2.),
      Point<3>(right,left , -len/2.),
      Point<3>(right,rl2  , -len/2.),
      Point<3>(rl2,  right, -len/2.),
      Point<3>(left, right, -len/2.),
      Point<3>(right,right, -len/2.),
      Point<3>(rl2,  left , -len/2.),
      Point<3>(left, left , len/2.),
      Point<3>(rl2,  left , len/2.),
      Point<3>(rl2,  rl2  , len/2.),
      Point<3>(left, rl2  , len/2.),
      Point<3>(right,left , len/2.),
      Point<3>(right,rl2  , len/2.),
      Point<3>(rl2,  right, len/2.),
      Point<3>(left, right, len/2.),
      Point<3>(right,right, len/2.),
      Point<3>(rl2,  left , len/2.)
    };
    const int cell_vertices[4][8] = { { 0,1,3,2, 10, 11, 13, 12 },
      { 9,4,2,5, 19,14, 12, 15 },
      { 3,2,7,6,13,12,17,16 },
      { 2,5,6,8,12,15,16,18 }
    };
    std::vector<CellData<3> > cells (4, CellData<3>());
    for (unsigned int i=0; i<4; ++i)
      {
        for (unsigned int j=0; j<8; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };
    tria.create_triangulation (
      std::vector<Point<3> >(&vertices[0], &vertices[20]),
      cells,
      SubCellData());       // no boundary information

    if (colorize)
      {
        Assert(false, ExcNotImplemented());
        Triangulation<3>::cell_iterator cell = tria.begin();
        cell->face(1)->set_boundary_indicator(1);
        ++cell;
        cell->face(3)->set_boundary_indicator(2);
      }
  }



// Implementation for 3D only
  template <>
  void enclosed_hyper_cube (Triangulation<3> &tria,
                            const double        left,
                            const double        right,
                            const double        thickness,
                            const bool          colorize)
  {
    Assert(left<right,
           ExcMessage ("Invalid left-to-right bounds of enclosed hypercube"));

    std::vector<Point<3> > vertices(64);
    double coords[4];
    coords[0] = left-thickness;
    coords[1] = left;
    coords[2] = right;
    coords[3] = right+thickness;

    unsigned int k=0;
    for (unsigned int z=0; z<4; ++z)
      for (unsigned int y=0; y<4; ++y)
        for (unsigned int x=0; x<4; ++x)
          vertices[k++] = Point<3>(coords[x], coords[y], coords[z]);

    const types::material_id materials[27] =
    {
      21,20,22,
      17,16,18,
      25,24,26,
      5 , 4, 6,
      1 , 0, 2,
      9 , 8,10,
      37,36,38,
      33,32,34,
      41,40,42
    };

    std::vector<CellData<3> > cells(27);
    k = 0;
    for (unsigned int z=0; z<3; ++z)
      for (unsigned int y=0; y<3; ++y)
        for (unsigned int x=0; x<3; ++x)
          {
            cells[k].vertices[0] = x+4*y+16*z;
            cells[k].vertices[1] = x+4*y+16*z+1;
            cells[k].vertices[2] = x+4*y+16*z+4;
            cells[k].vertices[3] = x+4*y+16*z+5;
            cells[k].vertices[4] = x+4*y+16*z+16;
            cells[k].vertices[5] = x+4*y+16*z+17;
            cells[k].vertices[6] = x+4*y+16*z+20;
            cells[k].vertices[7] = x+4*y+16*z+21;
            if (colorize)
              cells[k].material_id = materials[k];
            ++k;
          }
    tria.create_triangulation (
      vertices,
      cells,
      SubCellData());       // no boundary information
  }



  template <>
  void truncated_cone (Triangulation<3> &triangulation,
                       const double radius_0,
                       const double radius_1,
                       const double half_length)
  {
    // Determine number of cells and vertices
    const unsigned int
    n_cells = static_cast<unsigned int>(std::ceil (half_length /
                                                   std::max (radius_0,
                                                             radius_1)));
    const unsigned int n_vertices = 4 * (n_cells + 1);
    std::vector<Point<3> > vertices_tmp(n_vertices);

    vertices_tmp[0] = Point<3> (-half_length, 0, -radius_0);
    vertices_tmp[1] = Point<3> (-half_length, radius_0, 0);
    vertices_tmp[2] = Point<3> (-half_length, -radius_0, 0);
    vertices_tmp[3] = Point<3> (-half_length, 0, radius_0);

    const double dx = 2 * half_length / n_cells;

    for (unsigned int i = 0; i < n_cells; ++i)
      {
        vertices_tmp[4 * (i + 1)]
          = vertices_tmp[4 * i] +
            Point<3> (dx, 0, 0.5 * (radius_0 - radius_1) * dx / half_length);
        vertices_tmp[4 * i + 5]
          = vertices_tmp[4 * i + 1] +
            Point<3> (dx, 0.5 * (radius_1 - radius_0) * dx / half_length, 0);
        vertices_tmp[4 * i + 6]
          = vertices_tmp[4 * i + 2] +
            Point<3> (dx, 0.5 * (radius_0 - radius_1) * dx / half_length, 0);
        vertices_tmp[4 * i + 7]
          = vertices_tmp[4 * i + 3] +
            Point<3> (dx, 0, 0.5 * (radius_1 - radius_0) * dx / half_length);
      }

    const std::vector<Point<3> > vertices (vertices_tmp.begin(),
                                           vertices_tmp.end());
    Table<2,unsigned int> cell_vertices(n_cells,GeometryInfo<3>::vertices_per_cell);

    for (unsigned int i = 0; i < n_cells; ++i)
      for (unsigned int j = 0; j < GeometryInfo<3>::vertices_per_cell; ++j)
        cell_vertices[i][j] = 4 * i + j;

    std::vector<CellData<3> > cells (n_cells, CellData<3> ());

    for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < GeometryInfo<3>::vertices_per_cell; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];

        cells[i].material_id = 0;
      }

    triangulation.create_triangulation (vertices, cells, SubCellData ());

    for (Triangulation<3>::cell_iterator cell = triangulation.begin ();
         cell != triangulation.end (); ++cell)
      {
        if (cell->vertex (0) (0) == -half_length)
          {
            cell->face (4)->set_boundary_indicator (1);

            for (unsigned int i = 0; i < 4; ++i)
              cell->line (i)->set_boundary_indicator (0);
          }

        if (cell->vertex (4) (0) == half_length)
          {
            cell->face (5)->set_boundary_indicator (2);

            for (unsigned int i = 4; i < 8; ++i)
              cell->line (i)->set_boundary_indicator (0);
          }

        for (unsigned int i = 0; i < 4; ++i)
          cell->face (i)->set_boundary_indicator (0);
      }
  }


// Implementation for 3D only
  template <>
  void
  hyper_L (Triangulation<3> &tria,
           const double      a,
           const double      b)
  {
    // we slice out the top back right
    // part of the cube
    const Point<3> vertices[26]
    =
    {
      // front face of the big cube
      Point<3> (a,      a,a),
      Point<3> ((a+b)/2,a,a),
      Point<3> (b,      a,a),
      Point<3> (a,      a,(a+b)/2),
      Point<3> ((a+b)/2,a,(a+b)/2),
      Point<3> (b,      a,(a+b)/2),
      Point<3> (a,      a,b),
      Point<3> ((a+b)/2,a,b),
      Point<3> (b,      a,b),
      // middle face of the big cube
      Point<3> (a,      (a+b)/2,a),
      Point<3> ((a+b)/2,(a+b)/2,a),
      Point<3> (b,      (a+b)/2,a),
      Point<3> (a,      (a+b)/2,(a+b)/2),
      Point<3> ((a+b)/2,(a+b)/2,(a+b)/2),
      Point<3> (b,      (a+b)/2,(a+b)/2),
      Point<3> (a,      (a+b)/2,b),
      Point<3> ((a+b)/2,(a+b)/2,b),
      Point<3> (b,      (a+b)/2,b),
      // back face of the big cube
      // last (top right) point is missing
      Point<3> (a,      b,a),
      Point<3> ((a+b)/2,b,a),
      Point<3> (b,      b,a),
      Point<3> (a,      b,(a+b)/2),
      Point<3> ((a+b)/2,b,(a+b)/2),
      Point<3> (b,      b,(a+b)/2),
      Point<3> (a,      b,b),
      Point<3> ((a+b)/2,b,b)
    };
    const int cell_vertices[7][8] = {{0, 1, 9, 10, 3, 4, 12, 13},
      {1, 2, 10, 11, 4, 5, 13, 14},
      {3, 4, 12, 13, 6, 7, 15, 16},
      {4, 5, 13, 14, 7, 8, 16, 17},
      {9, 10, 18, 19, 12, 13, 21, 22},
      {10, 11, 19, 20, 13, 14, 22, 23},
      {12, 13, 21, 22, 15, 16, 24, 25}
    };

    std::vector<CellData<3> > cells (7, CellData<3>());

    for (unsigned int i=0; i<7; ++i)
      {
        for (unsigned int j=0; j<8; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<3> >(&vertices[0], &vertices[26]),
      cells,
      SubCellData());       // no boundary information
  }



// Implementation for 3D only
  template <>
  void
  hyper_ball (Triangulation<3> &tria,
              const Point<3>   &p,
              const double radius)
  {
    const double a = 1./(1+std::sqrt(3.0)); // equilibrate cell sizes at transition
    // from the inner part to the radial
    // cells
    const unsigned int n_vertices = 16;
    const Point<3> vertices[n_vertices]
    =
    {
      // first the vertices of the inner
      // cell
      p+Point<3>(-1,-1,-1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(+1,-1,-1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(+1,-1,+1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(-1,-1,+1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(-1,+1,-1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(+1,+1,-1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(+1,+1,+1) *(radius/std::sqrt(3.0)*a),
      p+Point<3>(-1,+1,+1) *(radius/std::sqrt(3.0)*a),
      // now the eight vertices at
      // the outer sphere
      p+Point<3>(-1,-1,-1) *(radius/std::sqrt(3.0)),
      p+Point<3>(+1,-1,-1) *(radius/std::sqrt(3.0)),
      p+Point<3>(+1,-1,+1) *(radius/std::sqrt(3.0)),
      p+Point<3>(-1,-1,+1) *(radius/std::sqrt(3.0)),
      p+Point<3>(-1,+1,-1) *(radius/std::sqrt(3.0)),
      p+Point<3>(+1,+1,-1) *(radius/std::sqrt(3.0)),
      p+Point<3>(+1,+1,+1) *(radius/std::sqrt(3.0)),
      p+Point<3>(-1,+1,+1) *(radius/std::sqrt(3.0)),
    };

    // one needs to draw the seven cubes to
    // understand what's going on here
    const unsigned int n_cells = 7;
    const int cell_vertices[n_cells][8] = {{0, 1, 4, 5, 3, 2, 7, 6}, // center
      {8, 9, 12, 13, 0, 1, 4, 5}, // bottom
      {9, 13, 1, 5, 10, 14, 2, 6}, // right
      {11, 10, 3, 2, 15, 14, 7, 6}, // top
      {8, 0, 12, 4, 11, 3, 15, 7}, // left
      {8, 9, 0, 1, 11, 10, 3, 2}, // front
      {12, 4, 13, 5, 15, 7, 14, 6}
    }; // back

    std::vector<CellData<3> > cells (n_cells, CellData<3>());

    for (unsigned int i=0; i<n_cells; ++i)
      {
        for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<3> >(&vertices[0], &vertices[n_vertices]),
      cells,
      SubCellData());       // no boundary information
  }



// Implementation for 3D only
  template <>
  void
  cylinder (Triangulation<3> &tria,
            const double radius,
            const double half_length)
  {
    // Copy the base from hyper_ball<3>
    // and transform it to yz
    const double d = radius/std::sqrt(2.0);
    const double a = d/(1+std::sqrt(2.0));
    Point<3> vertices[24] =
    {
      Point<3>(-d, -half_length,-d),
      Point<3>( d, -half_length,-d),
      Point<3>(-a, -half_length,-a),
      Point<3>( a, -half_length,-a),
      Point<3>(-a, -half_length, a),
      Point<3>( a, -half_length, a),
      Point<3>(-d, -half_length, d),
      Point<3>( d, -half_length, d),
      Point<3>(-d, 0,-d),
      Point<3>( d, 0,-d),
      Point<3>(-a, 0,-a),
      Point<3>( a, 0,-a),
      Point<3>(-a, 0, a),
      Point<3>( a, 0, a),
      Point<3>(-d, 0, d),
      Point<3>( d, 0, d),
      Point<3>(-d, half_length,-d),
      Point<3>( d, half_length,-d),
      Point<3>(-a, half_length,-a),
      Point<3>( a, half_length,-a),
      Point<3>(-a, half_length, a),
      Point<3>( a, half_length, a),
      Point<3>(-d, half_length, d),
      Point<3>( d, half_length, d),
    };
    // Turn cylinder such that y->x
    for (unsigned int i=0; i<24; ++i)
      {
        const double h = vertices[i](1);
        vertices[i](1) = -vertices[i](0);
        vertices[i](0) = h;
      }

    int cell_vertices[10][8] =
    {
      {0, 1, 8, 9, 2, 3, 10, 11},
      {0, 2, 8, 10, 6, 4, 14, 12},
      {2, 3, 10, 11, 4, 5, 12, 13},
      {1, 7, 9, 15, 3, 5, 11, 13},
      {6, 4, 14, 12, 7, 5, 15, 13}
    };
    for (unsigned int i=0; i<5; ++i)
      for (unsigned int j=0; j<8; ++j)
        cell_vertices[i+5][j] = cell_vertices[i][j]+8;

    std::vector<CellData<3> > cells (10, CellData<3>());

    for (unsigned int i=0; i<10; ++i)
      {
        for (unsigned int j=0; j<8; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<3> >(&vertices[0], &vertices[24]),
      cells,
      SubCellData());       // no boundary information

    // set boundary indicators for the
    // faces at the ends to 1 and 2,
    // respectively. note that we also
    // have to deal with those lines
    // that are purely in the interior
    // of the ends. we determine whether
    // an edge is purely in the
    // interior if one of its vertices
    // is at coordinates '+-a' as set
    // above
    Triangulation<3>::cell_iterator cell = tria.begin();
    Triangulation<3>::cell_iterator end = tria.end();

    for (; cell != end; ++cell)
      for (unsigned int i=0; i<GeometryInfo<3>::faces_per_cell; ++i)
        if (cell->at_boundary(i))
          {
            if (cell->face(i)->center()(0) > half_length-1.e-5)
              {
                cell->face(i)->set_boundary_indicator(2);

                for (unsigned int e=0; e<GeometryInfo<3>::lines_per_face; ++e)
                  if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a) ||
                      (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a) ||
                      (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a) ||
                      (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a))
                    cell->face(i)->line(e)->set_boundary_indicator(2);
              }
            else if (cell->face(i)->center()(0) < -half_length+1.e-5)
              {
                cell->face(i)->set_boundary_indicator(1);

                for (unsigned int e=0; e<GeometryInfo<3>::lines_per_face; ++e)
                  if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a) ||
                      (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a) ||
                      (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a) ||
                      (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a))
                    cell->face(i)->line(e)->set_boundary_indicator(1);
              }
          }
  }



// Implementation for 3D only
  template <>
  void
  half_hyper_ball (Triangulation<3> &tria,
                   const Point<3> &center,
                   const double radius)
  {
    // These are for the two lower squares
    const double d = radius/std::sqrt(2.0);
    const double a = d/(1+std::sqrt(2.0));
    // These are for the two upper square
    const double b = a/2.0;
    const double c = d/2.0;
    // And so are these
    const double hb = radius*std::sqrt(3.0)/4.0;
    const double hc = radius*std::sqrt(3.0)/2.0;

    Point<3> vertices[16] =
    {
      center+Point<3>( 0,  d, -d),
      center+Point<3>( 0, -d, -d),
      center+Point<3>( 0,  a, -a),
      center+Point<3>( 0, -a, -a),
      center+Point<3>( 0,  a,  a),
      center+Point<3>( 0, -a,  a),
      center+Point<3>( 0,  d,  d),
      center+Point<3>( 0, -d,  d),

      center+Point<3>(hc,  c, -c),
      center+Point<3>(hc, -c, -c),
      center+Point<3>(hb,  b, -b),
      center+Point<3>(hb, -b, -b),
      center+Point<3>(hb,  b,  b),
      center+Point<3>(hb, -b,  b),
      center+Point<3>(hc,  c,  c),
      center+Point<3>(hc, -c,  c),
    };

    int cell_vertices[6][8] =
    {
      {0, 1, 8, 9, 2, 3, 10, 11},
      {0, 2, 8, 10, 6, 4, 14, 12},
      {2, 3, 10, 11, 4, 5, 12, 13},
      {1, 7, 9, 15, 3, 5, 11, 13},
      {6, 4, 14, 12, 7, 5, 15, 13},
      {8, 10, 9, 11, 14, 12, 15, 13}
    };

    std::vector<CellData<3> > cells (6, CellData<3>());

    for (unsigned int i=0; i<6; ++i)
      {
        for (unsigned int j=0; j<8; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    tria.create_triangulation (
      std::vector<Point<3> >(&vertices[0], &vertices[16]),
      cells,
      SubCellData());       // no boundary information

    Triangulation<3>::cell_iterator cell = tria.begin();
    Triangulation<3>::cell_iterator end = tria.end();

    // go over all faces. for the ones on the flat face, set boundary
    // indicator for face and edges to one; the rest will remain at
    // zero but we have to pay attention to those edges that are
    // at the perimeter of the flat face since they should not be
    // set to one
    while (cell != end)
      {
        for (unsigned int i=0; i<GeometryInfo<3>::faces_per_cell; ++i)
          {
            if (!cell->at_boundary(i))
              continue;

            // If the center is on the plane x=0, this is a planar element. set
            // its boundary indicator. also set the boundary indicators of the
            // bounding faces unless both vertices are on the perimeter
            if (cell->face(i)->center()(0) < center(0)+1.e-5*radius)
              {
                cell->face(i)->set_boundary_indicator(1);
                for (unsigned int j=0; j<GeometryInfo<3>::lines_per_face; ++j)
                  {
                    const Point<3> line_vertices[2]
                      = { cell->face(i)->line(j)->vertex(0),
                          cell->face(i)->line(j)->vertex(1)
                        };
                    if ((std::fabs(line_vertices[0].distance(center)-radius) >
                         1e-5*radius)
                        ||
                        (std::fabs(line_vertices[1].distance(center)-radius) >
                         1e-5*radius))
                      cell->face(i)->line(j)->set_boundary_indicator(1);
                  }
              }
          }
        ++cell;
      }
  }


  template <>
  void
  hyper_shell (Triangulation<3> &tria,
               const Point<3> &p,
               const double inner_radius,
               const double outer_radius,
               const unsigned int n,
               const bool colorize)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());

    const double irad = inner_radius/std::sqrt(3.0);
    const double orad = outer_radius/std::sqrt(3.0);
    std::vector<Point<3> > vertices;
    std::vector<CellData<3> > cells;

    // Start with the shell bounded by
    // two nested cubes
    if (n == 6)
      {
        for (unsigned int i=0; i<8; ++i)
          vertices.push_back(p+hexahedron[i]*irad);
        for (unsigned int i=0; i<8; ++i)
          vertices.push_back(p+hexahedron[i]*orad);

        const unsigned int n_cells = 6;
        const int cell_vertices[n_cells][8] =
        {
          {8, 9, 10, 11, 0, 1, 2, 3}, // bottom
          {9, 11, 1, 3, 13, 15, 5, 7}, // right
          {12, 13, 4, 5, 14, 15, 6, 7}, // top
          {8, 0, 10, 2, 12, 4, 14, 6}, // left
          {8, 9, 0, 1, 12, 13, 4, 5}, // front
          {10, 2, 11, 3, 14, 6, 15, 7}
        }; // back

        cells.resize(n_cells, CellData<3>());

        for (unsigned int i=0; i<n_cells; ++i)
          {
            for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
              cells[i].vertices[j] = cell_vertices[i][j];
            cells[i].material_id = 0;
          }

        tria.create_triangulation (vertices, cells, SubCellData());
      }
    // A more regular subdivision can
    // be obtained by two nested
    // rhombic dodecahedra
    else if (n == 12)
      {
        for (unsigned int i=0; i<8; ++i)
          vertices.push_back(p+hexahedron[i]*irad);
        for (unsigned int i=0; i<6; ++i)
          vertices.push_back(p+octahedron[i]*inner_radius);
        for (unsigned int i=0; i<8; ++i)
          vertices.push_back(p+hexahedron[i]*orad);
        for (unsigned int i=0; i<6; ++i)
          vertices.push_back(p+octahedron[i]*outer_radius);

        const unsigned int n_cells = 12;
        const unsigned int rhombi[n_cells][4] =
        {
          { 10,  4,  0,  8},
          {  4, 13,  8,  6},
          { 10,  5,  4, 13},
          {  1,  9, 10,  5},
          {  9,  7,  5, 13},
          {  7, 11, 13,  6},
          {  9,  3,  7, 11},
          {  1, 12,  9,  3},
          { 12,  2,  3, 11},
          {  2,  8, 11,  6},
          { 12,  0,  2,  8},
          {  1, 10, 12,  0}
        };

        cells.resize(n_cells, CellData<3>());

        for (unsigned int i=0; i<n_cells; ++i)
          {
            for (unsigned int j=0; j<4; ++j)
              {
                cells[i].vertices[j  ] = rhombi[i][j];
                cells[i].vertices[j+4] = rhombi[i][j] + 14;
              }
            cells[i].material_id = 0;
          }

        tria.create_triangulation (vertices, cells, SubCellData());
      }
    else if (n == 96)
      {
        // create a triangulation based on the
        // 12-cell one where we refine the mesh
        // once and then re-arrange all
        // interior nodes so that the mesh is
        // the least distorted
        HyperShellBoundary<3> boundary (p);
        Triangulation<3> tmp;
        hyper_shell (tmp, p, inner_radius, outer_radius, 12);
        tmp.set_boundary(0, boundary);
        tmp.set_boundary(1, boundary);
        tmp.refine_global (1);

        // let's determine the distance at
        // which the interior nodes should be
        // from the center. let's say we
        // measure distances in multiples of
        // outer_radius and call
        // r=inner_radius.
        //
        // then note
        // that we now have 48 faces on the
        // inner and 48 on the outer sphere,
        // each with an area of approximately
        // 4*pi/48*r^2 and 4*pi/48, for
        // a face edge length of approximately
        // sqrt(pi/12)*r and sqrt(pi/12)
        //
        // let's say we put the interior nodes
        // at a distance rho, then a measure of
        // deformation for the inner cells
        // would be
        //   di=max(sqrt(pi/12)*r/(rho-r),
        //          (rho-r)/sqrt(pi/12)/r)
        // and for the outer cells
        //   do=max(sqrt(pi/12)/(1-rho),
        //          (1-rho)/sqrt(pi/12))
        //
        // we now seek a rho so that the
        // deformation of cells on the inside
        // and outside is equal. there are in
        // principle four possibilities for one
        // of the branches of do== one of the
        // branches of di, though not all of
        // them satisfy do==di, of
        // course. however, we are not
        // interested in cases where the inner
        // cell is long and skinny and the
        // outer one tall -- yes, they have the
        // same aspect ratio, but in different
        // space directions.
        //
        // so it only boils down to the
        // following two possibilities: the
        // first branch of each max(.,.)
        // functions are equal, or the second
        // one are. on the other hand, since
        // they two branches are reciprocals of
        // each other, if one pair of branches
        // is equal, so is the other
        //
        // this yields the following equation
        // for rho:
        //   sqrt(pi/12)*r/(rho-r)
        //   == sqrt(pi/12)/(1-rho)
        // with solution rho=2r/(1+r)
        const double r = inner_radius / outer_radius;
        const double rho = 2*r/(1+r);

        // then this is the distance of the
        // interior nodes from the center:
        const double middle_radius = rho * outer_radius;

        // mark vertices we've already moved or
        // that we want to ignore: we don't
        // want to move vertices at the inner
        // or outer boundaries
        std::vector<bool> vertex_already_treated (tmp.n_vertices(), false);
        for (Triangulation<3>::active_cell_iterator cell = tmp.begin_active();
             cell != tmp.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
              for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_face; ++v)
                vertex_already_treated[cell->face(f)->vertex_index(v)] = true;

        // now move the remaining vertices
        for (Triangulation<3>::active_cell_iterator cell = tmp.begin_active();
             cell != tmp.end(); ++cell)
          for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
            if (vertex_already_treated[cell->vertex_index(v)] == false)
              {
                // this is a new interior
                // vertex. mesh refinement may
                // have placed it at a number
                // of places in radial
                // direction and oftentimes not
                // in a particularly good
                // one. move it to halfway
                // between inner and outer
                // sphere
                const Point<3> old_distance = cell->vertex(v) - p;
                const double old_radius = cell->vertex(v).distance(p);
                cell->vertex(v) = p + old_distance * (middle_radius / old_radius);

                vertex_already_treated[cell->vertex_index(v)] = true;
              }

        // now copy the resulting level 1 cells
        // into the new triangulation,
        cells.resize(tmp.n_active_cells(), CellData<3>());

        unsigned int index = 0;
        for (Triangulation<3>::active_cell_iterator cell = tmp.begin_active();
             cell != tmp.end(); ++cell, ++index)
          {
            for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
              cells[index].vertices[v] = cell->vertex_index(v);
            cells[index].material_id = 0;
          }

        tria.create_triangulation (tmp.get_vertices(), cells, SubCellData());
      }
    else
      {
        Assert(false, ExcMessage ("Invalid number of coarse mesh cells."));
      }

    if (colorize)
      colorize_hyper_shell(tria, p, inner_radius, outer_radius);
  }




// Implementation for 3D only
  template <>
  void
  half_hyper_shell (Triangulation<3> &tria,
                    const Point<3> &center,
                    const double inner_radius,
                    const double outer_radius,
                    const unsigned int n,
                    const bool colorize)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());
    Assert(colorize == false, ExcNotImplemented());

    if (n <= 5)
      {
        // These are for the two lower squares
        const double d = outer_radius/std::sqrt(2.0);
        const double a = inner_radius/std::sqrt(2.0);
        // These are for the two upper square
        const double b = a/2.0;
        const double c = d/2.0;
        // And so are these
        const double hb = inner_radius*std::sqrt(3.0)/2.0;
        const double hc = outer_radius*std::sqrt(3.0)/2.0;

        Point<3> vertices[16] =
        {
          center+Point<3>( 0,  d, -d),
          center+Point<3>( 0, -d, -d),
          center+Point<3>( 0,  a, -a),
          center+Point<3>( 0, -a, -a),
          center+Point<3>( 0,  a,  a),
          center+Point<3>( 0, -a,  a),
          center+Point<3>( 0,  d,  d),
          center+Point<3>( 0, -d,  d),

          center+Point<3>(hc,  c, -c),
          center+Point<3>(hc, -c, -c),
          center+Point<3>(hb,  b, -b),
          center+Point<3>(hb, -b, -b),
          center+Point<3>(hb,  b,  b),
          center+Point<3>(hb, -b,  b),
          center+Point<3>(hc,  c,  c),
          center+Point<3>(hc, -c,  c),
        };

        int cell_vertices[5][8] =
        {
          {0, 1, 8, 9, 2, 3, 10, 11},
          {0, 2, 8, 10, 6, 4, 14, 12},
          {1, 7, 9, 15, 3, 5, 11, 13},
          {6, 4, 14, 12, 7, 5, 15, 13},
          {8, 10, 9, 11, 14, 12, 15, 13}
        };

        std::vector<CellData<3> > cells (5, CellData<3>());

        for (unsigned int i=0; i<5; ++i)
          {
            for (unsigned int j=0; j<8; ++j)
              cells[i].vertices[j] = cell_vertices[i][j];
            cells[i].material_id = 0;
          };

        tria.create_triangulation (
          std::vector<Point<3> >(&vertices[0], &vertices[16]),
          cells,
          SubCellData());       // no boundary information
      }
    else
      {
        Assert(false, ExcIndexRange(n, 0, 5));
      }

  }


// Implementation for 3D only
  template <>
  void quarter_hyper_shell (Triangulation<3> &tria,
                            const Point<3> &center,
                            const double inner_radius,
                            const double outer_radius,
                            const unsigned int n,
                            const bool colorize)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());
    if (n == 0 || n == 3)
      {
        const double a = inner_radius*std::sqrt(2.0)/2e0;
        const double b = outer_radius*std::sqrt(2.0)/2e0;
        const double c = a*std::sqrt(3.0)/2e0;
        const double d = b*std::sqrt(3.0)/2e0;
        const double e = outer_radius/2e0;
        const double h = inner_radius/2e0;

        std::vector<Point<3> > vertices;

        vertices.push_back (center+Point<3>( 0,  inner_radius, 0)); //0
        vertices.push_back (center+Point<3>( a,  a, 0));                  //1
        vertices.push_back (center+Point<3>( b,  b, 0));                  //2
        vertices.push_back (center+Point<3>( 0, outer_radius, 0));        //3
        vertices.push_back (center+Point<3>( 0, a , a));                  //4
        vertices.push_back (center+Point<3>( c, c , h));                  //5
        vertices.push_back (center+Point<3>( d, d , e));                  //6
        vertices.push_back (center+Point<3>( 0, b , b));                  //7
        vertices.push_back (center+Point<3>( inner_radius, 0 , 0));       //8
        vertices.push_back (center+Point<3>( outer_radius, 0 , 0));       //9
        vertices.push_back (center+Point<3>( a, 0 , a));  //10
        vertices.push_back (center+Point<3>( b, 0 , b));  //11
        vertices.push_back (center+Point<3>( 0, 0 , inner_radius));       //12
        vertices.push_back (center+Point<3>( 0, 0 , outer_radius));       //13

        const int cell_vertices[3][8] =
        {
          {0, 1, 3, 2, 4, 5, 7, 6},
          {1, 8, 2, 9, 5, 10, 6, 11},
          {4, 5, 7, 6, 12, 10, 13, 11},
        };
        std::vector<CellData<3> > cells(3);

        for (unsigned int i=0; i<3; ++i)
          {
            for (unsigned int j=0; j<8; ++j)
              cells[i].vertices[j] = cell_vertices[i][j];
            cells[i].material_id = 0;
          }

        tria.create_triangulation ( vertices, cells, SubCellData());       // no boundary information
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }

    if (colorize)
      colorize_quarter_hyper_shell(tria, center, inner_radius, outer_radius);
  }


// Implementation for 3D only
  template <>
  void cylinder_shell (Triangulation<3> &tria,
                       const double      length,
                       const double      inner_radius,
                       const double      outer_radius,
                       const unsigned int   n_radial_cells,
                       const unsigned int   n_axial_cells)
  {
    Assert ((inner_radius > 0) && (inner_radius < outer_radius),
            ExcInvalidRadii ());

    const double pi = numbers::PI;

    // determine the number of cells
    // for the grid. if not provided by
    // the user determine it such that
    // the length of each cell on the
    // median (in the middle between
    // the two circles) is equal to its
    // radial extent (which is the
    // difference between the two
    // radii)
    const unsigned int N_r = (n_radial_cells == 0 ?
                              static_cast<unsigned int>
                              (std::ceil((2*pi* (outer_radius + inner_radius)/2) /
                                         (outer_radius - inner_radius))) :
                              n_radial_cells);
    const unsigned int N_z = (n_axial_cells == 0 ?
                              static_cast<unsigned int>
                              (std::ceil (length /
                                          (2*pi*(outer_radius + inner_radius)/2/N_r))) :
                              n_axial_cells);

    // set up N vertices on the
    // outer and N vertices on
    // the inner circle. the
    // first N ones are on the
    // outer one, and all are
    // numbered counter-clockwise
    std::vector<Point<2> > vertices_2d(2*N_r);
    for (unsigned int i=0; i<N_r; ++i)
      {
        vertices_2d[i] = Point<2>( std::cos(2*pi * i/N_r),
                                   std::sin(2*pi * i/N_r)) * outer_radius;
        vertices_2d[i+N_r] = vertices_2d[i] * (inner_radius/outer_radius);
      };

    std::vector<Point<3> > vertices_3d;
    vertices_3d.reserve (2*N_r*(N_z+1));
    for (unsigned int j=0; j<=N_z; ++j)
      for (unsigned int i=0; i<2*N_r; ++i)
        {
          const Point<3> v (vertices_2d[i][0],
                            vertices_2d[i][1],
                            j*length/N_z);
          vertices_3d.push_back (v);
        }

    std::vector<CellData<3> > cells (N_r*N_z, CellData<3>());

    for (unsigned int j=0; j<N_z; ++j)
      for (unsigned int i=0; i<N_r; ++i)
        {
          cells[i+j*N_r].vertices[0] = i + (j+1)*2*N_r;
          cells[i+j*N_r].vertices[1] = (i+1)%N_r + (j+1)*2*N_r;
          cells[i+j*N_r].vertices[2] = i + j*2*N_r;
          cells[i+j*N_r].vertices[3] = (i+1)%N_r + j*2*N_r;

          cells[i+j*N_r].vertices[4] = N_r+i + (j+1)*2*N_r;
          cells[i+j*N_r].vertices[5] = N_r+((i+1)%N_r) + (j+1)*2*N_r;
          cells[i+j*N_r].vertices[6] = N_r+i + j*2*N_r;
          cells[i+j*N_r].vertices[7] = N_r+((i+1)%N_r) + j*2*N_r;

          cells[i+j*N_r].material_id = 0;
        }

    tria.create_triangulation (
      vertices_3d, cells, SubCellData());
  }



  template <int dim, int spacedim>
  void
  merge_triangulations (const Triangulation<dim, spacedim> &triangulation_1,
                        const Triangulation<dim, spacedim> &triangulation_2,
                        Triangulation<dim, spacedim>       &result)
  {
    Assert (triangulation_1.n_levels() == 1,
            ExcMessage ("The input triangulations must be coarse meshes."));
    Assert (triangulation_2.n_levels() == 1,
            ExcMessage ("The input triangulations must be coarse meshes."));

    // get the union of the set of vertices
    std::vector<Point<spacedim> > vertices = triangulation_1.get_vertices();
    vertices.insert (vertices.end(),
                     triangulation_2.get_vertices().begin(),
                     triangulation_2.get_vertices().end());

    // now form the union of the set of cells
    std::vector<CellData<dim> > cells;
    cells.reserve (triangulation_1.n_cells() + triangulation_2.n_cells());
    for (typename Triangulation<dim,spacedim>::cell_iterator
         cell = triangulation_1.begin(); cell != triangulation_1.end(); ++cell)
      {
        CellData<dim> this_cell;
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          this_cell.vertices[v] = cell->vertex_index(v);
        this_cell.material_id = cell->material_id();
        cells.push_back (this_cell);
      }

    // now do the same for the other other mesh. note that we have to
    // translate the vertex indices
    for (typename Triangulation<dim,spacedim>::cell_iterator
         cell = triangulation_2.begin(); cell != triangulation_2.end(); ++cell)
      {
        CellData<dim> this_cell;
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          this_cell.vertices[v] = cell->vertex_index(v) + triangulation_1.n_vertices();
        this_cell.material_id = cell->material_id();
        cells.push_back (this_cell);
      }

    // throw out duplicated vertices from the two meshes, reorder vertices as
    // necessary and create the triangulation
    SubCellData subcell_data;
    std::vector<unsigned int> considered_vertices;
    GridTools::delete_duplicated_vertices (vertices, cells,
                                           subcell_data,
                                           considered_vertices);

    // reorder the cells to ensure that they satisfy the convention for
    // edge and face directions
    GridReordering<dim, spacedim>::reorder_cells(cells, true);
    result.clear ();
    result.create_triangulation (vertices, cells, subcell_data);
  }


  template <int dim, int spacedim>
  void
  create_union_triangulation (const Triangulation<dim, spacedim> &triangulation_1,
                              const Triangulation<dim, spacedim> &triangulation_2,
                              Triangulation<dim, spacedim>       &result)
  {
    Assert (GridTools::have_same_coarse_mesh (triangulation_1, triangulation_2),
            ExcMessage ("The two input triangulations are not derived from "
                        "the same coarse mesh as required."));

    // first copy triangulation_1, and
    // then do as many iterations as
    // there are levels in
    // triangulation_2 to refine
    // additional cells. since this is
    // the maximum number of
    // refinements to get from the
    // coarse grid to triangulation_2,
    // it is clear that this is also
    // the maximum number of
    // refinements to get from any cell
    // on triangulation_1 to
    // triangulation_2
    result.clear ();
    result.copy_triangulation (triangulation_1);
    for (unsigned int iteration=0; iteration<triangulation_2.n_levels();
         ++iteration)
      {
        InterGridMap<Triangulation<dim, spacedim> > intergrid_map;
        intergrid_map.make_mapping (result, triangulation_2);

        bool any_cell_flagged = false;
        for (typename Triangulation<dim, spacedim>::active_cell_iterator
             result_cell = result.begin_active();
             result_cell != result.end(); ++result_cell)
          if (intergrid_map[result_cell]->has_children())
            {
              any_cell_flagged = true;
              result_cell->set_refine_flag ();
            }

        if (any_cell_flagged == false)
          break;
        else
          result.execute_coarsening_and_refinement();
      }
  }



  void
  extrude_triangulation(const Triangulation<2, 2> &input,
                        const unsigned int n_slices,
                        const double height,
                        Triangulation<3,3> &result)
  {
    Assert (input.n_levels() == 1,
            ExcMessage ("The input triangulations must be coarse meshes."));
    Assert(result.n_cells()==0, ExcMessage("resultin Triangulation need to be empty upon calling extrude_triangulation."));
    Assert(height>0, ExcMessage("The height in extrude_triangulation needs to be positive."));
    Assert(n_slices>=2, ExcMessage("The number of slices in extrude_triangulation needs to be at least 2."));

    std::vector<Point<3> > points(n_slices*input.n_vertices());
    std::vector<CellData<3> > cells;
    cells.reserve((n_slices-1)*input.n_active_cells());

    for (unsigned int slice=0; slice<n_slices; ++slice)
      {
        for (unsigned int i=0; i<input.n_vertices(); ++i)

          {
            const Point<2> &v = input.get_vertices()[i];
            points[i+slice*input.n_vertices()](0) = v(0);
            points[i+slice*input.n_vertices()](1) = v(1);
            points[i+slice*input.n_vertices()](2) = height * slice / (n_slices-1);
          }
      }

    for (Triangulation<2,2>::cell_iterator
         cell = input.begin(); cell != input.end(); ++cell)
      {
        for (unsigned int slice=0; slice<n_slices-1; ++slice)
          {
            CellData<3> this_cell;
            for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
              {
                this_cell.vertices[v]
                  = cell->vertex_index(v)+slice*input.n_vertices();
                this_cell.vertices[v+GeometryInfo<2>::vertices_per_cell]
                  = cell->vertex_index(v)+(slice+1)*input.n_vertices();
              }

            this_cell.material_id = cell->material_id();
            cells.push_back(this_cell);
          }
      }

    SubCellData s;
    types::boundary_id bid=0;
    s.boundary_quads.reserve(input.n_active_lines()*(n_slices-1) + input.n_active_cells()*2);
    for (Triangulation<2,2>::cell_iterator
         cell = input.begin(); cell != input.end(); ++cell)
      {
        CellData<2> quad;
        for (unsigned int f=0; f<4; ++f)
          if (cell->at_boundary(f))
            {
              quad.boundary_id = cell->face(f)->boundary_indicator();
              bid = std::max(bid, quad.boundary_id);
              for (unsigned int slice=0; slice<n_slices-1; ++slice)
                {
                  quad.vertices[0] = cell->face(f)->vertex_index(0)+slice*input.n_vertices();
                  quad.vertices[1] = cell->face(f)->vertex_index(1)+slice*input.n_vertices();
                  quad.vertices[2] = cell->face(f)->vertex_index(0)+(slice+1)*input.n_vertices();
                  quad.vertices[3] = cell->face(f)->vertex_index(1)+(slice+1)*input.n_vertices();
                  s.boundary_quads.push_back(quad);
                }
            }
      }

    for (Triangulation<2,2>::cell_iterator
         cell = input.begin(); cell != input.end(); ++cell)
      {
        CellData<2> quad;
        quad.boundary_id = bid + 1;
        quad.vertices[0] = cell->vertex_index(0);
        quad.vertices[1] = cell->vertex_index(1);
        quad.vertices[2] = cell->vertex_index(2);
        quad.vertices[3] = cell->vertex_index(3);
        s.boundary_quads.push_back(quad);

        quad.boundary_id = bid + 2;
        for (int i=0; i<4; ++i)
          quad.vertices[i] += (n_slices-1)*input.n_vertices();
        s.boundary_quads.push_back(quad);
      }

    result.create_triangulation (points,
                                 cells,
                                 s);
  }


  /**
   * Solve the Laplace equation for @p laplace_transformation function for one
   * of the @p dim space dimensions. Factorized into a function of its own
   * in order to allow parallel execution.
   */
  void laplace_solve (const SparseMatrix<double> &S,
                      const std::map<unsigned int,double> &m,
                      Vector<double> &u)
  {
    const unsigned int n_dofs=S.n();
    FilteredMatrix<Vector<double> > SF (S);
    PreconditionJacobi<SparseMatrix<double> > prec;
    prec.initialize(S, 1.2);
    FilteredMatrix<Vector<double> > PF (prec);

    SolverControl control (n_dofs, 1.e-10, false, false);
    GrowingVectorMemory<Vector<double> > mem;
    SolverCG<Vector<double> > solver (control, mem);

    Vector<double> f(n_dofs);

    SF.add_constraints(m);
    SF.apply_constraints (f, true);
    solver.solve(SF, u, f, PF);
  }


// Implementation for 1D only
  template <>
  void laplace_transformation (Triangulation<1> &,
                               const std::map<unsigned int,Point<1> > &,
                               const Function<1> *)
  {
    Assert(false, ExcNotImplemented());
  }


// Implementation for dimensions except 1
  template <int dim>
  void laplace_transformation (Triangulation<dim> &tria,
                               const std::map<unsigned int,Point<dim> > &new_points,
                               const Function<dim> *coefficient)
  {
    // first provide everything that is
    // needed for solving a Laplace
    // equation.
    MappingQ1<dim> mapping_q1;
    FE_Q<dim> q1(1);

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(q1);

    CompressedSparsityPattern c_sparsity_pattern (dof_handler.n_dofs (),
                                                  dof_handler.n_dofs ());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity_pattern);
    c_sparsity_pattern.compress ();

    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from (c_sparsity_pattern);
    sparsity_pattern.compress ();

    SparseMatrix<double> S(sparsity_pattern);

    QGauss<dim> quadrature(4);

    MatrixCreator::create_laplace_matrix(mapping_q1, dof_handler, quadrature, S,coefficient);

    // set up the boundary values for
    // the laplace problem
    std::vector<std::map<unsigned int,double> > m(dim);
    typename std::map<unsigned int,Point<dim> >::const_iterator map_end=new_points.end();

    // fill these maps using the data
    // given by new_points
    typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin_active(),
                                            endc=dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face=cell->face(face_no);

            // loop over all vertices of the cell and see if it is listed in the map
            // given as first argument of the function
            for (unsigned int vertex_no=0;
                 vertex_no<GeometryInfo<dim>::vertices_per_face; ++vertex_no)
              {
                const unsigned int vertex_index=face->vertex_index(vertex_no);

                const typename std::map<unsigned int,Point<dim> >::const_iterator map_iter
                  = new_points.find(vertex_index);

                if (map_iter!=map_end)
                  for (unsigned int i=0; i<dim; ++i)
                    m[i].insert(std::pair<unsigned int,double> (
                                  face->vertex_dof_index(vertex_no, 0), map_iter->second(i)));
              }
          }
      }

    // solve the dim problems with
    // different right hand sides.
    Vector<double> us[dim];
    for (unsigned int i=0; i<dim; ++i)
      us[i].reinit (dof_handler.n_dofs());

    // solve linear systems in parallel
    Threads::TaskGroup<> tasks;
    for (unsigned int i=0; i<dim; ++i)
      tasks += Threads::new_task (&laplace_solve,
                                  S, m[i], us[i]);
    tasks.join_all ();

    // change the coordinates of the
    // points of the triangulation
    // according to the computed values
    for (cell=dof_handler.begin_active(); cell!=endc; ++cell)
      for (unsigned int vertex_no=0;
           vertex_no<GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
        {
          Point<dim> &v=cell->vertex(vertex_no);
          const unsigned int dof_index=cell->vertex_dof_index(vertex_no, 0);
          for (unsigned int i=0; i<dim; ++i)
            v(i)=us[i](dof_index);
        }
  }



  template <>
  void hyper_cube_with_cylindrical_hole (Triangulation<1> &,
                                         const double,
                                         const double,
                                         const double,
                                         const unsigned int,
                                         bool)
  {
    Assert(false, ExcNotImplemented());
  }



  template <>
  void
  hyper_cube_with_cylindrical_hole (Triangulation<2> &triangulation,
                                    const double inner_radius,
                                    const double outer_radius,
                                    const double, // width,
                                    const unsigned int, // width_repetition,
                                    bool colorize)
  {
    const int dim = 2;

    Assert(inner_radius < outer_radius,
           ExcMessage("outer_radius has to be bigger than inner_radius."));

    Point<dim> center;
    // We create an hyper_shell in two dimensions, and then we modify it.
    hyper_shell (triangulation,
                 center, inner_radius, outer_radius,
                 8);
    Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    std::vector<bool> treated_vertices(triangulation.n_vertices(), false);
    for (; cell != endc; ++cell)
      {
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            {
              for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v)
                {
                  unsigned int vv = cell->face(f)->vertex_index(v);
                  if (treated_vertices[vv] == false)
                    {
                      treated_vertices[vv] = true;
                      switch (vv)
                        {
                        case 1:
                          cell->face(f)->vertex(v) = center+Point<dim>(outer_radius,outer_radius);
                          break;
                        case 3:
                          cell->face(f)->vertex(v) = center+Point<dim>(-outer_radius,outer_radius);
                          break;
                        case 5:
                          cell->face(f)->vertex(v) = center+Point<dim>(-outer_radius,-outer_radius);
                          break;
                        case 7:
                          cell->face(f)->vertex(v) = center+Point<dim>(outer_radius,-outer_radius);
                        default:
                          break;
                        }
                    }
                }
            }
      }
    double eps = 1e-3 * outer_radius;
    cell = triangulation.begin_active();
    for (; cell != endc; ++cell)
      {
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            {
              double dx = cell->face(f)->center()(0) - center(0);
              double dy = cell->face(f)->center()(1) - center(1);
              if (colorize)
                {
                  if (std::abs(dx + outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(0);
                  else if (std::abs(dx - outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(1);
                  else if (std::abs(dy + outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(2);
                  else if (std::abs(dy - outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(3);
                  else
                    cell->face(f)->set_boundary_indicator(4);
                }
              else
                {
                  double d = (cell->face(f)->center() - center).norm();
                  if (d-inner_radius < 0)
                    cell->face(f)->set_boundary_indicator(1);
                  else
                    cell->face(f)->set_boundary_indicator(0);
                }
            }
      }
  }



  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<3> &triangulation,
                                        const double inner_radius,
                                        const double outer_radius,
                                        const double L,
                                        const unsigned int Nz,
                                        bool colorize)
  {
    const int dim = 3;

    Assert(inner_radius < outer_radius,
           ExcMessage("outer_radius has to be bigger than inner_radius."));
    Assert(L > 0,
           ExcMessage("Must give positive extension L"));
    Assert(Nz >= 1, ExcLowerRange(1, Nz));

    cylinder_shell (triangulation,
                    L, inner_radius, outer_radius,
                    8,
                    Nz);

    Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    std::vector<bool> treated_vertices(triangulation.n_vertices(), false);
    for (; cell != endc; ++cell)
      {
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            {
              for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v)
                {
                  unsigned int vv = cell->face(f)->vertex_index(v);
                  if (treated_vertices[vv] == false)
                    {
                      treated_vertices[vv] = true;
                      for (unsigned int i=0; i<=Nz; ++i)
                        {
                          double d = ((double) i)*L/((double) Nz);
                          switch (vv-i*16)
                            {
                            case 1:
                              cell->face(f)->vertex(v) = Point<dim>(outer_radius,outer_radius,d);
                              break;
                            case 3:
                              cell->face(f)->vertex(v) = Point<dim>(-outer_radius,outer_radius,d);
                              break;
                            case 5:
                              cell->face(f)->vertex(v) = Point<dim>(-outer_radius,-outer_radius,d);
                              break;
                            case 7:
                              cell->face(f)->vertex(v) = Point<dim>(outer_radius,-outer_radius,d);
                              break;
                            default:
                              break;
                            }
                        }
                    }
                }
            }
      }
    double eps = 1e-3 * outer_radius;
    cell = triangulation.begin_active();
    for (; cell != endc; ++cell)
      {
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            {
              double dx = cell->face(f)->center()(0);
              double dy = cell->face(f)->center()(1);
              double dz = cell->face(f)->center()(2);

              if (colorize)
                {
                  if (std::abs(dx + outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(0);

                  else if (std::abs(dx - outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(1);

                  else if (std::abs(dy + outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(2);

                  else if (std::abs(dy - outer_radius) < eps)
                    cell->face(f)->set_boundary_indicator(3);

                  else if (std::abs(dz) < eps)
                    cell->face(f)->set_boundary_indicator(4);

                  else if (std::abs(dz - L) < eps)
                    cell->face(f)->set_boundary_indicator(5);

                  else
                    {
                      cell->face(f)->set_boundary_indicator(6);
                      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
                        cell->face(f)->line(l)->set_boundary_indicator(6);
                    }

                }
              else
                {
                  Point<dim> c = cell->face(f)->center();
                  c(2) = 0;
                  double d = c.norm();
                  if (d-inner_radius < 0)
                    {
                      cell->face(f)->set_boundary_indicator(1);
                      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
                        cell->face(f)->line(l)->set_boundary_indicator(1);
                    }
                  else
                    cell->face(f)->set_boundary_indicator(0);
                }
            }
      }
  }

  template <int dim, int spacedim1, int spacedim2>
  void flatten_triangulation(const Triangulation<dim, spacedim1> &in_tria,
                             Triangulation<dim,spacedim2> &out_tria)
  {
    const parallel::distributed::Triangulation<dim, spacedim1> *pt =
      dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim1> *>(&in_tria);

    Assert (pt == NULL,
            ExcMessage("Cannot use this function on parallel::distributed::Triangulation."));

    std::vector<Point<spacedim2> > v;
    std::vector<CellData<dim> > cells;
    SubCellData subcelldata;

    const unsigned int spacedim = std::min(spacedim1,spacedim2);
    const std::vector<Point<spacedim1> > &in_vertices = in_tria.get_vertices();

    v.resize(in_vertices.size());
    for (unsigned int i=0; i<in_vertices.size(); ++i)
      for (unsigned int d=0; d<spacedim; ++d)
        v[i][d] = in_vertices[i][d];

    cells.resize(in_tria.n_active_cells());
    typename Triangulation<dim,spacedim1>::active_cell_iterator
    cell = in_tria.begin_active(),
    endc = in_tria.end();

    for (unsigned int id=0; cell != endc; ++cell, ++id)
      {
        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          cells[id].vertices[i] = cell->vertex_index(i);
        cells[id].material_id = cell->material_id();
        cells[id].manifold_id = cell->manifold_id();
      }

    if (dim>1)
      {
        typename Triangulation<dim,spacedim1>::active_face_iterator
        face = in_tria.begin_active_face(),
        endf = in_tria.end_face();

        // Face counter for both dim == 2 and dim == 3
        unsigned int f=0;
        switch (dim)
          {
          case 2:
          {
            subcelldata.boundary_lines.resize(in_tria.n_active_faces());
            for (; face != endf; ++face)
              if (face->at_boundary())
                {
                  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
                    subcelldata.boundary_lines[f].vertices[i] = face->vertex_index(i);
                  subcelldata.boundary_lines[f].boundary_id = face->boundary_indicator();
                  subcelldata.boundary_lines[f].manifold_id = face->manifold_id();
                  ++f;
                }
            subcelldata.boundary_lines.resize(f);
          }
          break;
          case 3:
          {
            subcelldata.boundary_quads.resize(in_tria.n_active_faces());
            for (; face != endf; ++face)
              if (face->at_boundary())
                {
                  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
                    subcelldata.boundary_quads[f].vertices[i] = face->vertex_index(i);
                  subcelldata.boundary_quads[f].boundary_id = face->boundary_indicator();
                  subcelldata.boundary_quads[f].manifold_id = face->manifold_id();
                  ++f;
                }
            subcelldata.boundary_quads.resize(f);
          }
          break;
          default:
            Assert(false, ExcInternalError());
          }
      }
    out_tria.create_triangulation(v, cells, subcelldata);
  }



  // This anonymous namespace contains utility functions to extract the
  // triangulation from any container such as DoFHandler
  // and the like
  namespace
  {
    template<int dim, int spacedim>
    const Triangulation<dim, spacedim> &
    get_tria(const Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, int spacedim>
    const Triangulation<dim, spacedim> &
    get_tria(const parallel::distributed::Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, template<int, int> class Container, int spacedim>
    const Triangulation<dim,spacedim> &
    get_tria(const Container<dim,spacedim> &container)
    {
      return container.get_tria();
    }


    template<int dim, int spacedim>
    Triangulation<dim, spacedim> &
    get_tria(Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, int spacedim>
    Triangulation<dim, spacedim> &
    get_tria(parallel::distributed::Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, template<int, int> class Container, int spacedim>
    const Triangulation<dim,spacedim> &
    get_tria(Container<dim,spacedim> &container)
    {
      return container.get_tria();
    }
  }



  template <template <int,int> class Container, int dim, int spacedim>
#ifndef _MSC_VER
  std::map<typename Container<dim-1,spacedim>::cell_iterator,
      typename Container<dim,spacedim>::face_iterator>
#else
  typename ExtractBoundaryMesh<Container,dim,spacedim>::return_type
#endif
      extract_boundary_mesh (const Container<dim,spacedim> &volume_mesh,
                             Container<dim-1,spacedim>     &surface_mesh,
                             const std::set<types::boundary_id> &boundary_ids)
  {
// This function works using the following assumption:
//    Triangulation::create_triangulation(...) will create cells that preserve
//    the order of cells passed in using the CellData argument; also,
//    that it will not reorder the vertices.

    std::map<typename Container<dim-1,spacedim>::cell_iterator,
        typename Container<dim,spacedim>::face_iterator>
        surface_to_volume_mapping;

    const unsigned int boundary_dim = dim-1; //dimension of the boundary mesh

    // First create surface mesh and mapping
    // from only level(0) cells of volume_mesh
    std::vector<typename Container<dim,spacedim>::face_iterator>
    mapping;  // temporary map for level==0


    std::vector< bool > touched (get_tria(volume_mesh).n_vertices(), false);
    std::vector< CellData< boundary_dim > > cells;
    std::vector< Point<spacedim> >      vertices;

    std::map<unsigned int,unsigned int> map_vert_index; //volume vertex indices to surf ones

    unsigned int v_index;
    CellData< boundary_dim > c_data;

    for (typename Container<dim,spacedim>::cell_iterator
         cell = volume_mesh.begin(0);
         cell != volume_mesh.end(0);
         ++cell)
      for (unsigned int i=0; i < GeometryInfo<dim>::faces_per_cell; ++i)
        {
          const typename Container<dim,spacedim>::face_iterator
          face = cell->face(i);

          if ( face->at_boundary()
               &&
               (boundary_ids.empty() ||
                ( boundary_ids.find(face->boundary_indicator()) != boundary_ids.end())) )
            {
              for (unsigned int j=0;
                   j<GeometryInfo<boundary_dim>::vertices_per_cell; ++j)
                {
                  v_index = face->vertex_index(j);

                  if ( !touched[v_index] )
                    {
                      vertices.push_back(face->vertex(j));
                      map_vert_index[v_index] = vertices.size() - 1;
                      touched[v_index] = true;
                    }

                  c_data.vertices[j] = map_vert_index[v_index];
                  c_data.material_id = static_cast<types::material_id>(face->boundary_indicator());
                }

              cells.push_back(c_data);
              mapping.push_back(face);
            }
        }

    // create level 0 surface triangulation
    Assert (cells.size() > 0, ExcMessage ("No boundary faces selected"));
    const_cast<Triangulation<dim-1,spacedim>&>(get_tria(surface_mesh))
    .create_triangulation (vertices, cells, SubCellData());

    // Make the actual mapping
    for (typename Container<dim-1,spacedim>::active_cell_iterator
         cell = surface_mesh.begin(0);
         cell!=surface_mesh.end(0); ++cell)
      surface_to_volume_mapping[cell] = mapping.at(cell->index());

    do
      {
        bool changed = false;

        for (typename Container<dim-1,spacedim>::active_cell_iterator
             cell = surface_mesh.begin_active(); cell!=surface_mesh.end(); ++cell)
          if (surface_to_volume_mapping[cell]->has_children() == true )
            {
              cell->set_refine_flag ();
              changed = true;
            }

        if (changed)
          {
            const_cast<Triangulation<dim-1,spacedim>&>(get_tria(surface_mesh))
            .execute_coarsening_and_refinement();

            for (typename Container<dim-1,spacedim>::cell_iterator
                 surface_cell = surface_mesh.begin(); surface_cell!=surface_mesh.end(); ++surface_cell)
              for (unsigned int c=0; c<surface_cell->n_children(); c++)
                if (surface_to_volume_mapping.find(surface_cell->child(c)) == surface_to_volume_mapping.end())
                  surface_to_volume_mapping[surface_cell->child(c)]
                    = surface_to_volume_mapping[surface_cell]->child(c);
          }
        else
          break;
      }
    while (true);

    return surface_to_volume_mapping;
  }

}

// explicit instantiations
#include "grid_generator.inst"

DEAL_II_NAMESPACE_CLOSE

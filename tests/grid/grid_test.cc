// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/numbers.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


// 1: continuous refinement of the unit square always in the middle
// 2: refinement of the circle at the boundary
// 2: refinement of a wiggled area at the boundary

template <int dim>
class Ball : public FlatManifold<dim>
{
public:
  virtual Point<dim>
  get_new_point(const ArrayView<const Point<dim>> &surrounding_points,
                const ArrayView<const double>     &weights) const override
  {
    Point<dim> middle =
      FlatManifold<dim>::get_new_point(surrounding_points, weights);

    for (int i = 0; i < dim; ++i)
      middle[i] -= .5;
    middle *=
      std::sqrt(static_cast<double>(dim)) / (std::sqrt(middle.square()) * 2);
    for (int i = 0; i < dim; ++i)
      middle[i] += .5;

    return middle;
  }

  virtual std::unique_ptr<Manifold<dim, dim>>
  clone() const override
  {
    return std::unique_ptr<Manifold<dim, dim>>(new Ball<dim>());
  }
};


template <int dim>
class CurvedLine : public FlatManifold<dim>
{
public:
  virtual Point<dim>
  get_new_point_on_line(
    const typename Triangulation<dim>::line_iterator &line) const;

  virtual Point<dim>
  get_new_point_on_quad(
    const typename Triangulation<dim>::quad_iterator &quad) const;
};


template <int dim>
Point<dim>
CurvedLine<dim>::get_new_point_on_line(
  const typename Triangulation<dim>::line_iterator &line) const
{
  Point<dim> middle = FlatManifold<dim>::get_new_point_on_line(line);

  // if the line is at the top or bottom
  // face: do a special treatment on
  // this line. Note that if the
  // z-value of the midpoint is either
  // 0 or 1, then the z-values of all
  // vertices of the line is like that
  if (dim >= 3)
    if (((middle[2] == 0) || (middle[2] == 1))
        // find out, if the line is in the
        // interior of the top or bottom face
        // of the domain, or at the edge.
        // lines at the edge need to undergo
        // the usual treatment, while for
        // interior lines taking the midpoint
        // is sufficient
        //
        // note: the trick with the boundary
        // id was invented after the above was
        // written, so we are not very strict
        // here with using these flags
        && (line->manifold_id() == 1))
      return middle;


  double x = middle[0], y = middle[1];

  if (y < x)
    if (y < 1 - x)
      middle[1] = 0.04 * std::sin(6 * numbers::PI * middle[0]);
    else
      middle[0] = 1 + 0.04 * std::sin(6 * numbers::PI * middle[1]);

  else if (y < 1 - x)
    middle[0] = 0.04 * std::sin(6 * numbers::PI * middle[1]);
  else
    middle[1] = 1 + 0.04 * std::sin(6 * numbers::PI * middle[0]);

  return middle;
}


template <int dim>
Point<dim>
CurvedLine<dim>::get_new_point_on_quad(
  const typename Triangulation<dim>::quad_iterator &quad) const
{
  Point<dim> middle = FlatManifold<dim>::get_new_point_on_quad(quad);

  // if the face is at the top of bottom
  // face: do not move the midpoint in
  // x/y direction. Note that if the
  // z-value of the midpoint is either
  // 0 or 1, then the z-values of all
  // vertices of the quad is like that
  if (dim == 3 && ((middle[dim - 1] == 0) || (middle[dim - 1] == 1)))
    return middle;

  double x = middle[0], y = middle[1];

  if (y < x)
    if (y < 1 - x)
      middle[1] = 0.04 * std::sin(6 * numbers::PI * middle[0]);
    else
      middle[0] = 1 + 0.04 * std::sin(6 * numbers::PI * middle[1]);

  else if (y < 1 - x)
    middle[0] = 0.04 * std::sin(6 * numbers::PI * middle[1]);
  else
    middle[1] = 1 + 0.04 * std::sin(6 * numbers::PI * middle[0]);

  return middle;
}


template <int dim>
void
test(const int test_case)
{
  char testname[100];
  sprintf(testname, "Test%d.dim%d", test_case, dim);

  deallog.push(testname);
  deallog << "Start" << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_material_id(3);

  if ((dim == 1) && ((test_case == 2) || (test_case == 3)))
    {
      deallog << "Impossible for this dimension." << std::endl;
      return;
    };


  switch (test_case)
    {
      case 1:
        {
          // refine first cell
          tria.begin_active()->set_refine_flag();
          tria.execute_coarsening_and_refinement();

          // refine first active cell
          // on coarsest level
          tria.begin_active()->set_refine_flag();
          tria.execute_coarsening_and_refinement();

          typename Triangulation<dim>::active_cell_iterator cell;
          for (int i = 0; i < (dim == 2 ? 3 : 2); ++i)
            {
              // refine the presently
              // last cell several
              // times
              cell = tria.last_active();
              cell->set_refine_flag();
              tria.execute_coarsening_and_refinement();
            };

          break;
        }

      case 2:
      case 3:
        {
          if (dim == 3)
            {
              tria.begin_active()->face(4)->set_manifold_id(1);
              tria.begin_active()->face(5)->set_manifold_id(1);
            };


          // set the manifold function
          Ball<dim>       ball;
          CurvedLine<dim> curved_line;
          if (test_case == 2)
            tria.set_manifold(1, ball);
          else
            tria.set_manifold(1, curved_line);

          // refine once
          tria.begin_active()->set_refine_flag();
          tria.execute_coarsening_and_refinement();

          typename Triangulation<dim>::active_cell_iterator cell, endc;
          const unsigned int steps[4] = {0, 2, 2, 2};
          for (unsigned int i = 0; i < steps[dim]; ++i)
            {
              cell = tria.begin_active();
              endc = tria.end();

              // refine all
              // boundary cells
              for (; cell != endc; ++cell)
                if (cell->at_boundary())
                  cell->set_refine_flag();

              tria.execute_coarsening_and_refinement();
            };

          tria.reset_manifold(1);
          break;
        }
    }

  GridOut go;
  go.set_flags(GridOutFlags::Ucd(true));
  go.write_ucd(tria, deallog.get_file_stream());

  deallog << "     Total number of cells        = " << tria.n_cells()
          << std::endl
          << "     Total number of active cells = " << tria.n_active_cells()
          << std::endl;

  deallog.pop();
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(8);

  for (unsigned int i = 1; i <= 3; ++i)
    test<2>(i);
  for (unsigned int i = 1; i <= 3; ++i)
    test<3>(i);

  return 0;
}

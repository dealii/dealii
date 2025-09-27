#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

template <int dim>
void
test(const ReferenceCell &reference_cell,
     const Point<dim>    &p,
     const bool           print = true)
{
  const Point<dim> projected       = reference_cell.closest_point(p);
  const double     distance_square = projected.distance_square(p);
  // its difficult to place points exactly on sloped surfaces (e.g., for
  // pyramids), but we should always be close
  AssertThrow(reference_cell.contains_point(
                projected, 2.0 * std::numeric_limits<double>::epsilon()),
              ExcInternalError());
  if (print)
    deallog << "point = " << p << std::endl
            << "  d^2 = " << distance_square << std::endl;

  // The distance we get should be at least as good as the distance to the
  // nearest vertex
  for (const unsigned int vertex_no : reference_cell.vertex_indices())
    AssertThrow(distance_square <=
                  p.distance_square(
                    reference_cell.template vertex<dim>(vertex_no)) *
                    (1. + 2. * std::numeric_limits<double>::epsilon()),
                ExcInternalError());

  // Also generate a lot of points in the unit cell - we should be closer than
  // any generated point
  const unsigned int n_points_1d = 100;
  if (dim == 2)
    {
      for (unsigned int i = 0; i < n_points_1d; ++i)
        for (unsigned int j = 0; j < n_points_1d; ++j)
          {
            const Point<dim> grid_point(i / double(n_points_1d),
                                        j / double(n_points_1d));
            if (reference_cell.contains_point(grid_point))
              AssertThrow(distance_square <= grid_point.distance_square(p),
                          ExcInternalError());
          }
    }
  else if (dim == 3)
    {
      for (unsigned int i = 0; i < n_points_1d; ++i)
        for (unsigned int j = 0; j < n_points_1d; ++j)
          for (unsigned int k = 0; j < n_points_1d; ++j)
            {
              // scale things so that we cover the bounding box of both pyramids
              // and wedges
              const Point<dim> grid_point(-1.0 + 2.0 * i / double(n_points_1d),
                                          -1.0 + 2.0 * j / double(n_points_1d),
                                          k / double(n_points_1d));
              if (reference_cell.contains_point(grid_point))
                AssertThrow(distance_square <= grid_point.distance_square(p),
                            ExcInternalError());
            }
    }
}

int
main()
{
  initlog();

  for (const auto &reference_cell : {ReferenceCells::Line})
    {
      deallog << "reference cell = " << reference_cell.to_string() << std::endl;
      test(reference_cell, Point<1>(-1));
      test(reference_cell, Point<1>(-0.5));
      test(reference_cell, Point<1>(0.5));
      test(reference_cell, Point<1>(1.0));
      test(reference_cell, Point<1>(1.25));
    }

  for (const auto &reference_cell :
       {ReferenceCells::Triangle, ReferenceCells::Quadrilateral})
    {
      deallog << "reference cell = " << reference_cell.to_string() << std::endl;
      test(reference_cell, Point<2>(-0.5, 0.5));
      test(reference_cell, Point<2>(-1, -1));
      test(reference_cell, Point<2>(0.5, -0.5));
      test(reference_cell, Point<2>(1.0, 0.0));
      test(reference_cell, Point<2>(1.0, 1.0));
      test(reference_cell, Point<2>(0.25, 0.25));

      for (unsigned int i = 0; i < 1e3; ++i)
        test(reference_cell,
             Point<2>(random_value(-4.0, 4.0), random_value(-4.0, 4.0)),
             false);
    }

  for (const auto &reference_cell : {ReferenceCells::Tetrahedron,
                                     ReferenceCells::Pyramid,
                                     ReferenceCells::Wedge,
                                     ReferenceCells::Hexahedron})
    {
      deallog << "reference cell = " << reference_cell.to_string() << std::endl;
      test(reference_cell, Point<3>(-0.5, 0.5, 0.5));
      test(reference_cell, Point<3>(-1, -1, -1));
      test(reference_cell, Point<3>(0.5, -0.5, -0.5));
      test(reference_cell, Point<3>(1.0, 0.0, 0.0));
      test(reference_cell, Point<3>(1.0, 1.0, 1.0));
      test(reference_cell, Point<3>(-0.75, 0.5, 1.5));
      test(reference_cell, Point<3>(2.0, 2.0, 1.0));
      test(reference_cell, Point<3>(0.25, 0.25, 0.25));
      test(reference_cell, Point<3>(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
      test(reference_cell, Point<3>(0.7, 0.9, 0.1));
      test(reference_cell, Point<3>(0.08, 0.18, -0.07));
      test(reference_cell, Point<3>(0.0, 0.0, 0.0));

      for (unsigned int i = 0; i < 1e3; ++i)
        test(reference_cell,
             Point<3>(random_value(-4.0, 4.0),
                      random_value(-4.0, 4.0),
                      random_value(-4.0, 4.0)),
             false);
    }
  deallog << "OK!" << std::endl;
}

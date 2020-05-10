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

// Short test to validate GridTools::compute_maximum_aspect_ratio()

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


/*
 * The parameters left, right, refinements follow the hyper-rectangle
 * notation.
 *
 * degree denotes the polynomial degree of the mapping.
 *
 * n_q_points denotes the number of 1D quadrature points used during
 * the computation of the aspect ratio.
 *
 * If deform is true, the first vertex of the first element of the
 * hyper_rectangle is shifted, with factor being a scaling factor to
 * be able to generate inverted elements.
 */
template <int dim>
double
compute_aspect_ratio_hyper_rectangle(
  Point<dim> const &               left,
  Point<dim> const &               right,
  std::vector<unsigned int> const &refinements,
  unsigned int                     degree     = 1,
  unsigned int                     n_q_points = 2,
  bool                             deform     = false,
  double                           factor     = 1.0)
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_rectangle(
    tria, refinements, left, right, false);

  if (deform)
    {
      Tensor<1, dim> diag = right - left;
      double         l    = diag.norm();
      Point<dim>     shift;
      for (unsigned int d = 0; d < dim; ++d)
        shift[d] = l * factor * (0.05 + d * 0.01);
      tria.begin_active()->vertex(0) += shift;
    }

  MappingQGeneric<dim> const mapping(degree);
  QGauss<dim> const          gauss(n_q_points);

  Vector<double> ratios =
    GridTools::compute_aspect_ratio_of_cells(mapping, tria, gauss);

  deallog << std::endl
          << "Parameters:"
          << " d = " << dim << ", degree = " << degree
          << ", n_q_points = " << n_q_points << ", deform = " << deform
          << ", factor = " << factor << std::endl;

  deallog << "Aspect ratio vector = " << std::endl;
  ratios.print(deallog.get_file_stream());

  return GridTools::compute_maximum_aspect_ratio(mapping, tria, gauss);
}

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(6);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

      // 1D
      {
        Point<1> left  = Point<1>(0.0);
        Point<1> right = Point<1>(1.0);

        double ar = 0.0;

        std::vector<unsigned int> refine(1, 2);

        ar = compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        refine[0] = 5;

        ar = compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;
      }

      // 2D
      {
        Point<2> left  = Point<2>(0.0, 0.0);
        Point<2> right = Point<2>(1.0, 1.0);

        double ar = 0.0;

        std::vector<unsigned int> refine(2, 2);

        ar = compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar =
          compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar = compute_aspect_ratio_hyper_rectangle(
          left, right, refine, 1, 2, true, 10);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        refine[0] = 5;

        ar = compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar =
          compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar =
          compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 5, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar = compute_aspect_ratio_hyper_rectangle(
          left, right, refine, 1, 15, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar = compute_aspect_ratio_hyper_rectangle(
          left, right, refine, 1, 2, true, 10);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;
      }

      // 3D
      {
        Point<3> left  = Point<3>(0.0, 0.0, 0.0);
        Point<3> right = Point<3>(1.0, 1.0, 1.0);

        double ar = 0.0;

        std::vector<unsigned int> refine(3, 2);

        ar = compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar =
          compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar = compute_aspect_ratio_hyper_rectangle(
          left, right, refine, 1, 2, true, 10);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        refine[0] = 5;
        refine[1] = 3;

        ar = compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar =
          compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 2, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar =
          compute_aspect_ratio_hyper_rectangle(left, right, refine, 1, 5, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar = compute_aspect_ratio_hyper_rectangle(
          left, right, refine, 1, 15, true);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;

        ar = compute_aspect_ratio_hyper_rectangle(
          left, right, refine, 1, 2, true, 10);
        deallog << "aspect ratio max    = " << ar << std::endl << std::endl;
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}

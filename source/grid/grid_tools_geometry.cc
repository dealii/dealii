// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools_geometry.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/numerics/vector_tools_integrate_difference.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN


namespace GridTools
{
  template <int dim, int spacedim>
  double
  diameter(const Triangulation<dim, spacedim> &tria)
  {
    // we can't deal with distributed meshes since we don't have all
    // vertices locally. there is one exception, however: if the mesh has
    // never been refined. the way to test this is not to ask
    // tria.n_levels()==1, since this is something that can happen on one
    // processor without being true on all. however, we can ask for the
    // global number of active cells and use that
    if constexpr (running_in_debug_mode())
      {
        if (const auto *p_tria = dynamic_cast<
              const parallel::DistributedTriangulationBase<dim, spacedim> *>(
              &tria))
          Assert(p_tria->n_global_active_cells() == tria.n_cells(0),
                 ExcNotImplemented());
      }

    // the algorithm used simply traverses all cells and picks out the
    // boundary vertices. it may or may not be faster to simply get all
    // vectors, don't mark boundary vertices, and compute the distances
    // thereof, but at least as the mesh is refined, it seems better to
    // first mark boundary nodes, as marking is O(N) in the number of
    // cells/vertices, while computing the maximal distance is O(N*N)
    const std::vector<Point<spacedim>> &vertices = tria.get_vertices();
    std::vector<bool> boundary_vertices(vertices.size(), false);

    typename Triangulation<dim, spacedim>::active_cell_iterator cell =
      tria.begin_active();
    const typename Triangulation<dim, spacedim>::active_cell_iterator endc =
      tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int face : cell->face_indices())
        if (cell->face(face)->at_boundary())
          for (unsigned int i = 0; i < cell->face(face)->n_vertices(); ++i)
            boundary_vertices[cell->face(face)->vertex_index(i)] = true;

    // now traverse the list of boundary vertices and check distances.
    // since distances are symmetric, we only have to check one half
    double                            max_distance_sqr = 0;
    std::vector<bool>::const_iterator pi = boundary_vertices.begin();
    const unsigned int                N  = boundary_vertices.size();
    for (unsigned int i = 0; i < N; ++i, ++pi)
      {
        std::vector<bool>::const_iterator pj = pi + 1;
        for (unsigned int j = i + 1; j < N; ++j, ++pj)
          if ((*pi == true) && (*pj == true) &&
              ((vertices[i] - vertices[j]).norm_square() > max_distance_sqr))
            max_distance_sqr = (vertices[i] - vertices[j]).norm_square();
      }

    return std::sqrt(max_distance_sqr);
  }



  template <int dim, int spacedim>
  double
  volume(const Triangulation<dim, spacedim> &triangulation)
  {
    Assert(triangulation.get_reference_cells().size() == 1,
           ExcNotImplemented());
    const ReferenceCell reference_cell = triangulation.get_reference_cells()[0];
    return volume(
      triangulation,
      reference_cell.template get_default_linear_mapping<dim, spacedim>());
  }



  template <int dim, int spacedim>
  double
  volume(const Triangulation<dim, spacedim> &triangulation,
         const Mapping<dim, spacedim>       &mapping)
  {
    // get the degree of the mapping if possible. if not, just assume 1
    unsigned int mapping_degree = 1;
    if (const auto *p = dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping))
      mapping_degree = p->get_degree();
    else if (const auto *p =
               dynamic_cast<const MappingFE<dim, spacedim> *>(&mapping))
      mapping_degree = p->get_degree();

    // then initialize an appropriate quadrature formula
    Assert(triangulation.get_reference_cells().size() == 1,
           ExcNotImplemented());
    const ReferenceCell reference_cell = triangulation.get_reference_cells()[0];
    const Quadrature<dim> quadrature_formula =
      reference_cell.template get_gauss_type_quadrature<dim>(mapping_degree +
                                                             1);
    const unsigned int n_q_points = quadrature_formula.size();

    // we really want the JxW values from the FEValues object, but it
    // wants a finite element. create a cheap element as a dummy
    // element
    FE_Nothing<dim, spacedim> dummy_fe(reference_cell);
    FEValues<dim, spacedim>   fe_values(mapping,
                                      dummy_fe,
                                      quadrature_formula,
                                      update_JxW_values);

    double local_volume = 0;

    // compute the integral quantities by quadrature
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          for (unsigned int q = 0; q < n_q_points; ++q)
            local_volume += fe_values.JxW(q);
        }

    const double global_volume =
      Utilities::MPI::sum(local_volume, triangulation.get_mpi_communicator());

    return global_volume;
  }



  template <int dim, int spacedim>
  std::pair<unsigned int, double>
  get_longest_direction(
    typename Triangulation<dim, spacedim>::active_cell_iterator cell)
  {
    double       max_ratio = 1;
    unsigned int index     = 0;

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i + 1; j < dim; ++j)
        {
          unsigned int ax      = i % dim;
          unsigned int next_ax = j % dim;

          double ratio =
            cell->extent_in_direction(ax) / cell->extent_in_direction(next_ax);

          if (ratio > max_ratio)
            {
              max_ratio = ratio;
              index     = ax;
            }
          else if (1.0 / ratio > max_ratio)
            {
              max_ratio = 1.0 / ratio;
              index     = next_ax;
            }
        }
    return std::make_pair(index, max_ratio);
  }



  namespace
  {
    /**
     * The algorithm to compute the affine approximation to a cell finds an
     * affine map A x_hat + b from the reference cell to the real space.
     *
     * Some details about how we compute the least square plane. We look
     * for a spacedim x (dim + 1) matrix X such that X * M = Y where M is
     * a (dim+1) x n_vertices matrix and Y a spacedim x n_vertices.  And:
     * The i-th column of M is unit_vertex[i] and the last row all
     * 1's. The i-th column of Y is real_vertex[i].  If we split X=[A|b],
     * the least square approx is A x_hat+b Classically X = Y * (M^t (M
     * M^t)^{-1}) Let K = M^t * (M M^t)^{-1} = [KA Kb] this can be
     * precomputed, and that is exactly what we do.  Finally A = Y*KA and
     * b = Y*Kb.
     */
    template <int dim>
    struct TransformR2UAffine
    {
      static const double KA[GeometryInfo<dim>::vertices_per_cell][dim];
      static const double Kb[GeometryInfo<dim>::vertices_per_cell];
    };


    /*
      Octave code:
      M=[0 1; 1 1];
      K1 = transpose(M) * inverse (M*transpose(M));
      printf ("{%f, %f},\n", K1' );
    */
    template <>
    const double TransformR2UAffine<1>::KA[GeometryInfo<1>::vertices_per_cell]
                                          [1] = {{-1.000000}, {1.000000}};

    template <>
    const double TransformR2UAffine<1>::Kb[GeometryInfo<1>::vertices_per_cell] =
      {1.000000, 0.000000};


    /*
      Octave code:
      M=[0 1 0 1;0 0 1 1;1 1 1 1];
      K2 = transpose(M) * inverse (M*transpose(M));
      printf ("{%f, %f, %f},\n", K2' );
    */
    template <>
    const double TransformR2UAffine<2>::KA[GeometryInfo<2>::vertices_per_cell]
                                          [2] = {{-0.500000, -0.500000},
                                                 {0.500000, -0.500000},
                                                 {-0.500000, 0.500000},
                                                 {0.500000, 0.500000}};

    /*
      Octave code:
      M=[0 1 0 1 0 1 0 1;0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1; 1 1 1 1 1 1 1 1];
      K3 = transpose(M) * inverse (M*transpose(M))
      printf ("{%f, %f, %f, %f},\n", K3' );
    */
    template <>
    const double TransformR2UAffine<2>::Kb[GeometryInfo<2>::vertices_per_cell] =
      {0.750000, 0.250000, 0.250000, -0.250000};


    template <>
    const double TransformR2UAffine<3>::KA[GeometryInfo<3>::vertices_per_cell]
                                          [3] = {
                                            {-0.250000, -0.250000, -0.250000},
                                            {0.250000, -0.250000, -0.250000},
                                            {-0.250000, 0.250000, -0.250000},
                                            {0.250000, 0.250000, -0.250000},
                                            {-0.250000, -0.250000, 0.250000},
                                            {0.250000, -0.250000, 0.250000},
                                            {-0.250000, 0.250000, 0.250000},
                                            {0.250000, 0.250000, 0.250000}

    };


    template <>
    const double TransformR2UAffine<3>::Kb[GeometryInfo<3>::vertices_per_cell] =
      {0.500000,
       0.250000,
       0.250000,
       0.000000,
       0.250000,
       0.000000,
       0.000000,
       -0.250000};
  } // namespace



  template <int dim, int spacedim>
  std::pair<DerivativeForm<1, dim, spacedim>, Tensor<1, spacedim>>
  affine_cell_approximation(const ArrayView<const Point<spacedim>> &vertices)
  {
    AssertDimension(vertices.size(), GeometryInfo<dim>::vertices_per_cell);

    // A = vertex * KA
    DerivativeForm<1, dim, spacedim> A;

    for (unsigned int d = 0; d < spacedim; ++d)
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        for (unsigned int e = 0; e < dim; ++e)
          A[d][e] += vertices[v][d] * TransformR2UAffine<dim>::KA[v][e];

    // b = vertex * Kb
    Tensor<1, spacedim> b;
    for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
      b += vertices[v] * TransformR2UAffine<dim>::Kb[v];

    return std::make_pair(A, b);
  }



  template <int dim>
  Vector<double>
  compute_aspect_ratio_of_cells(const Mapping<dim>       &mapping,
                                const Triangulation<dim> &triangulation,
                                const Quadrature<dim>    &quadrature)
  {
    FE_Nothing<dim> fe;
    FEValues<dim>   fe_values(mapping, fe, quadrature, update_jacobians);

    Vector<double> aspect_ratio_vector(triangulation.n_active_cells());

    // loop over cells of processor
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            double aspect_ratio_cell = 0.0;

            fe_values.reinit(cell);

            // loop over quadrature points
            for (unsigned int q = 0; q < quadrature.size(); ++q)
              {
                const Tensor<2, dim, double> jacobian =
                  Tensor<2, dim, double>(fe_values.jacobian(q));

                // We intentionally do not want to throw an exception in case of
                // inverted elements since this is not the task of this
                // function. Instead, inf is written into the vector in case of
                // inverted elements.
                if (determinant(jacobian) <= 0)
                  {
                    aspect_ratio_cell = std::numeric_limits<double>::infinity();
                  }
                else
                  {
                    LAPACKFullMatrix<double> J = LAPACKFullMatrix<double>(dim);
                    for (unsigned int i = 0; i < dim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        J(i, j) = jacobian[i][j];

                    J.compute_svd();

                    const double max_sv = J.singular_value(0);
                    const double min_sv = J.singular_value(dim - 1);
                    const double ar     = max_sv / min_sv;

                    // Take the max between the previous and the current
                    // aspect ratio value; if we had previously encountered
                    // an inverted cell, we will have placed an infinity
                    // in the aspect_ratio_cell variable, and that value
                    // will survive this max operation.
                    aspect_ratio_cell = std::max(aspect_ratio_cell, ar);
                  }
              }

            // fill vector
            aspect_ratio_vector(cell->active_cell_index()) = aspect_ratio_cell;
          }
      }

    return aspect_ratio_vector;
  }



  template <int dim>
  double
  compute_maximum_aspect_ratio(const Mapping<dim>       &mapping,
                               const Triangulation<dim> &triangulation,
                               const Quadrature<dim>    &quadrature)
  {
    Vector<double> aspect_ratio_vector =
      compute_aspect_ratio_of_cells(mapping, triangulation, quadrature);

    return VectorTools::compute_global_error(triangulation,
                                             aspect_ratio_vector,
                                             VectorTools::Linfty_norm);
  }



  template <int dim, int spacedim>
  BoundingBox<spacedim>
  compute_bounding_box(const Triangulation<dim, spacedim> &tria)
  {
    using iterator =
      typename Triangulation<dim, spacedim>::active_cell_iterator;
    const auto predicate = [](const iterator &) { return true; };

    return compute_bounding_box(
      tria, std::function<bool(const iterator &)>(predicate));
  }



  template <int dim, int spacedim>
  double
  minimal_cell_diameter(const Triangulation<dim, spacedim> &triangulation,
                        const Mapping<dim, spacedim>       &mapping)
  {
    double min_diameter = std::numeric_limits<double>::max();
    for (const auto &cell : triangulation.active_cell_iterators())
      if (!cell->is_artificial())
        min_diameter = std::min(min_diameter, cell->diameter(mapping));

    const double global_min_diameter =
      Utilities::MPI::min(min_diameter, triangulation.get_mpi_communicator());
    return global_min_diameter;
  }



  template <int dim, int spacedim>
  double
  maximal_cell_diameter(const Triangulation<dim, spacedim> &triangulation,
                        const Mapping<dim, spacedim>       &mapping)
  {
    double max_diameter = 0.;
    for (const auto &cell : triangulation.active_cell_iterators())
      if (!cell->is_artificial())
        max_diameter = std::max(max_diameter, cell->diameter(mapping));

    const double global_max_diameter =
      Utilities::MPI::max(max_diameter, triangulation.get_mpi_communicator());
    return global_max_diameter;
  }
} /* namespace GridTools */


// explicit instantiations
#include "grid/grid_tools_geometry.inst"

DEAL_II_NAMESPACE_CLOSE

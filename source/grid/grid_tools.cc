// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/constrained_linear_operator.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <array>
#include <cmath>
#include <iostream>
#include <list>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_map>

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
#if defined(DEAL_II_WITH_P4EST) && defined(DEBUG)
    if (const parallel::distributed::Triangulation<dim, spacedim> *p_tria =
          dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> *>(&tria))
      Assert(p_tria->n_global_active_cells() == tria.n_cells(0),
             ExcNotImplemented());
#endif

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
  volume(const Triangulation<dim, spacedim> &triangulation,
         const Mapping<dim, spacedim> &      mapping)
  {
    // get the degree of the mapping if possible. if not, just assume 1
    unsigned int mapping_degree = 1;
    if (const auto *p =
          dynamic_cast<const MappingQGeneric<dim, spacedim> *>(&mapping))
      mapping_degree = p->get_degree();
    else if (const auto *p =
               dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping))
      mapping_degree = p->get_degree();

    // then initialize an appropriate quadrature formula
    const QGauss<dim>  quadrature_formula(mapping_degree + 1);
    const unsigned int n_q_points = quadrature_formula.size();

    // we really want the JxW values from the FEValues object, but it
    // wants a finite element. create a cheap element as a dummy
    // element
    FE_Nothing<dim, spacedim> dummy_fe;
    FEValues<dim, spacedim>   fe_values(mapping,
                                      dummy_fe,
                                      quadrature_formula,
                                      update_JxW_values);

    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();

    double local_volume = 0;

    // compute the integral quantities by quadrature
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          for (unsigned int q = 0; q < n_q_points; ++q)
            local_volume += fe_values.JxW(q);
        }

    double global_volume = 0;

#ifdef DEAL_II_WITH_MPI
    if (const parallel::TriangulationBase<dim, spacedim> *p_tria =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &triangulation))
      global_volume =
        Utilities::MPI::sum(local_volume, p_tria->get_communicator());
    else
#endif
      global_volume = local_volume;

    return global_volume;
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
  compute_aspect_ratio_of_cells(const Mapping<dim> &      mapping,
                                const Triangulation<dim> &triangulation,
                                const Quadrature<dim> &   quadrature)
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
                    for (unsigned int i = 0; i < dim; i++)
                      for (unsigned int j = 0; j < dim; j++)
                        J(i, j) = jacobian[i][j];

                    J.compute_svd();

                    double const max_sv = J.singular_value(0);
                    double const min_sv = J.singular_value(dim - 1);
                    double const ar     = max_sv / min_sv;

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
  compute_maximum_aspect_ratio(const Mapping<dim> &      mapping,
                               const Triangulation<dim> &triangulation,
                               const Quadrature<dim> &   quadrature)
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



  // Generic functions for appending face data in 2D or 3D. TODO: we can
  // remove these once we have 'if constexpr'.
  namespace internal
  {
    inline void
    append_face_data(const CellData<1> &face_data, SubCellData &subcell_data)
    {
      subcell_data.boundary_lines.push_back(face_data);
    }



    inline void
    append_face_data(const CellData<2> &face_data, SubCellData &subcell_data)
    {
      subcell_data.boundary_quads.push_back(face_data);
    }



    // Lexical comparison for sorting CellData objects.
    template <int structdim>
    struct CellDataComparator
    {
      bool
      operator()(const CellData<structdim> &a,
                 const CellData<structdim> &b) const
      {
        // Check vertices:
        if (std::lexicographical_compare(std::begin(a.vertices),
                                         std::end(a.vertices),
                                         std::begin(b.vertices),
                                         std::end(b.vertices)))
          return true;
          // it should never be necessary to check the material or manifold
          // ids as a 'tiebreaker' (since they must be equal if the vertex
          // indices are equal). Assert it anyway:
#ifdef DEBUG
        if (std::equal(std::begin(a.vertices),
                       std::end(a.vertices),
                       std::begin(b.vertices)))
          {
            Assert(a.material_id == b.material_id &&
                     a.manifold_id == b.manifold_id,
                   ExcMessage(
                     "Two CellData objects with equal vertices must "
                     "have the same material/boundary ids and manifold "
                     "ids."));
          }
#endif
        return false;
      }
    };


    /**
     * get_coarse_mesh_description() needs to store face data for dim>1, but
     * we can not have this code in the function, as this requires either an
     * instantiation of CellData<dim-1>, or constexpr if. We use a class with
     * specialization instead for now.
     *
     * Data on faces is added with insert_face_data() and then retrieved with
     * get().
     */
    template <int dim>
    class FaceDataHelper
    {
    public:
      /**
       * Store the data about the given face @p face.
       */
      template <class FaceIteratorType>
      void
      insert_face_data(const FaceIteratorType &face)
      {
        CellData<dim - 1> face_cell_data;
        for (unsigned int vertex_n = 0; vertex_n < face->n_vertices();
             ++vertex_n)
          face_cell_data.vertices[vertex_n] = face->vertex_index(vertex_n);
        face_cell_data.boundary_id = face->boundary_id();
        face_cell_data.manifold_id = face->manifold_id();

        face_data.insert(face_cell_data);
      }

      /**
       * Return  the @p subcell_data with the stored information.
       */
      SubCellData
      get()
      {
        SubCellData subcell_data;

        for (const CellData<dim - 1> &face_cell_data : face_data)
          internal::append_face_data(face_cell_data, subcell_data);
        return subcell_data;
      }


    private:
      std::set<CellData<dim - 1>, internal::CellDataComparator<dim - 1>>
        face_data;
    };


    // Do nothing for dim=1:
    template <>
    class FaceDataHelper<1>
    {
    public:
      template <class FaceIteratorType>
      void
      insert_face_data(const FaceIteratorType &)
      {}

      SubCellData
      get()
      {
        return SubCellData();
      }
    };
  } // namespace internal



  template <int dim, int spacedim>
  std::
    tuple<std::vector<Point<spacedim>>, std::vector<CellData<dim>>, SubCellData>
    get_coarse_mesh_description(const Triangulation<dim, spacedim> &tria)
  {
    Assert(1 <= tria.n_levels(),
           ExcMessage("The input triangulation must be non-empty."));

    std::vector<Point<spacedim>> vertices;
    std::vector<CellData<dim>>   cells;

    unsigned int max_level_0_vertex_n = 0;
    for (const auto &cell : tria.cell_iterators_on_level(0))
      for (const unsigned int cell_vertex_n : cell->vertex_indices())
        max_level_0_vertex_n =
          std::max(cell->vertex_index(cell_vertex_n), max_level_0_vertex_n);
    vertices.resize(max_level_0_vertex_n + 1);

    internal::FaceDataHelper<dim> face_data;
    std::set<CellData<1>, internal::CellDataComparator<1>>
      line_data; // only used in 3D

    for (const auto &cell : tria.cell_iterators_on_level(0))
      {
        // Save cell data
        CellData<dim> cell_data(cell->n_vertices());
        for (const unsigned int cell_vertex_n : cell->vertex_indices())
          {
            Assert(cell->vertex_index(cell_vertex_n) < vertices.size(),
                   ExcInternalError());
            vertices[cell->vertex_index(cell_vertex_n)] =
              cell->vertex(cell_vertex_n);
            cell_data.vertices[cell_vertex_n] =
              cell->vertex_index(cell_vertex_n);
          }
        cell_data.material_id = cell->material_id();
        cell_data.manifold_id = cell->manifold_id();
        cells.push_back(cell_data);

        // Save face data
        if (dim > 1)
          {
            for (const unsigned int face_n : cell->face_indices())
              face_data.insert_face_data(cell->face(face_n));
          }
        // Save line data
        if (dim == 3)
          {
            for (unsigned int line_n = 0; line_n < cell->n_lines(); ++line_n)
              {
                const auto  line = cell->line(line_n);
                CellData<1> line_cell_data;
                for (unsigned int vertex_n = 0; vertex_n < line->n_vertices();
                     ++vertex_n)
                  line_cell_data.vertices[vertex_n] =
                    line->vertex_index(vertex_n);
                line_cell_data.boundary_id = line->boundary_id();
                line_cell_data.manifold_id = line->manifold_id();

                line_data.insert(line_cell_data);
              }
          }
      }

      // Double-check that there are no unused vertices:
#ifdef DEBUG
    {
      std::vector<bool> used_vertices(vertices.size());
      for (const CellData<dim> &cell_data : cells)
        for (const auto v : cell_data.vertices)
          used_vertices[v] = true;
      Assert(std::find(used_vertices.begin(), used_vertices.end(), false) ==
               used_vertices.end(),
             ExcMessage("The level zero vertices should form a contiguous "
                        "range."));
    }
#endif

    SubCellData subcell_data = face_data.get();

    if (dim == 3)
      for (const CellData<1> &face_line_data : line_data)
        subcell_data.boundary_lines.push_back(face_line_data);

    return std::tuple<std::vector<Point<spacedim>>,
                      std::vector<CellData<dim>>,
                      SubCellData>(std::move(vertices),
                                   std::move(cells),
                                   std::move(subcell_data));
  }



  template <int dim, int spacedim>
  void
  delete_unused_vertices(std::vector<Point<spacedim>> &vertices,
                         std::vector<CellData<dim>> &  cells,
                         SubCellData &                 subcelldata)
  {
    Assert(
      subcelldata.check_consistency(dim),
      ExcMessage(
        "Invalid SubCellData supplied according to ::check_consistency(). "
        "This is caused by data containing objects for the wrong dimension."));

    // first check which vertices are actually used
    std::vector<bool> vertex_used(vertices.size(), false);
    for (unsigned int c = 0; c < cells.size(); ++c)
      for (unsigned int v = 0; v < cells[c].vertices.size(); ++v)
        {
          Assert(cells[c].vertices[v] < vertices.size(),
                 ExcMessage("Invalid vertex index encountered! cells[" +
                            Utilities::int_to_string(c) + "].vertices[" +
                            Utilities::int_to_string(v) + "]=" +
                            Utilities::int_to_string(cells[c].vertices[v]) +
                            " is invalid, because only " +
                            Utilities::int_to_string(vertices.size()) +
                            " vertices were supplied."));
          vertex_used[cells[c].vertices[v]] = true;
        }


    // then renumber the vertices that are actually used in the same order as
    // they were beforehand
    const unsigned int        invalid_vertex = numbers::invalid_unsigned_int;
    std::vector<unsigned int> new_vertex_numbers(vertices.size(),
                                                 invalid_vertex);
    unsigned int              next_free_number = 0;
    for (unsigned int i = 0; i < vertices.size(); ++i)
      if (vertex_used[i] == true)
        {
          new_vertex_numbers[i] = next_free_number;
          ++next_free_number;
        }

    // next replace old vertex numbers by the new ones
    for (unsigned int c = 0; c < cells.size(); ++c)
      for (auto &v : cells[c].vertices)
        v = new_vertex_numbers[v];

    // same for boundary data
    for (unsigned int c = 0; c < subcelldata.boundary_lines.size(); // NOLINT
         ++c)
      for (unsigned int v = 0;
           v < subcelldata.boundary_lines[c].vertices.size();
           ++v)
        {
          Assert(subcelldata.boundary_lines[c].vertices[v] <
                   new_vertex_numbers.size(),
                 ExcMessage(
                   "Invalid vertex index in subcelldata.boundary_lines. "
                   "subcelldata.boundary_lines[" +
                   Utilities::int_to_string(c) + "].vertices[" +
                   Utilities::int_to_string(v) + "]=" +
                   Utilities::int_to_string(
                     subcelldata.boundary_lines[c].vertices[v]) +
                   " is invalid, because only " +
                   Utilities::int_to_string(vertices.size()) +
                   " vertices were supplied."));
          subcelldata.boundary_lines[c].vertices[v] =
            new_vertex_numbers[subcelldata.boundary_lines[c].vertices[v]];
        }

    for (unsigned int c = 0; c < subcelldata.boundary_quads.size(); // NOLINT
         ++c)
      for (unsigned int v = 0;
           v < subcelldata.boundary_quads[c].vertices.size();
           ++v)
        {
          Assert(subcelldata.boundary_quads[c].vertices[v] <
                   new_vertex_numbers.size(),
                 ExcMessage(
                   "Invalid vertex index in subcelldata.boundary_quads. "
                   "subcelldata.boundary_quads[" +
                   Utilities::int_to_string(c) + "].vertices[" +
                   Utilities::int_to_string(v) + "]=" +
                   Utilities::int_to_string(
                     subcelldata.boundary_quads[c].vertices[v]) +
                   " is invalid, because only " +
                   Utilities::int_to_string(vertices.size()) +
                   " vertices were supplied."));

          subcelldata.boundary_quads[c].vertices[v] =
            new_vertex_numbers[subcelldata.boundary_quads[c].vertices[v]];
        }

    // finally copy over the vertices which we really need to a new array and
    // replace the old one by the new one
    std::vector<Point<spacedim>> tmp;
    tmp.reserve(std::count(vertex_used.begin(), vertex_used.end(), true));
    for (unsigned int v = 0; v < vertices.size(); ++v)
      if (vertex_used[v] == true)
        tmp.push_back(vertices[v]);
    swap(vertices, tmp);
  }



  template <int dim, int spacedim>
  void
  delete_duplicated_vertices(std::vector<Point<spacedim>> &vertices,
                             std::vector<CellData<dim>> &  cells,
                             SubCellData &                 subcelldata,
                             std::vector<unsigned int> &   considered_vertices,
                             const double                  tol)
  {
    AssertIndexRange(2, vertices.size());
    // create a vector of vertex indices. initialize it to the identity, later
    // on change that if necessary.
    std::vector<unsigned int> new_vertex_numbers(vertices.size());
    std::iota(new_vertex_numbers.begin(), new_vertex_numbers.end(), 0);

    // if the considered_vertices vector is empty, consider all vertices
    if (considered_vertices.size() == 0)
      considered_vertices = new_vertex_numbers;
    Assert(considered_vertices.size() <= vertices.size(), ExcInternalError());

    // The algorithm below improves upon the naive O(n^2) algorithm by first
    // sorting vertices by their value in one component and then only
    // comparing vertices for equality which are nearly equal in that
    // component. For example, if @p vertices form a cube, then we will only
    // compare points that have the same x coordinate when we try to find
    // duplicated vertices.

    // Start by finding the longest coordinate direction. This minimizes the
    // number of points that need to be compared against each-other in a
    // single set for typical geometries.
    const BoundingBox<spacedim> bbox(vertices);
    const auto &                min = bbox.get_boundary_points().first;
    const auto &                max = bbox.get_boundary_points().second;

    unsigned int longest_coordinate_direction = 0;
    double       longest_coordinate_length    = max[0] - min[0];
    for (unsigned int d = 1; d < spacedim; ++d)
      {
        const double coordinate_length = max[d] - min[d];
        if (longest_coordinate_length < coordinate_length)
          {
            longest_coordinate_length    = coordinate_length;
            longest_coordinate_direction = d;
          }
      }

    // Sort vertices (while preserving their vertex numbers) along that
    // coordinate direction:
    std::vector<std::pair<unsigned int, Point<spacedim>>> sorted_vertices;
    sorted_vertices.reserve(vertices.size());
    for (const unsigned int vertex_n : considered_vertices)
      {
        AssertIndexRange(vertex_n, vertices.size());
        sorted_vertices.emplace_back(vertex_n, vertices[vertex_n]);
      }
    std::sort(sorted_vertices.begin(),
              sorted_vertices.end(),
              [&](const std::pair<unsigned int, Point<spacedim>> &a,
                  const std::pair<unsigned int, Point<spacedim>> &b) {
                return a.second[longest_coordinate_direction] <
                       b.second[longest_coordinate_direction];
              });

    auto within_tolerance = [=](const Point<spacedim> &a,
                                const Point<spacedim> &b) {
      for (unsigned int d = 0; d < spacedim; ++d)
        if (std::abs(a[d] - b[d]) > tol)
          return false;
      return true;
    };

    // Find a range of numbers that have the same component in the longest
    // coordinate direction:
    auto range_start = sorted_vertices.begin();
    while (range_start != sorted_vertices.end())
      {
        auto range_end = range_start + 1;
        while (range_end != sorted_vertices.end() &&
               std::abs(range_end->second[longest_coordinate_direction] -
                        range_start->second[longest_coordinate_direction]) <
                 tol)
          ++range_end;

        // preserve behavior with older versions of this function by replacing
        // higher vertex numbers by lower vertex numbers
        std::sort(range_start,
                  range_end,
                  [](const std::pair<unsigned int, Point<spacedim>> &a,
                     const std::pair<unsigned int, Point<spacedim>> &b) {
                    return a.first < b.first;
                  });

        // Now de-duplicate [range_start, range_end)
        //
        // We have identified all points that are within a strip of width 'tol'
        // in one coordinate direction. Now we need to figure out which of these
        // are also close in other coordinate directions. If two are close, we
        // can mark the second one for deletion.
        for (auto reference = range_start; reference != range_end; ++reference)
          {
            if (reference->first != numbers::invalid_unsigned_int)
              for (auto it = reference + 1; it != range_end; ++it)
                {
                  if (within_tolerance(reference->second, it->second))
                    {
                      new_vertex_numbers[it->first] = reference->first;
                      // skip the replaced vertex in the future
                      it->first = numbers::invalid_unsigned_int;
                    }
                }
          }
        range_start = range_end;
      }

    // now we got a renumbering list. simply renumber all vertices
    // (non-duplicate vertices get renumbered to themselves, so nothing bad
    // happens). after that, the duplicate vertices will be unused, so call
    // delete_unused_vertices() to do that part of the job.
    for (auto &cell : cells)
      for (auto &vertex_index : cell.vertices)
        vertex_index = new_vertex_numbers[vertex_index];
    for (auto &quad : subcelldata.boundary_quads)
      for (auto &vertex_index : quad.vertices)
        vertex_index = new_vertex_numbers[vertex_index];
    for (auto &line : subcelldata.boundary_lines)
      for (auto &vertex_index : line.vertices)
        vertex_index = new_vertex_numbers[vertex_index];

    delete_unused_vertices(vertices, cells, subcelldata);
  }



  template <int dim, int spacedim>
  void
  invert_all_negative_measure_cells(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>> &        cells)
  {
    if (dim == 1)
      return;
    if (dim == 2 && spacedim == 3)
      Assert(false, ExcNotImplemented());

    std::size_t n_negative_cells = 0;
    for (auto &cell : cells)
      {
        Assert(cell.vertices.size() ==
                 ReferenceCells::get_hypercube<dim>().n_vertices(),
               ExcNotImplemented());
        const ArrayView<const unsigned int> vertices(cell.vertices);
        if (GridTools::cell_measure<dim>(all_vertices, vertices) < 0)
          {
            ++n_negative_cells;

            // TODO: this only works for quads and hexes
            if (dim == 2)
              {
                // flip the cell across the y = x line in 2D
                std::swap(cell.vertices[1], cell.vertices[2]);
              }
            else if (dim == 3)
              {
                // swap the front and back faces in 3D
                std::swap(cell.vertices[0], cell.vertices[2]);
                std::swap(cell.vertices[1], cell.vertices[3]);
                std::swap(cell.vertices[4], cell.vertices[6]);
                std::swap(cell.vertices[5], cell.vertices[7]);
              }

            // Check whether the resulting cell is now ok.
            // If not, then the grid is seriously broken and
            // we just give up.
            AssertThrow(GridTools::cell_measure<dim>(all_vertices, vertices) >
                          0,
                        ExcInternalError());
          }
      }

    // We assume that all cells of a grid have
    // either positive or negative volumes but
    // not both mixed. Although above reordering
    // might work also on single cells, grids
    // with both kind of cells are very likely to
    // be broken. Check for this here.
    AssertThrow(n_negative_cells == 0 || n_negative_cells == cells.size(),
                ExcMessage(
                  std::string(
                    "This function assumes that either all cells have positive "
                    "volume, or that all cells have been specified in an "
                    "inverted vertex order so that their volume is negative. "
                    "(In the latter case, this class automatically inverts "
                    "every cell.) However, the mesh you have specified "
                    "appears to have both cells with positive and cells with "
                    "negative volume. You need to check your mesh which "
                    "cells these are and how they got there.\n"
                    "As a hint, of the total ") +
                  std::to_string(cells.size()) + " cells in the mesh, " +
                  std::to_string(n_negative_cells) +
                  " appear to have a negative volume."));
  }



  // Functions and classes for consistently_order_cells
  namespace
  {
    /**
     * A simple data structure denoting an edge, i.e., the ordered pair
     * of its vertex indices. This is only used in the is_consistent()
     * function.
     */
    struct CheapEdge
    {
      /**
       * Construct an edge from the global indices of its two vertices.
       */
      CheapEdge(const unsigned int v0, const unsigned int v1)
        : v0(v0)
        , v1(v1)
      {}

      /**
       * Comparison operator for edges. It compares based on the
       * lexicographic ordering of the two vertex indices.
       */
      bool
      operator<(const CheapEdge &e) const
      {
        return ((v0 < e.v0) || ((v0 == e.v0) && (v1 < e.v1)));
      }

    private:
      /**
       * The global indices of the vertices that define the edge.
       */
      const unsigned int v0, v1;
    };


    /**
     * A function that determines whether the edges in a mesh are
     * already consistently oriented. It does so by adding all edges
     * of all cells into a set (which automatically eliminates
     * duplicates) but before that checks whether the reverse edge is
     * already in the set -- which would imply that a neighboring cell
     * is inconsistently oriented.
     */
    template <int dim>
    bool
    is_consistent(const std::vector<CellData<dim>> &cells)
    {
      std::set<CheapEdge> edges;

      for (typename std::vector<CellData<dim>>::const_iterator c =
             cells.begin();
           c != cells.end();
           ++c)
        {
          // construct the edges in reverse order. for each of them,
          // ensure that the reverse edge is not yet in the list of
          // edges (return false if the reverse edge already *is* in
          // the list) and then add the actual edge to it; std::set
          // eliminates duplicates automatically
          for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
            {
              const CheapEdge reverse_edge(
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 1)],
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 0)]);
              if (edges.find(reverse_edge) != edges.end())
                return false;


              // ok, not. insert edge in correct order
              const CheapEdge correct_edge(
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 0)],
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 1)]);
              edges.insert(correct_edge);
            }
        }

      // no conflicts found, so return true
      return true;
    }


    /**
     * A structure that describes some properties of parallel edges
     * such as what starter edges are (i.e., representative elements
     * of the sets of parallel edges within a cell) and what the set
     * of parallel edges to each edge is.
     */
    template <int dim>
    struct ParallelEdges
    {
      /**
       * An array that contains the indices of dim edges that can
       * serve as (arbitrarily chosen) starting points for the
       * dim sets of parallel edges within each cell.
       */
      static const unsigned int starter_edges[dim];

      /**
       * Number and indices of all of those edges parallel to each of the
       * edges in a cell.
       */
      static const unsigned int n_other_parallel_edges = (1 << (dim - 1)) - 1;
      static const unsigned int
        parallel_edges[GeometryInfo<dim>::lines_per_cell]
                      [n_other_parallel_edges];
    };

    template <>
    const unsigned int ParallelEdges<2>::starter_edges[2] = {0, 2};

    template <>
    const unsigned int ParallelEdges<2>::parallel_edges[4][1] = {{1},
                                                                 {0},
                                                                 {3},
                                                                 {2}};

    template <>
    const unsigned int ParallelEdges<3>::starter_edges[3] = {0, 2, 8};

    template <>
    const unsigned int ParallelEdges<3>::parallel_edges[12][3] = {
      {1, 4, 5},   // line 0
      {0, 4, 5},   // line 1
      {3, 6, 7},   // line 2
      {2, 6, 7},   // line 3
      {0, 1, 5},   // line 4
      {0, 1, 4},   // line 5
      {2, 3, 7},   // line 6
      {2, 3, 6},   // line 7
      {9, 10, 11}, // line 8
      {8, 10, 11}, // line 9
      {8, 9, 11},  // line 10
      {8, 9, 10}   // line 11
    };


    /**
     * A structure that store the index of a cell and, crucially, how a
     * given edge relates to this cell.
     */
    struct AdjacentCell
    {
      /**
       * Default constructor. Initialize the fields with invalid values.
       */
      AdjacentCell()
        : cell_index(numbers::invalid_unsigned_int)
        , edge_within_cell(numbers::invalid_unsigned_int)
      {}

      /**
       * Constructor. Initialize the fields with the given values.
       */
      AdjacentCell(const unsigned int cell_index,
                   const unsigned int edge_within_cell)
        : cell_index(cell_index)
        , edge_within_cell(edge_within_cell)
      {}


      unsigned int cell_index;
      unsigned int edge_within_cell;
    };



    template <int dim>
    class AdjacentCells;

    /**
     * A class that represents all of the cells adjacent to a given edge.
     * This class corresponds to the 2d case where each edge has at most
     * two adjacent cells.
     */
    template <>
    class AdjacentCells<2>
    {
    public:
      /**
       * An iterator that allows iterating over all cells adjacent
       * to the edge represented by the current object.
       */
      using const_iterator = const AdjacentCell *;

      /**
       * Add the given cell to the collection of cells adjacent to
       * the edge this object corresponds to. Since we are covering
       * the 2d case, the set of adjacent cells currently
       * represented by this object must have either zero or
       * one element already, since we can not add more than two
       * adjacent cells for each edge.
       */
      void
      push_back(const AdjacentCell &adjacent_cell)
      {
        if (adjacent_cells[0].cell_index == numbers::invalid_unsigned_int)
          adjacent_cells[0] = adjacent_cell;
        else
          {
            Assert(adjacent_cells[1].cell_index ==
                     numbers::invalid_unsigned_int,
                   ExcInternalError());
            adjacent_cells[1] = adjacent_cell;
          }
      }


      /**
       * Return an iterator to the first valid cell stored as adjacent to the
       * edge represented by the current object.
       */
      const_iterator
      begin() const
      {
        return adjacent_cells;
      }


      /**
       * Return an iterator to the element past the last valid cell stored
       * as adjacent to the edge represented by the current object.
       * @return
       */
      const_iterator
      end() const
      {
        // check whether the current object stores zero, one, or two
        // adjacent cells, and use this to point to the element past the
        // last valid one
        if (adjacent_cells[0].cell_index == numbers::invalid_unsigned_int)
          return adjacent_cells;
        else if (adjacent_cells[1].cell_index == numbers::invalid_unsigned_int)
          return adjacent_cells + 1;
        else
          return adjacent_cells + 2;
      }

    private:
      /**
       * References to the (at most) two cells that are adjacent to
       * the edge this object corresponds to. Unused elements are
       * default-initialized and have invalid values; in particular,
       * their cell_index field equals numbers::invalid_unsigned_int.
       */
      AdjacentCell adjacent_cells[2];
    };



    /**
     * A class that represents all of the cells adjacent to a given edge.
     * This class corresponds to the 3d case where each edge can have an
     * arbitrary number of adjacent cells. We represent this as a
     * std::vector<AdjacentCell>, from which class the current one is
     * derived and from which it inherits all of its member functions.
     */
    template <>
    class AdjacentCells<3> : public std::vector<AdjacentCell>
    {};


    /**
     * A class that describes all of the relevant properties of an
     * edge. For the purpose of what we do here, that includes the
     * indices of the two vertices, and the indices of the adjacent
     * cells (together with a description *where* in each of the
     * adjacent cells the edge is located). It also includes the
     * (global) direction of the edge: either from the first vertex to
     * the second, the other way around, or so far undetermined.
     */
    template <int dim>
    class Edge
    {
    public:
      /**
       * Constructor. Create the edge based on the information given
       * in @p cell, and selecting the edge with number @p edge_number
       * within this cell. Initialize the edge as unoriented.
       */
      Edge(const CellData<dim> &cell, const unsigned int edge_number)
        : orientation_status(not_oriented)
      {
        Assert(edge_number < GeometryInfo<dim>::lines_per_cell,
               ExcInternalError());

        // copy vertices for this particular line
        vertex_indices[0] =
          cell
            .vertices[GeometryInfo<dim>::line_to_cell_vertices(edge_number, 0)];
        vertex_indices[1] =
          cell
            .vertices[GeometryInfo<dim>::line_to_cell_vertices(edge_number, 1)];

        // bring them into standard orientation
        if (vertex_indices[0] > vertex_indices[1])
          std::swap(vertex_indices[0], vertex_indices[1]);
      }

      /**
       * Comparison operator for edges. It compares based on the
       * lexicographic ordering of the two vertex indices.
       */
      bool
      operator<(const Edge<dim> &e) const
      {
        return ((vertex_indices[0] < e.vertex_indices[0]) ||
                ((vertex_indices[0] == e.vertex_indices[0]) &&
                 (vertex_indices[1] < e.vertex_indices[1])));
      }

      /**
       * Compare two edges for equality based on their vertex indices.
       */
      bool
      operator==(const Edge<dim> &e) const
      {
        return ((vertex_indices[0] == e.vertex_indices[0]) &&
                (vertex_indices[1] == e.vertex_indices[1]));
      }

      /**
       * The global indices of the two vertices that bound this edge. These
       * will be ordered so that the first index is less than the second.
       */
      unsigned int vertex_indices[2];

      /**
       * An enum that indicates the direction of this edge with
       * regard to the two vertices that bound it.
       */
      enum OrientationStatus
      {
        not_oriented,
        forward,
        backward
      };

      OrientationStatus orientation_status;

      /**
       * Store the set of cells adjacent to this edge (these cells then
       * also store *where* in the cell the edge is located).
       */
      AdjacentCells<dim> adjacent_cells;
    };



    /**
     * A data structure that represents a cell with all of its vertices
     * and edges.
     */
    template <int dim>
    struct Cell
    {
      /**
       * Construct a Cell object from a CellData object. Also take a
       * (sorted) list of edges and to point the edges of the current
       * object into this list of edges.
       */
      Cell(const CellData<dim> &c, const std::vector<Edge<dim>> &edge_list)
      {
        for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
          vertex_indices[i] = c.vertices[i];

        // now for each of the edges of this cell, find the location inside the
        // given edge_list array and store than index
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          {
            const Edge<dim> e(c, l);
            edge_indices[l] =
              (std::lower_bound(edge_list.begin(), edge_list.end(), e) -
               edge_list.begin());
            Assert(edge_indices[l] < edge_list.size(), ExcInternalError());
            Assert(edge_list[edge_indices[l]] == e, ExcInternalError())
          }
      }

      /**
       * A list of global indices for the vertices that bound this cell.
       */
      unsigned int vertex_indices[GeometryInfo<dim>::vertices_per_cell];

      /**
       * A list of indices into the 'edge_list' array passed to the constructor
       * for the edges of the current cell.
       */
      unsigned int edge_indices[GeometryInfo<dim>::lines_per_cell];
    };



    template <int dim>
    class EdgeDeltaSet;

    /**
     * A class that represents by how much the set of parallel edges
     * grows in each step. In the graph orientation paper, this set is
     * called $\Delta_k$, thus the name.
     *
     * In 2d, this set can only include zero, one, or two elements.
     * Consequently, the appropriate data structure is one in which
     * we store at most 2 elements in a fixed sized data structure.
     */
    template <>
    class EdgeDeltaSet<2>
    {
    public:
      /**
       * Iterator type for the elements of the set.
       */
      using const_iterator = const unsigned int *;

      /**
       * Default constructor. Initialize both slots as unused, corresponding
       * to an empty set.
       */
      EdgeDeltaSet()
      {
        edge_indices[0] = edge_indices[1] = numbers::invalid_unsigned_int;
      }


      /**
       * Delete the elements of the set by marking both slots as unused.
       */
      void
      clear()
      {
        edge_indices[0] = edge_indices[1] = numbers::invalid_unsigned_int;
      }

      /**
       * Insert one element into the set. This will fail if the set already
       * has two elements.
       */
      void
      insert(const unsigned int edge_index)
      {
        if (edge_indices[0] == numbers::invalid_unsigned_int)
          edge_indices[0] = edge_index;
        else
          {
            Assert(edge_indices[1] == numbers::invalid_unsigned_int,
                   ExcInternalError());
            edge_indices[1] = edge_index;
          }
      }


      /**
       * Return an iterator pointing to the first element of the set.
       */
      const_iterator
      begin() const
      {
        return edge_indices;
      }


      /**
       * Return an iterator pointing to the element past the last used one.
       */
      const_iterator
      end() const
      {
        // check whether the current object stores zero, one, or two
        // indices, and use this to point to the element past the
        // last valid one
        if (edge_indices[0] == numbers::invalid_unsigned_int)
          return edge_indices;
        else if (edge_indices[1] == numbers::invalid_unsigned_int)
          return edge_indices + 1;
        else
          return edge_indices + 2;
      }

    private:
      /**
       * Storage space to store the indices of at most two edges.
       */
      unsigned int edge_indices[2];
    };



    /**
     * A class that represents by how much the set of parallel edges
     * grows in each step. In the graph orientation paper, this set is
     * called $\Delta_k$, thus the name.
     *
     * In 3d, this set can have arbitrarily many elements, unlike the
     * 2d case specialized above. Consequently, we simply represent
     * the data structure with a std::set. Class derivation ensures
     * that we simply inherit all of the member functions of the
     * base class.
     */
    template <>
    class EdgeDeltaSet<3> : public std::set<unsigned int>
    {};



    /**
     * From a list of cells, build a sorted vector that contains all of the
     * edges that exist in the mesh.
     */
    template <int dim>
    std::vector<Edge<dim>>
    build_edges(const std::vector<CellData<dim>> &cells)
    {
      // build the edge list for all cells. because each cell has
      // GeometryInfo<dim>::lines_per_cell edges, the total number
      // of edges is this many times the number of cells. of course
      // some of them will be duplicates, and we throw them out below
      std::vector<Edge<dim>> edge_list;
      edge_list.reserve(cells.size() * GeometryInfo<dim>::lines_per_cell);
      for (unsigned int i = 0; i < cells.size(); ++i)
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          edge_list.emplace_back(cells[i], l);

      // next sort the edge list and then remove duplicates
      std::sort(edge_list.begin(), edge_list.end());
      edge_list.erase(std::unique(edge_list.begin(), edge_list.end()),
                      edge_list.end());

      return edge_list;
    }



    /**
     * Build the cell list. Update the edge array to let edges know
     * which cells are adjacent to them.
     */
    template <int dim>
    std::vector<Cell<dim>>
    build_cells_and_connect_edges(const std::vector<CellData<dim>> &cells,
                                  std::vector<Edge<dim>> &          edges)
    {
      std::vector<Cell<dim>> cell_list;
      cell_list.reserve(cells.size());
      for (unsigned int i = 0; i < cells.size(); ++i)
        {
          // create our own data structure for the cells and let it
          // connect to the edges array
          cell_list.emplace_back(cells[i], edges);

          // then also inform the edges that they are adjacent
          // to the current cell, and where within this cell
          for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
            edges[cell_list.back().edge_indices[l]].adjacent_cells.push_back(
              AdjacentCell(i, l));
        }
      Assert(cell_list.size() == cells.size(), ExcInternalError());

      return cell_list;
    }



    /**
     * Return the index within 'cells' of the first cell that has at least one
     * edge that is not yet oriented.
     */
    template <int dim>
    unsigned int
    get_next_unoriented_cell(const std::vector<Cell<dim>> &cells,
                             const std::vector<Edge<dim>> &edges,
                             const unsigned int            current_cell)
    {
      for (unsigned int c = current_cell; c < cells.size(); ++c)
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          if (edges[cells[c].edge_indices[l]].orientation_status ==
              Edge<dim>::not_oriented)
            return c;

      return numbers::invalid_unsigned_int;
    }



    /**
     * Given a set of cells and edges, orient all edges that are
     * (global) parallel to the one identified by the @p cell and
     * within it the one with index @p local_edge.
     */
    template <int dim>
    void
    orient_one_set_of_parallel_edges(const std::vector<Cell<dim>> &cells,
                                     std::vector<Edge<dim>> &      edges,
                                     const unsigned int            cell,
                                     const unsigned int            local_edge)
    {
      // choose the direction of the first edge. we have free choice
      // here and could simply choose "forward" if that's what pleases
      // us. however, for backward compatibility with the previous
      // implementation used till 2016, let us just choose the
      // direction so that it matches what we have in the given cell.
      //
      // in fact, in what can only be assumed to be a bug in the
      // original implementation, after orienting all edges, the code
      // that rotates the cells so that they match edge orientations
      // (see the rotate_cell() function below) rotated the cell two
      // more times by 90 degrees. this is ok -- it simply flips all
      // edge orientations, which leaves them valid. rather than do
      // the same in the current implementation, we can achieve the
      // same effect by modifying the rule above to choose the
      // direction of the starting edge of this parallel set
      // *opposite* to what it looks like in the current cell
      //
      // this bug only existed in the 2d implementation since there
      // were different implementations for 2d and 3d. consequently,
      // only replicate it for the 2d case and be "intuitive" in 3d.
      if (edges[cells[cell].edge_indices[local_edge]].vertex_indices[0] ==
          cells[cell].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(
            local_edge, 0)])
        // orient initial edge *opposite* to the way it is in the cell
        // (see above for the reason)
        edges[cells[cell].edge_indices[local_edge]].orientation_status =
          (dim == 2 ? Edge<dim>::backward : Edge<dim>::forward);
      else
        {
          Assert(
            edges[cells[cell].edge_indices[local_edge]].vertex_indices[0] ==
              cells[cell].vertex_indices
                [GeometryInfo<dim>::line_to_cell_vertices(local_edge, 1)],
            ExcInternalError());
          Assert(
            edges[cells[cell].edge_indices[local_edge]].vertex_indices[1] ==
              cells[cell].vertex_indices
                [GeometryInfo<dim>::line_to_cell_vertices(local_edge, 0)],
            ExcInternalError());

          // orient initial edge *opposite* to the way it is in the cell
          // (see above for the reason)
          edges[cells[cell].edge_indices[local_edge]].orientation_status =
            (dim == 2 ? Edge<dim>::forward : Edge<dim>::backward);
        }

      // walk outward from the given edge as described in
      // the algorithm in the paper that documents all of
      // this
      //
      // note that in 2d, each of the Deltas can at most
      // contain two elements, whereas in 3d it can be arbitrarily many
      EdgeDeltaSet<dim> Delta_k;
      EdgeDeltaSet<dim> Delta_k_minus_1;
      Delta_k_minus_1.insert(cells[cell].edge_indices[local_edge]);

      while (Delta_k_minus_1.begin() !=
             Delta_k_minus_1.end()) // while set is not empty
        {
          Delta_k.clear();

          for (typename EdgeDeltaSet<dim>::const_iterator delta =
                 Delta_k_minus_1.begin();
               delta != Delta_k_minus_1.end();
               ++delta)
            {
              Assert(edges[*delta].orientation_status !=
                       Edge<dim>::not_oriented,
                     ExcInternalError());

              // now go through the cells adjacent to this edge
              for (typename AdjacentCells<dim>::const_iterator adjacent_cell =
                     edges[*delta].adjacent_cells.begin();
                   adjacent_cell != edges[*delta].adjacent_cells.end();
                   ++adjacent_cell)
                {
                  const unsigned int K = adjacent_cell->cell_index;
                  const unsigned int delta_is_edge_in_K =
                    adjacent_cell->edge_within_cell;

                  // figure out the direction of delta with respect to the cell
                  // K (in the orientation in which the user has given it to us)
                  const unsigned int first_edge_vertex =
                    (edges[*delta].orientation_status == Edge<dim>::forward ?
                       edges[*delta].vertex_indices[0] :
                       edges[*delta].vertex_indices[1]);
                  const unsigned int first_edge_vertex_in_K =
                    cells[K]
                      .vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(
                        delta_is_edge_in_K, 0)];
                  Assert(
                    first_edge_vertex == first_edge_vertex_in_K ||
                      first_edge_vertex ==
                        cells[K].vertex_indices[GeometryInfo<
                          dim>::line_to_cell_vertices(delta_is_edge_in_K, 1)],
                    ExcInternalError());

                  // now figure out which direction the each of the "opposite"
                  // edges needs to be oriented into.
                  for (unsigned int o_e = 0;
                       o_e < ParallelEdges<dim>::n_other_parallel_edges;
                       ++o_e)
                    {
                      // get the index of the opposite edge and select which its
                      // first vertex needs to be based on how the current edge
                      // is oriented in the current cell
                      const unsigned int opposite_edge =
                        cells[K].edge_indices[ParallelEdges<
                          dim>::parallel_edges[delta_is_edge_in_K][o_e]];
                      const unsigned int first_opposite_edge_vertex =
                        cells[K].vertex_indices
                          [GeometryInfo<dim>::line_to_cell_vertices(
                            ParallelEdges<
                              dim>::parallel_edges[delta_is_edge_in_K][o_e],
                            (first_edge_vertex == first_edge_vertex_in_K ? 0 :
                                                                           1))];

                      // then determine the orientation of the edge based on
                      // whether the vertex we want to be the edge's first
                      // vertex is already the first vertex of the edge, or
                      // whether it points in the opposite direction
                      const typename Edge<dim>::OrientationStatus
                        opposite_edge_orientation =
                          (edges[opposite_edge].vertex_indices[0] ==
                               first_opposite_edge_vertex ?
                             Edge<dim>::forward :
                             Edge<dim>::backward);

                      // see if the opposite edge (there is only one in 2d) has
                      // already been oriented.
                      if (edges[opposite_edge].orientation_status ==
                          Edge<dim>::not_oriented)
                        {
                          // the opposite edge is not yet oriented. do orient it
                          // and add it to Delta_k
                          edges[opposite_edge].orientation_status =
                            opposite_edge_orientation;
                          Delta_k.insert(opposite_edge);
                        }
                      else
                        {
                          // this opposite edge has already been oriented. it
                          // should be consistent with the current one in 2d,
                          // while in 3d it may in fact be mis-oriented, and in
                          // that case the mesh will not be orientable. indicate
                          // this by throwing an exception that we can catch
                          // further up; this has the advantage that we can
                          // propagate through a couple of functions without
                          // having to do error checking and without modifying
                          // the 'cells' array that the user gave us
                          if (dim == 2)
                            {
                              Assert(edges[opposite_edge].orientation_status ==
                                       opposite_edge_orientation,
                                     ExcMeshNotOrientable());
                            }
                          else if (dim == 3)
                            {
                              if (edges[opposite_edge].orientation_status !=
                                  opposite_edge_orientation)
                                throw ExcMeshNotOrientable();
                            }
                          else
                            Assert(false, ExcNotImplemented());
                        }
                    }
                }
            }

          // finally copy the new set to the previous one
          // (corresponding to increasing 'k' by one in the
          // algorithm)
          Delta_k_minus_1 = Delta_k;
        }
    }


    /**
     * Given data structures @p cell_list and @p edge_list, where
     * all edges are already oriented, rotate the cell with
     * index @p cell_index in such a way that its local coordinate
     * system matches the ones of the adjacent edges. Store the
     * rotated order of vertices in <code>raw_cells[cell_index]</code>.
     */
    template <int dim>
    void
    rotate_cell(const std::vector<Cell<dim>> &cell_list,
                const std::vector<Edge<dim>> &edge_list,
                const unsigned int            cell_index,
                std::vector<CellData<dim>> &  raw_cells)
    {
      // find the first vertex of the cell. this is the vertex where dim edges
      // originate, so for each of the edges record which the starting vertex is
      unsigned int starting_vertex_of_edge[GeometryInfo<dim>::lines_per_cell];
      for (unsigned int e = 0; e < GeometryInfo<dim>::lines_per_cell; ++e)
        {
          Assert(edge_list[cell_list[cell_index].edge_indices[e]]
                     .orientation_status != Edge<dim>::not_oriented,
                 ExcInternalError());
          if (edge_list[cell_list[cell_index].edge_indices[e]]
                .orientation_status == Edge<dim>::forward)
            starting_vertex_of_edge[e] =
              edge_list[cell_list[cell_index].edge_indices[e]]
                .vertex_indices[0];
          else
            starting_vertex_of_edge[e] =
              edge_list[cell_list[cell_index].edge_indices[e]]
                .vertex_indices[1];
        }

      // find the vertex number that appears dim times. this will then be
      // the vertex at which we want to locate the origin of the cell's
      // coordinate system (i.e., vertex 0)
      unsigned int origin_vertex_of_cell = numbers::invalid_unsigned_int;
      switch (dim)
        {
          case 2:
            {
              // in 2d, we can simply enumerate the possibilities where the
              // origin may be located because edges zero and one don't share
              // any vertices, and the same for edges two and three
              if ((starting_vertex_of_edge[0] == starting_vertex_of_edge[2]) ||
                  (starting_vertex_of_edge[0] == starting_vertex_of_edge[3]))
                origin_vertex_of_cell = starting_vertex_of_edge[0];
              else if ((starting_vertex_of_edge[1] ==
                        starting_vertex_of_edge[2]) ||
                       (starting_vertex_of_edge[1] ==
                        starting_vertex_of_edge[3]))
                origin_vertex_of_cell = starting_vertex_of_edge[1];
              else
                Assert(false, ExcInternalError());

              break;
            }

          case 3:
            {
              // one could probably do something similar in 3d, but that seems
              // more complicated than one wants to write down. just go
              // through the list of possible starting vertices and check
              for (origin_vertex_of_cell = 0;
                   origin_vertex_of_cell < GeometryInfo<dim>::vertices_per_cell;
                   ++origin_vertex_of_cell)
                if (std::count(starting_vertex_of_edge,
                               starting_vertex_of_edge +
                                 GeometryInfo<dim>::lines_per_cell,
                               cell_list[cell_index]
                                 .vertex_indices[origin_vertex_of_cell]) == dim)
                  break;
              Assert(origin_vertex_of_cell <
                       GeometryInfo<dim>::vertices_per_cell,
                     ExcInternalError());

              break;
            }

          default:
            Assert(false, ExcNotImplemented());
        }

      // now rotate raw_cells[cell_index] in such a way that its orientation
      // matches that of cell_list[cell_index]
      switch (dim)
        {
          case 2:
            {
              // in 2d, we can literally rotate the cell until its origin
              // matches the one that we have determined above should be
              // the origin vertex
              //
              // when doing a rotation, take into account the ordering of
              // vertices (not in clockwise or counter-clockwise sense)
              while (raw_cells[cell_index].vertices[0] != origin_vertex_of_cell)
                {
                  const unsigned int tmp = raw_cells[cell_index].vertices[0];
                  raw_cells[cell_index].vertices[0] =
                    raw_cells[cell_index].vertices[1];
                  raw_cells[cell_index].vertices[1] =
                    raw_cells[cell_index].vertices[3];
                  raw_cells[cell_index].vertices[3] =
                    raw_cells[cell_index].vertices[2];
                  raw_cells[cell_index].vertices[2] = tmp;
                }
              break;
            }

          case 3:
            {
              // in 3d, the situation is a bit more complicated. from above, we
              // now know which vertex is at the origin (because 3 edges
              // originate from it), but that still leaves 3 possible rotations
              // of the cube. the important realization is that we can choose
              // any of them: in all 3 rotations, all edges originate from the
              // one vertex, and that fixes the directions of all 12 edges in
              // the cube because these 3 cover all 3 equivalence classes!
              // consequently, we can select an arbitrary one among the
              // permutations -- for example the following ones:
              static const unsigned int cube_permutations[8][8] = {
                {0, 1, 2, 3, 4, 5, 6, 7},
                {1, 5, 3, 7, 0, 4, 2, 6},
                {2, 6, 0, 4, 3, 7, 1, 5},
                {3, 2, 1, 0, 7, 6, 5, 4},
                {4, 0, 6, 2, 5, 1, 7, 3},
                {5, 4, 7, 6, 1, 0, 3, 2},
                {6, 7, 4, 5, 2, 3, 0, 1},
                {7, 3, 5, 1, 6, 2, 4, 0}};

              unsigned int
                temp_vertex_indices[GeometryInfo<dim>::vertices_per_cell];
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                temp_vertex_indices[v] =
                  raw_cells[cell_index]
                    .vertices[cube_permutations[origin_vertex_of_cell][v]];
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                raw_cells[cell_index].vertices[v] = temp_vertex_indices[v];

              break;
            }

          default:
            {
              Assert(false, ExcNotImplemented());
            }
        }
    }


    /**
     * Given a set of cells, find globally unique edge orientations
     * and then rotate cells so that the coordinate system of the cell
     * coincides with the coordinate systems of the adjacent edges.
     */
    template <int dim>
    void
    reorient(std::vector<CellData<dim>> &cells)
    {
      // first build the arrays that connect cells to edges and the other
      // way around
      std::vector<Edge<dim>> edge_list = build_edges(cells);
      std::vector<Cell<dim>> cell_list =
        build_cells_and_connect_edges(cells, edge_list);

      // then loop over all cells and start orienting parallel edge sets
      // of cells that still have non-oriented edges
      unsigned int next_cell_with_unoriented_edge = 0;
      while ((next_cell_with_unoriented_edge = get_next_unoriented_cell(
                cell_list, edge_list, next_cell_with_unoriented_edge)) !=
             numbers::invalid_unsigned_int)
        {
          // see which edge sets are still not oriented
          //
          // we do not need to look at each edge because if we orient edge
          // 0, we will end up with edge 1 also oriented (in 2d; in 3d, there
          // will be 3 other edges that are also oriented). there are only
          // dim independent sets of edges, so loop over these.
          //
          // we need to check whether each one of these starter edges may
          // already be oriented because the line (sheet) that connects
          // globally parallel edges may be self-intersecting in the
          // current cell
          for (unsigned int l = 0; l < dim; ++l)
            if (edge_list[cell_list[next_cell_with_unoriented_edge]
                            .edge_indices[ParallelEdges<dim>::starter_edges[l]]]
                  .orientation_status == Edge<dim>::not_oriented)
              orient_one_set_of_parallel_edges(
                cell_list,
                edge_list,
                next_cell_with_unoriented_edge,
                ParallelEdges<dim>::starter_edges[l]);

          // ensure that we have really oriented all edges now, not just
          // the starter edges
          for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
            Assert(edge_list[cell_list[next_cell_with_unoriented_edge]
                               .edge_indices[l]]
                       .orientation_status != Edge<dim>::not_oriented,
                   ExcInternalError());
        }

      // now that we have oriented all edges, we need to rotate cells
      // so that the edges point in the right direction with the now
      // rotated coordinate system
      for (unsigned int c = 0; c < cells.size(); ++c)
        rotate_cell(cell_list, edge_list, c, cells);
    }


    // overload of the function above for 1d -- there is nothing
    // to orient in that case
    void reorient(std::vector<CellData<1>> &)
    {}
  } // namespace

  template <int dim>
  void
  consistently_order_cells(std::vector<CellData<dim>> &cells)
  {
    Assert(cells.size() != 0,
           ExcMessage(
             "List of elements to orient must have at least one cell"));

    // there is nothing for us to do in 1d
    if (dim == 1)
      return;

    // check if grids are already consistent. if so, do
    // nothing. if not, then do the reordering
    if (!is_consistent(cells))
      try
        {
          reorient(cells);
        }
      catch (const ExcMeshNotOrientable &)
        {
          // the mesh is not orientable. this is acceptable if we are in 3d,
          // as class Triangulation knows how to handle this, but it is
          // not in 2d; in that case, re-throw the exception
          if (dim < 3)
            throw;
        }
  }


  // define some transformations
  namespace internal
  {
    template <int spacedim>
    class Shift
    {
    public:
      explicit Shift(const Tensor<1, spacedim> &shift)
        : shift(shift)
      {}
      Point<spacedim>
      operator()(const Point<spacedim> p) const
      {
        return p + shift;
      }

    private:
      const Tensor<1, spacedim> shift;
    };


    // Transformation to rotate around one of the cartesian axes.
    class Rotate3d
    {
    public:
      Rotate3d(const double angle, const unsigned int axis)
        : angle(angle)
        , axis(axis)
      {}

      Point<3>
      operator()(const Point<3> &p) const
      {
        if (axis == 0)
          return {p(0),
                  std::cos(angle) * p(1) - std::sin(angle) * p(2),
                  std::sin(angle) * p(1) + std::cos(angle) * p(2)};
        else if (axis == 1)
          return {std::cos(angle) * p(0) + std::sin(angle) * p(2),
                  p(1),
                  -std::sin(angle) * p(0) + std::cos(angle) * p(2)};
        else
          return {std::cos(angle) * p(0) - std::sin(angle) * p(1),
                  std::sin(angle) * p(0) + std::cos(angle) * p(1),
                  p(2)};
      }

    private:
      const double       angle;
      const unsigned int axis;
    };

    template <int spacedim>
    class Scale
    {
    public:
      explicit Scale(const double factor)
        : factor(factor)
      {}
      Point<spacedim>
      operator()(const Point<spacedim> p) const
      {
        return p * factor;
      }

    private:
      const double factor;
    };
  } // namespace internal


  template <int dim, int spacedim>
  void
  shift(const Tensor<1, spacedim> &   shift_vector,
        Triangulation<dim, spacedim> &triangulation)
  {
    transform(internal::Shift<spacedim>(shift_vector), triangulation);
  }


  template <int dim>
  void
  rotate(const double           angle,
         const unsigned int     axis,
         Triangulation<dim, 3> &triangulation)
  {
    Assert(axis < 3, ExcMessage("Invalid axis given!"));

    transform(internal::Rotate3d(angle, axis), triangulation);
  }

  template <int dim, int spacedim>
  void
  scale(const double                  scaling_factor,
        Triangulation<dim, spacedim> &triangulation)
  {
    Assert(scaling_factor > 0, ExcScalingFactorNotPositive(scaling_factor));
    transform(internal::Scale<spacedim>(scaling_factor), triangulation);
  }


  namespace internal
  {
    /**
     * Solve the Laplace equation for the @p laplace_transform function for one
     * of the @p dim space dimensions. Factorized into a function of its own
     * in order to allow parallel execution.
     */
    inline void
    laplace_solve(const SparseMatrix<double> &     S,
                  const AffineConstraints<double> &constraints,
                  Vector<double> &                 u)
    {
      const unsigned int n_dofs = S.n();
      const auto         op     = linear_operator(S);
      const auto         SF     = constrained_linear_operator(constraints, op);
      PreconditionJacobi<SparseMatrix<double>> prec;
      prec.initialize(S, 1.2);

      SolverControl                       control(n_dofs, 1.e-10, false, false);
      GrowingVectorMemory<Vector<double>> mem;
      SolverCG<Vector<double>>            solver(control, mem);

      Vector<double> f(n_dofs);

      const auto constrained_rhs =
        constrained_right_hand_side(constraints, op, f);
      solver.solve(SF, u, constrained_rhs, prec);

      constraints.distribute(u);
    }
  } // namespace internal


  // Implementation for dimensions except 1
  template <int dim>
  void
  laplace_transform(const std::map<unsigned int, Point<dim>> &new_points,
                    Triangulation<dim> &                      triangulation,
                    const Function<dim> *                     coefficient,
                    const bool solve_for_absolute_positions)
  {
    if (dim == 1)
      Assert(false, ExcNotImplemented());

    // first provide everything that is needed for solving a Laplace
    // equation.
    FE_Q<dim> q1(1);

    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(q1);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    dsp.compress();

    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.compress();

    SparseMatrix<double> S(sparsity_pattern);

    QGauss<dim> quadrature(4);

    MatrixCreator::create_laplace_matrix(
      StaticMappingQ1<dim>::mapping, dof_handler, quadrature, S, coefficient);

    // set up the boundary values for the laplace problem
    std::array<AffineConstraints<double>, dim>                  constraints;
    typename std::map<unsigned int, Point<dim>>::const_iterator map_end =
      new_points.end();

    // fill these maps using the data given by new_points
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        // loop over all vertices of the cell and see if it is listed in the map
        // given as first argument of the function
        for (const unsigned int vertex_no : cell->vertex_indices())
          {
            const unsigned int vertex_index = cell->vertex_index(vertex_no);
            const Point<dim> & vertex_point = cell->vertex(vertex_no);

            const typename std::map<unsigned int, Point<dim>>::const_iterator
              map_iter = new_points.find(vertex_index);

            if (map_iter != map_end)
              for (unsigned int i = 0; i < dim; ++i)
                {
                  constraints[i].add_line(cell->vertex_dof_index(vertex_no, 0));
                  constraints[i].set_inhomogeneity(
                    cell->vertex_dof_index(vertex_no, 0),
                    (solve_for_absolute_positions ?
                       map_iter->second(i) :
                       map_iter->second(i) - vertex_point[i]));
                }
          }
      }

    for (unsigned int i = 0; i < dim; ++i)
      constraints[i].close();

    // solve the dim problems with different right hand sides.
    Vector<double> us[dim];
    for (unsigned int i = 0; i < dim; ++i)
      us[i].reinit(dof_handler.n_dofs());

    // solve linear systems in parallel
    Threads::TaskGroup<> tasks;
    for (unsigned int i = 0; i < dim; ++i)
      tasks +=
        Threads::new_task(&internal::laplace_solve, S, constraints[i], us[i]);
    tasks.join_all();

    // change the coordinates of the points of the triangulation
    // according to the computed values
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
    for (const auto &cell : dof_handler.active_cell_iterators())
      for (const unsigned int vertex_no : cell->vertex_indices())
        if (vertex_touched[cell->vertex_index(vertex_no)] == false)
          {
            Point<dim> &v = cell->vertex(vertex_no);

            const types::global_dof_index dof_index =
              cell->vertex_dof_index(vertex_no, 0);
            for (unsigned int i = 0; i < dim; ++i)
              if (solve_for_absolute_positions)
                v(i) = us[i](dof_index);
              else
                v(i) += us[i](dof_index);

            vertex_touched[cell->vertex_index(vertex_no)] = true;
          }
  }

  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  get_all_vertices_at_boundary(const Triangulation<dim, spacedim> &tria)
  {
    std::map<unsigned int, Point<spacedim>> vertex_map;
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = tria.begin_active(),
      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        for (unsigned int i : cell->face_indices())
          {
            const typename Triangulation<dim, spacedim>::face_iterator &face =
              cell->face(i);
            if (face->at_boundary())
              {
                for (unsigned j = 0; j < face->n_vertices(); ++j)
                  {
                    const Point<spacedim> &vertex       = face->vertex(j);
                    const unsigned int     vertex_index = face->vertex_index(j);
                    vertex_map[vertex_index]            = vertex;
                  }
              }
          }
      }
    return vertex_map;
  }

  /**
   * Distort a triangulation in
   * some random way.
   */
  template <int dim, int spacedim>
  void
  distort_random(const double                  factor,
                 Triangulation<dim, spacedim> &triangulation,
                 const bool                    keep_boundary,
                 const unsigned int            seed)
  {
    // if spacedim>dim we need to make sure that we perturb
    // points but keep them on
    // the manifold. however, this isn't implemented right now
    Assert(spacedim == dim, ExcNotImplemented());


    // find the smallest length of the
    // lines adjacent to the
    // vertex. take the initial value
    // to be larger than anything that
    // might be found: the diameter of
    // the triangulation, here
    // estimated by adding up the
    // diameters of the coarse grid
    // cells.
    double almost_infinite_length = 0;
    for (typename Triangulation<dim, spacedim>::cell_iterator cell =
           triangulation.begin(0);
         cell != triangulation.end(0);
         ++cell)
      almost_infinite_length += cell->diameter();

    std::vector<double> minimal_length(triangulation.n_vertices(),
                                       almost_infinite_length);

    // also note if a vertex is at the boundary
    std::vector<bool> at_boundary(keep_boundary ? triangulation.n_vertices() :
                                                  0,
                                  false);
    // for parallel::shared::Triangulation we need to work on all vertices,
    // not just the ones related to locally owned cells;
    const bool is_parallel_shared =
      (dynamic_cast<parallel::shared::Triangulation<dim, spacedim> *>(
         &triangulation) != nullptr);
    for (const auto &cell : triangulation.active_cell_iterators())
      if (is_parallel_shared || cell->is_locally_owned())
        {
          if (dim > 1)
            {
              for (unsigned int i = 0; i < cell->n_lines(); ++i)
                {
                  const typename Triangulation<dim, spacedim>::line_iterator
                    line = cell->line(i);

                  if (keep_boundary && line->at_boundary())
                    {
                      at_boundary[line->vertex_index(0)] = true;
                      at_boundary[line->vertex_index(1)] = true;
                    }

                  minimal_length[line->vertex_index(0)] =
                    std::min(line->diameter(),
                             minimal_length[line->vertex_index(0)]);
                  minimal_length[line->vertex_index(1)] =
                    std::min(line->diameter(),
                             minimal_length[line->vertex_index(1)]);
                }
            }
          else // dim==1
            {
              if (keep_boundary)
                for (unsigned int vertex = 0; vertex < 2; ++vertex)
                  if (cell->at_boundary(vertex) == true)
                    at_boundary[cell->vertex_index(vertex)] = true;

              minimal_length[cell->vertex_index(0)] =
                std::min(cell->diameter(),
                         minimal_length[cell->vertex_index(0)]);
              minimal_length[cell->vertex_index(1)] =
                std::min(cell->diameter(),
                         minimal_length[cell->vertex_index(1)]);
            }
        }

    // create a random number generator for the interval [-1,1]
    boost::random::mt19937                     rng(seed);
    boost::random::uniform_real_distribution<> uniform_distribution(-1, 1);

    // If the triangulation is distributed, we need to
    // exchange the moved vertices across mpi processes
    if (auto distributed_triangulation =
          dynamic_cast<parallel::DistributedTriangulationBase<dim, spacedim> *>(
            &triangulation))
      {
        const std::vector<bool> locally_owned_vertices =
          get_locally_owned_vertices(triangulation);
        std::vector<bool> vertex_moved(triangulation.n_vertices(), false);

        // Next move vertices on locally owned cells
        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              for (const unsigned int vertex_no : cell->vertex_indices())
                {
                  const unsigned global_vertex_no =
                    cell->vertex_index(vertex_no);

                  // ignore this vertex if we shall keep the boundary and
                  // this vertex *is* at the boundary, if it is already moved
                  // or if another process moves this vertex
                  if ((keep_boundary && at_boundary[global_vertex_no]) ||
                      vertex_moved[global_vertex_no] ||
                      !locally_owned_vertices[global_vertex_no])
                    continue;

                  // first compute a random shift vector
                  Point<spacedim> shift_vector;
                  for (unsigned int d = 0; d < spacedim; ++d)
                    shift_vector(d) = uniform_distribution(rng);

                  shift_vector *= factor * minimal_length[global_vertex_no] /
                                  std::sqrt(shift_vector.square());

                  // finally move the vertex
                  cell->vertex(vertex_no) += shift_vector;
                  vertex_moved[global_vertex_no] = true;
                }
            }

        distributed_triangulation->communicate_locally_moved_vertices(
          locally_owned_vertices);
      }
    else
      // if this is a sequential triangulation, we could in principle
      // use the algorithm above, but we'll use an algorithm that we used
      // before the parallel::distributed::Triangulation was introduced
      // in order to preserve backward compatibility
      {
        // loop over all vertices and compute their new locations
        const unsigned int           n_vertices = triangulation.n_vertices();
        std::vector<Point<spacedim>> new_vertex_locations(n_vertices);
        const std::vector<Point<spacedim>> &old_vertex_locations =
          triangulation.get_vertices();

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          {
            // ignore this vertex if we will keep the boundary and
            // this vertex *is* at the boundary
            if (keep_boundary && at_boundary[vertex])
              new_vertex_locations[vertex] = old_vertex_locations[vertex];
            else
              {
                // compute a random shift vector
                Point<spacedim> shift_vector;
                for (unsigned int d = 0; d < spacedim; ++d)
                  shift_vector(d) = uniform_distribution(rng);

                shift_vector *= factor * minimal_length[vertex] /
                                std::sqrt(shift_vector.square());

                // record new vertex location
                new_vertex_locations[vertex] =
                  old_vertex_locations[vertex] + shift_vector;
              }
          }

        // now do the actual move of the vertices
        for (const auto &cell : triangulation.active_cell_iterators())
          for (const unsigned int vertex_no : cell->vertex_indices())
            cell->vertex(vertex_no) =
              new_vertex_locations[cell->vertex_index(vertex_no)];
      }

    // Correct hanging nodes if necessary
    if (dim >= 2)
      {
        // We do the same as in GridTools::transform
        //
        // exclude hanging nodes at the boundaries of artificial cells:
        // these may belong to ghost cells for which we know the exact
        // location of vertices, whereas the artificial cell may or may
        // not be further refined, and so we cannot know whether
        // the location of the hanging node is correct or not
        typename Triangulation<dim, spacedim>::active_cell_iterator
          cell = triangulation.begin_active(),
          endc = triangulation.end();
        for (; cell != endc; ++cell)
          if (!cell->is_artificial())
            for (const unsigned int face : cell->face_indices())
              if (cell->face(face)->has_children() &&
                  !cell->face(face)->at_boundary())
                {
                  // this face has hanging nodes
                  if (dim == 2)
                    cell->face(face)->child(0)->vertex(1) =
                      (cell->face(face)->vertex(0) +
                       cell->face(face)->vertex(1)) /
                      2;
                  else if (dim == 3)
                    {
                      cell->face(face)->child(0)->vertex(1) =
                        .5 * (cell->face(face)->vertex(0) +
                              cell->face(face)->vertex(1));
                      cell->face(face)->child(0)->vertex(2) =
                        .5 * (cell->face(face)->vertex(0) +
                              cell->face(face)->vertex(2));
                      cell->face(face)->child(1)->vertex(3) =
                        .5 * (cell->face(face)->vertex(1) +
                              cell->face(face)->vertex(3));
                      cell->face(face)->child(2)->vertex(3) =
                        .5 * (cell->face(face)->vertex(2) +
                              cell->face(face)->vertex(3));

                      // center of the face
                      cell->face(face)->child(0)->vertex(3) =
                        .25 * (cell->face(face)->vertex(0) +
                               cell->face(face)->vertex(1) +
                               cell->face(face)->vertex(2) +
                               cell->face(face)->vertex(3));
                    }
                }
      }
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  unsigned int
  find_closest_vertex(const MeshType<dim, spacedim> &mesh,
                      const Point<spacedim> &        p,
                      const std::vector<bool> &      marked_vertices)
  {
    // first get the underlying triangulation from the mesh and determine
    // vertices and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const std::vector<Point<spacedim>> &vertices = tria.get_vertices();

    Assert(tria.get_vertices().size() == marked_vertices.size() ||
             marked_vertices.size() == 0,
           ExcDimensionMismatch(tria.get_vertices().size(),
                                marked_vertices.size()));

    // marked_vertices is expected to be a subset of used_vertices. Thus,
    // comparing the range marked_vertices.begin() to marked_vertices.end() with
    // the range used_vertices.begin() to used_vertices.end() the element in the
    // second range must be valid if the element in the first range is valid.
    Assert(
      marked_vertices.size() == 0 ||
        std::equal(marked_vertices.begin(),
                   marked_vertices.end(),
                   tria.get_used_vertices().begin(),
                   [](bool p, bool q) { return !p || q; }),
      ExcMessage(
        "marked_vertices should be a subset of used vertices in the triangulation "
        "but marked_vertices contains one or more vertices that are not used vertices!"));

    // If marked_indices is empty, consider all used_vertices for finding the
    // closest vertex to the point. Otherwise, marked_indices is used.
    const std::vector<bool> &vertices_to_use = (marked_vertices.size() == 0) ?
                                                 tria.get_used_vertices() :
                                                 marked_vertices;

    // At the beginning, the first used vertex is considered to be the closest
    // one.
    std::vector<bool>::const_iterator first =
      std::find(vertices_to_use.begin(), vertices_to_use.end(), true);

    // Assert that at least one vertex is actually used
    Assert(first != vertices_to_use.end(), ExcInternalError());

    unsigned int best_vertex = std::distance(vertices_to_use.begin(), first);
    double       best_dist   = (p - vertices[best_vertex]).norm_square();

    // For all remaining vertices, test
    // whether they are any closer
    for (unsigned int j = best_vertex + 1; j < vertices.size(); j++)
      if (vertices_to_use[j])
        {
          const double dist = (p - vertices[j]).norm_square();
          if (dist < best_dist)
            {
              best_vertex = j;
              best_dist   = dist;
            }
        }

    return best_vertex;
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  unsigned int
  find_closest_vertex(const Mapping<dim, spacedim> & mapping,
                      const MeshType<dim, spacedim> &mesh,
                      const Point<spacedim> &        p,
                      const std::vector<bool> &      marked_vertices)
  {
    // Take a shortcut in the simple case.
    if (mapping.preserves_vertex_locations() == true)
      return find_closest_vertex(mesh, p, marked_vertices);

    // first get the underlying triangulation from the mesh and determine
    // vertices and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    auto vertices = extract_used_vertices(tria, mapping);

    Assert(tria.get_vertices().size() == marked_vertices.size() ||
             marked_vertices.size() == 0,
           ExcDimensionMismatch(tria.get_vertices().size(),
                                marked_vertices.size()));

    // marked_vertices is expected to be a subset of used_vertices. Thus,
    // comparing the range marked_vertices.begin() to marked_vertices.end()
    // with the range used_vertices.begin() to used_vertices.end() the element
    // in the second range must be valid if the element in the first range is
    // valid.
    Assert(
      marked_vertices.size() == 0 ||
        std::equal(marked_vertices.begin(),
                   marked_vertices.end(),
                   tria.get_used_vertices().begin(),
                   [](bool p, bool q) { return !p || q; }),
      ExcMessage(
        "marked_vertices should be a subset of used vertices in the triangulation "
        "but marked_vertices contains one or more vertices that are not used vertices!"));

    // Remove from the map unwanted elements.
    if (marked_vertices.size() != 0)
      for (auto it = vertices.begin(); it != vertices.end();)
        {
          if (marked_vertices[it->first] == false)
            {
              it = vertices.erase(it);
            }
          else
            {
              ++it;
            }
        }

    return find_closest_vertex(vertices, p);
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  std::vector<typename MeshType<dim, spacedim>::active_cell_iterator>
#else
  std::vector<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type>
#endif
  find_cells_adjacent_to_vertex(const MeshType<dim, spacedim> &mesh,
                                const unsigned int             vertex)
  {
    // make sure that the given vertex is
    // an active vertex of the underlying
    // triangulation
    AssertIndexRange(vertex, mesh.get_triangulation().n_vertices());
    Assert(mesh.get_triangulation().get_used_vertices()[vertex],
           ExcVertexNotUsed(vertex));

    // use a set instead of a vector
    // to ensure that cells are inserted only
    // once
    std::set<typename dealii::internal::
               ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type>
      adjacent_cells;

    // go through all active cells and look if the vertex is part of that cell
    //
    // in 1d, this is all we need to care about. in 2d/3d we also need to worry
    // that the vertex might be a hanging node on a face or edge of a cell; in
    // this case, we would want to add those cells as well on whose faces the
    // vertex is located but for which it is not a vertex itself.
    //
    // getting this right is a lot simpler in 2d than in 3d. in 2d, a hanging
    // node can only be in the middle of a face and we can query the neighboring
    // cell from the current cell. on the other hand, in 3d a hanging node
    // vertex can also be on an edge but there can be many other cells on
    // this edge and we can not access them from the cell we are currently
    // on.
    //
    // so, in the 3d case, if we run the algorithm as in 2d, we catch all
    // those cells for which the vertex we seek is on a *subface*, but we
    // miss the case of cells for which the vertex we seek is on a
    // sub-edge for which there is no corresponding sub-face (because the
    // immediate neighbor behind this face is not refined), see for example
    // the bits/find_cells_adjacent_to_vertex_6 testcase. thus, if we
    // haven't yet found the vertex for the current cell we also need to
    // look at the mid-points of edges
    //
    // as a final note, deciding whether a neighbor is actually coarser is
    // simple in the case of isotropic refinement (we just need to look at
    // the level of the current and the neighboring cell). however, this
    // isn't so simple if we have used anisotropic refinement since then
    // the level of a cell is not indicative of whether it is coarser or
    // not than the current cell. ultimately, we want to add all cells on
    // which the vertex is, independent of whether they are coarser or
    // finer and so in the 2d case below we simply add *any* *active* neighbor.
    // in the worst case, we add cells multiple times to the adjacent_cells
    // list, but std::set throws out those cells already entered
    for (const auto &cell : mesh.active_cell_iterators())
      {
        for (const unsigned int v : cell->vertex_indices())
          if (cell->vertex_index(v) == vertex)
            {
              // OK, we found a cell that contains
              // the given vertex. We add it
              // to the list.
              adjacent_cells.insert(cell);

              // as explained above, in 2+d we need to check whether
              // this vertex is on a face behind which there is a
              // (possibly) coarser neighbor. if this is the case,
              // then we need to also add this neighbor
              if (dim >= 2)
                for (const auto face :
                     cell->reference_cell().faces_for_given_vertex(v))
                  if (!cell->at_boundary(face) &&
                      cell->neighbor(face)->is_active())
                    {
                      // there is a (possibly) coarser cell behind a
                      // face to which the vertex belongs. the
                      // vertex we are looking at is then either a
                      // vertex of that coarser neighbor, or it is a
                      // hanging node on one of the faces of that
                      // cell. in either case, it is adjacent to the
                      // vertex, so add it to the list as well (if
                      // the cell was already in the list then the
                      // std::set makes sure that we get it only
                      // once)
                      adjacent_cells.insert(cell->neighbor(face));
                    }

              // in any case, we have found a cell, so go to the next cell
              goto next_cell;
            }

        // in 3d also loop over the edges
        if (dim >= 3)
          {
            for (unsigned int e = 0; e < cell->n_lines(); ++e)
              if (cell->line(e)->has_children())
                // the only place where this vertex could have been
                // hiding is on the mid-edge point of the edge we
                // are looking at
                if (cell->line(e)->child(0)->vertex_index(1) == vertex)
                  {
                    adjacent_cells.insert(cell);

                    // jump out of this tangle of nested loops
                    goto next_cell;
                  }
          }

        // in more than 3d we would probably have to do the same as
        // above also for even lower-dimensional objects
        Assert(dim <= 3, ExcNotImplemented());

        // move on to the next cell if we have found the
        // vertex on the current one
      next_cell:;
      }

    // if this was an active vertex then there needs to have been
    // at least one cell to which it is adjacent!
    Assert(adjacent_cells.size() > 0, ExcInternalError());

    // return the result as a vector, rather than the set we built above
    return std::vector<
      typename dealii::internal::
        ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type>(
      adjacent_cells.begin(), adjacent_cells.end());
  }



  template <int dim, int spacedim>
  std::vector<std::vector<Tensor<1, spacedim>>>
  vertex_to_cell_centers_directions(
    const Triangulation<dim, spacedim> &mesh,
    const std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      &vertex_to_cells)
  {
    const std::vector<Point<spacedim>> &vertices   = mesh.get_vertices();
    const unsigned int                  n_vertices = vertex_to_cells.size();

    AssertDimension(vertices.size(), n_vertices);


    std::vector<std::vector<Tensor<1, spacedim>>> vertex_to_cell_centers(
      n_vertices);
    for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
      if (mesh.vertex_used(vertex))
        {
          const unsigned int n_neighbor_cells = vertex_to_cells[vertex].size();
          vertex_to_cell_centers[vertex].resize(n_neighbor_cells);

          typename std::set<typename Triangulation<dim, spacedim>::
                              active_cell_iterator>::iterator it =
            vertex_to_cells[vertex].begin();
          for (unsigned int cell = 0; cell < n_neighbor_cells; ++cell, ++it)
            {
              vertex_to_cell_centers[vertex][cell] =
                (*it)->center() - vertices[vertex];
              vertex_to_cell_centers[vertex][cell] /=
                vertex_to_cell_centers[vertex][cell].norm();
            }
        }
    return vertex_to_cell_centers;
  }


  namespace internal
  {
    template <int spacedim>
    bool
    compare_point_association(
      const unsigned int                      a,
      const unsigned int                      b,
      const Tensor<1, spacedim> &             point_direction,
      const std::vector<Tensor<1, spacedim>> &center_directions)
    {
      const double scalar_product_a = center_directions[a] * point_direction;
      const double scalar_product_b = center_directions[b] * point_direction;

      // The function is supposed to return if a is before b. We are looking
      // for the alignment of point direction and center direction, therefore
      // return if the scalar product of a is larger.
      return (scalar_product_a > scalar_product_b);
    }
  } // namespace internal

  template <int dim, template <int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#endif
  find_active_cell_around_point(
    const Mapping<dim, spacedim> & mapping,
    const MeshType<dim, spacedim> &mesh,
    const Point<spacedim> &        p,
    const std::vector<
      std::set<typename MeshType<dim, spacedim>::active_cell_iterator>>
      &                                                  vertex_to_cells,
    const std::vector<std::vector<Tensor<1, spacedim>>> &vertex_to_cell_centers,
    const typename MeshType<dim, spacedim>::active_cell_iterator &cell_hint,
    const std::vector<bool> &                              marked_vertices,
    const RTree<std::pair<Point<spacedim>, unsigned int>> &used_vertices_rtree,
    const double                                           tolerance)
  {
    std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
              Point<dim>>
      cell_and_position;
    // To handle points at the border we keep track of points which are close to
    // the unit cell:
    std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
              Point<dim>>
      cell_and_position_approx;

    bool found_cell  = false;
    bool approx_cell = false;

    unsigned int        closest_vertex_index = 0;
    Tensor<1, spacedim> vertex_to_point;
    auto                current_cell = cell_hint;

    while (found_cell == false)
      {
        // First look at the vertices of the cell cell_hint. If it's an
        // invalid cell, then query for the closest global vertex
        if (current_cell.state() == IteratorState::valid)
          {
            const auto cell_vertices = mapping.get_vertices(current_cell);
            const unsigned int closest_vertex =
              find_closest_vertex_of_cell<dim, spacedim>(current_cell,
                                                         p,
                                                         mapping);
            vertex_to_point      = p - cell_vertices[closest_vertex];
            closest_vertex_index = current_cell->vertex_index(closest_vertex);
          }
        else
          {
            if (!used_vertices_rtree.empty())
              {
                // If we have an rtree at our disposal, use it.
                using ValueType = std::pair<Point<spacedim>, unsigned int>;
                std::function<bool(const ValueType &)> marked;
                if (marked_vertices.size() == mesh.n_vertices())
                  marked = [&marked_vertices](const ValueType &value) -> bool {
                    return marked_vertices[value.second];
                  };
                else
                  marked = [](const ValueType &) -> bool { return true; };

                std::vector<std::pair<Point<spacedim>, unsigned int>> res;
                used_vertices_rtree.query(
                  boost::geometry::index::nearest(p, 1) &&
                    boost::geometry::index::satisfies(marked),
                  std::back_inserter(res));

                // We should have one and only one result
                AssertDimension(res.size(), 1);
                closest_vertex_index = res[0].second;
              }
            else
              {
                closest_vertex_index = GridTools::find_closest_vertex(
                  mapping, mesh, p, marked_vertices);
              }
            vertex_to_point = p - mesh.get_vertices()[closest_vertex_index];
          }

        const double vertex_point_norm = vertex_to_point.norm();
        if (vertex_point_norm > 0)
          vertex_to_point /= vertex_point_norm;

        const unsigned int n_neighbor_cells =
          vertex_to_cells[closest_vertex_index].size();

        // Create a corresponding map of vectors from vertex to cell center
        std::vector<unsigned int> neighbor_permutation(n_neighbor_cells);

        for (unsigned int i = 0; i < n_neighbor_cells; ++i)
          neighbor_permutation[i] = i;

        auto comp = [&](const unsigned int a, const unsigned int b) -> bool {
          return internal::compare_point_association<spacedim>(
            a,
            b,
            vertex_to_point,
            vertex_to_cell_centers[closest_vertex_index]);
        };

        std::sort(neighbor_permutation.begin(),
                  neighbor_permutation.end(),
                  comp);
        // It is possible the vertex is close
        // to an edge, thus we add a tolerance
        // to keep also the "best" cell
        double best_distance = tolerance;

        // Search all of the cells adjacent to the closest vertex of the cell
        // hint Most likely we will find the point in them.
        for (unsigned int i = 0; i < n_neighbor_cells; ++i)
          {
            try
              {
                auto cell = vertex_to_cells[closest_vertex_index].begin();
                std::advance(cell, neighbor_permutation[i]);

                if (!(*cell)->is_artificial())
                  {
                    const Point<dim> p_unit =
                      mapping.transform_real_to_unit_cell(*cell, p);
                    if (GeometryInfo<dim>::is_inside_unit_cell(p_unit,
                                                               tolerance))
                      {
                        cell_and_position.first  = *cell;
                        cell_and_position.second = p_unit;
                        found_cell               = true;
                        approx_cell              = false;
                        break;
                      }
                    // The point is not inside this cell: checking how far
                    // outside it is and whether we want to use this cell as a
                    // backup if we can't find a cell within which the point
                    // lies.
                    const double dist =
                      GeometryInfo<dim>::distance_to_unit_cell(p_unit);
                    if (dist < best_distance)
                      {
                        best_distance                   = dist;
                        cell_and_position_approx.first  = *cell;
                        cell_and_position_approx.second = p_unit;
                        approx_cell                     = true;
                      }
                  }
              }
            catch (typename Mapping<dim>::ExcTransformationFailed &)
              {}
          }

        if (found_cell == true)
          return cell_and_position;
        else if (approx_cell == true)
          return cell_and_position_approx;

        // The first time around, we check for vertices in the hint_cell. If
        // that does not work, we set the cell iterator to an invalid one, and
        // look for a global vertex close to the point. If that does not work,
        // we are in trouble, and just throw an exception.
        //
        // If we got here, then we did not find the point. If the
        // current_cell.state() here is not IteratorState::valid, it means that
        // the user did not provide a hint_cell, and at the beginning of the
        // while loop we performed an actual global search on the mesh
        // vertices. Not finding the point then means the point is outside the
        // domain, or that we've had problems with the algorithm above. Try as a
        // last resort the other (simpler) algorithm.
        if (current_cell.state() != IteratorState::valid)
          return find_active_cell_around_point(
            mapping, mesh, p, marked_vertices, tolerance);

        current_cell = typename MeshType<dim, spacedim>::active_cell_iterator();
      }
    return cell_and_position;
  }



  template <int dim, int spacedim>
  unsigned int
  find_closest_vertex_of_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
    const Point<spacedim> &                                            position,
    const Mapping<dim, spacedim> &                                     mapping)
  {
    const auto   vertices         = mapping.get_vertices(cell);
    double       minimum_distance = position.distance_square(vertices[0]);
    unsigned int closest_vertex   = 0;

    for (unsigned int v = 1; v < cell->n_vertices(); ++v)
      {
        const double vertex_distance = position.distance_square(vertices[v]);
        if (vertex_distance < minimum_distance)
          {
            closest_vertex   = v;
            minimum_distance = vertex_distance;
          }
      }
    return closest_vertex;
  }



  namespace internal
  {
    namespace BoundingBoxPredicate
    {
      template <class MeshType>
      std::tuple<BoundingBox<MeshType::space_dimension>, bool>
      compute_cell_predicate_bounding_box(
        const typename MeshType::cell_iterator &parent_cell,
        const std::function<
          bool(const typename MeshType::active_cell_iterator &)> &predicate)
      {
        bool has_predicate =
          false; // Start assuming there's no cells with predicate inside
        std::vector<typename MeshType::active_cell_iterator> active_cells;
        if (parent_cell->is_active())
          active_cells = {parent_cell};
        else
          // Finding all active cells descendants of the current one (or the
          // current one if it is active)
          active_cells = get_active_child_cells<MeshType>(parent_cell);

        const unsigned int spacedim = MeshType::space_dimension;

        // Looking for the first active cell which has the property predicate
        unsigned int i = 0;
        while (i < active_cells.size() && !predicate(active_cells[i]))
          ++i;

        // No active cells or no active cells with property
        if (active_cells.size() == 0 || i == active_cells.size())
          {
            BoundingBox<spacedim> bbox;
            return std::make_tuple(bbox, has_predicate);
          }

        // The two boundary points defining the boundary box
        Point<spacedim> maxp = active_cells[i]->vertex(0);
        Point<spacedim> minp = active_cells[i]->vertex(0);

        for (; i < active_cells.size(); ++i)
          if (predicate(active_cells[i]))
            for (const unsigned int v : active_cells[i]->vertex_indices())
              for (unsigned int d = 0; d < spacedim; ++d)
                {
                  minp[d] = std::min(minp[d], active_cells[i]->vertex(v)[d]);
                  maxp[d] = std::max(maxp[d], active_cells[i]->vertex(v)[d]);
                }

        has_predicate = true;
        BoundingBox<spacedim> bbox(std::make_pair(minp, maxp));
        return std::make_tuple(bbox, has_predicate);
      }
    } // namespace BoundingBoxPredicate
  }   // namespace internal



  template <class MeshType>
  std::vector<BoundingBox<MeshType::space_dimension>>
  compute_mesh_predicate_bounding_box(
    const MeshType &mesh,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &                predicate,
    const unsigned int refinement_level,
    const bool         allow_merge,
    const unsigned int max_boxes)
  {
    // Algorithm brief description: begin with creating bounding boxes of all
    // cells at refinement_level (and coarser levels if there are active cells)
    // which have the predicate property. These are then merged

    Assert(
      refinement_level <= mesh.n_levels(),
      ExcMessage(
        "Error: refinement level is higher then total levels in the triangulation!"));

    const unsigned int                 spacedim = MeshType::space_dimension;
    std::vector<BoundingBox<spacedim>> bounding_boxes;

    // Creating a bounding box for all active cell on coarser level

    for (unsigned int i = 0; i < refinement_level; ++i)
      for (const typename MeshType::cell_iterator &cell :
           mesh.active_cell_iterators_on_level(i))
        {
          bool                  has_predicate = false;
          BoundingBox<spacedim> bbox;
          std::tie(bbox, has_predicate) =
            internal::BoundingBoxPredicate::compute_cell_predicate_bounding_box<
              MeshType>(cell, predicate);
          if (has_predicate)
            bounding_boxes.push_back(bbox);
        }

    // Creating a Bounding Box for all cells on the chosen refinement_level
    for (const typename MeshType::cell_iterator &cell :
         mesh.cell_iterators_on_level(refinement_level))
      {
        bool                  has_predicate = false;
        BoundingBox<spacedim> bbox;
        std::tie(bbox, has_predicate) =
          internal::BoundingBoxPredicate::compute_cell_predicate_bounding_box<
            MeshType>(cell, predicate);
        if (has_predicate)
          bounding_boxes.push_back(bbox);
      }

    if (!allow_merge)
      // If merging is not requested return the created bounding_boxes
      return bounding_boxes;
    else
      {
        // Merging part of the algorithm
        // Part 1: merging neighbors
        // This array stores the indices of arrays we have already merged
        std::vector<unsigned int> merged_boxes_idx;
        bool                      found_neighbors = true;

        // We merge only neighbors which can be expressed by a single bounding
        // box e.g. in 1d [0,1] and [1,2] can be described with [0,2] without
        // losing anything
        while (found_neighbors)
          {
            found_neighbors = false;
            for (unsigned int i = 0; i < bounding_boxes.size() - 1; ++i)
              {
                if (std::find(merged_boxes_idx.begin(),
                              merged_boxes_idx.end(),
                              i) == merged_boxes_idx.end())
                  for (unsigned int j = i + 1; j < bounding_boxes.size(); ++j)
                    if (std::find(merged_boxes_idx.begin(),
                                  merged_boxes_idx.end(),
                                  j) == merged_boxes_idx.end() &&
                        bounding_boxes[i].get_neighbor_type(
                          bounding_boxes[j]) ==
                          NeighborType::mergeable_neighbors)
                      {
                        bounding_boxes[i].merge_with(bounding_boxes[j]);
                        merged_boxes_idx.push_back(j);
                        found_neighbors = true;
                      }
              }
          }

        // Copying the merged boxes into merged_b_boxes
        std::vector<BoundingBox<spacedim>> merged_b_boxes;
        for (unsigned int i = 0; i < bounding_boxes.size(); ++i)
          if (std::find(merged_boxes_idx.begin(), merged_boxes_idx.end(), i) ==
              merged_boxes_idx.end())
            merged_b_boxes.push_back(bounding_boxes[i]);

        // Part 2: if there are too many bounding boxes, merging smaller boxes
        // This has sense only in dimension 2 or greater, since  in dimension 1,
        // neighboring intervals can always be merged without problems
        if ((merged_b_boxes.size() > max_boxes) && (spacedim > 1))
          {
            std::vector<double> volumes;
            for (unsigned int i = 0; i < merged_b_boxes.size(); ++i)
              volumes.push_back(merged_b_boxes[i].volume());

            while (merged_b_boxes.size() > max_boxes)
              {
                unsigned int min_idx =
                  std::min_element(volumes.begin(), volumes.end()) -
                  volumes.begin();
                volumes.erase(volumes.begin() + min_idx);
                // Finding a neighbor
                bool not_removed = true;
                for (unsigned int i = 0;
                     i < merged_b_boxes.size() && not_removed;
                     ++i)
                  // We merge boxes if we have "attached" or "mergeable"
                  // neighbors, even though mergeable should be dealt with in
                  // Part 1
                  if (i != min_idx && (merged_b_boxes[i].get_neighbor_type(
                                         merged_b_boxes[min_idx]) ==
                                         NeighborType::attached_neighbors ||
                                       merged_b_boxes[i].get_neighbor_type(
                                         merged_b_boxes[min_idx]) ==
                                         NeighborType::mergeable_neighbors))
                    {
                      merged_b_boxes[i].merge_with(merged_b_boxes[min_idx]);
                      merged_b_boxes.erase(merged_b_boxes.begin() + min_idx);
                      not_removed = false;
                    }
                Assert(!not_removed,
                       ExcMessage("Error: couldn't merge bounding boxes!"));
              }
          }
        Assert(merged_b_boxes.size() <= max_boxes,
               ExcMessage(
                 "Error: couldn't reach target number of bounding boxes!"));
        return merged_b_boxes;
      }
  }



  template <int spacedim>
#ifndef DOXYGEN
  std::tuple<std::vector<std::vector<unsigned int>>,
             std::map<unsigned int, unsigned int>,
             std::map<unsigned int, std::vector<unsigned int>>>
#else
  return_type
#endif
  guess_point_owner(
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const std::vector<Point<spacedim>> &                   points)
  {
    unsigned int                           n_procs = global_bboxes.size();
    std::vector<std::vector<unsigned int>> point_owners(n_procs);
    std::map<unsigned int, unsigned int>   map_owners_found;
    std::map<unsigned int, std::vector<unsigned int>> map_owners_guessed;

    unsigned int n_points = points.size();
    for (unsigned int pt = 0; pt < n_points; ++pt)
      {
        // Keep track of how many processes we guess to own the point
        std::vector<unsigned int> owners_found;
        // Check in which other processes the point might be
        for (unsigned int rk = 0; rk < n_procs; ++rk)
          {
            for (const BoundingBox<spacedim> &bbox : global_bboxes[rk])
              if (bbox.point_inside(points[pt]))
                {
                  point_owners[rk].emplace_back(pt);
                  owners_found.emplace_back(rk);
                  break; // We can check now the next process
                }
          }
        Assert(owners_found.size() > 0,
               ExcMessage("No owners found for the point " +
                          std::to_string(pt)));
        if (owners_found.size() == 1)
          map_owners_found[pt] = owners_found[0];
        else
          // Multiple owners
          map_owners_guessed[pt] = owners_found;
      }

    return std::make_tuple(std::move(point_owners),
                           std::move(map_owners_found),
                           std::move(map_owners_guessed));
  }

  template <int spacedim>
#ifndef DOXYGEN
  std::tuple<std::map<unsigned int, std::vector<unsigned int>>,
             std::map<unsigned int, unsigned int>,
             std::map<unsigned int, std::vector<unsigned int>>>
#else
  return_type
#endif
  guess_point_owner(
    const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &covering_rtree,
    const std::vector<Point<spacedim>> &                         points)
  {
    std::map<unsigned int, std::vector<unsigned int>> point_owners;
    std::map<unsigned int, unsigned int>              map_owners_found;
    std::map<unsigned int, std::vector<unsigned int>> map_owners_guessed;
    std::vector<std::pair<BoundingBox<spacedim>, unsigned int>> search_result;

    unsigned int n_points = points.size();
    for (unsigned int pt_n = 0; pt_n < n_points; ++pt_n)
      {
        search_result.clear(); // clearing last output

        // Running tree search
        covering_rtree.query(boost::geometry::index::intersects(points[pt_n]),
                             std::back_inserter(search_result));

        // Keep track of how many processes we guess to own the point
        std::set<unsigned int> owners_found;
        // Check in which other processes the point might be
        for (const auto &rank_bbox : search_result)
          {
            // Try to add the owner to the owners found,
            // and check if it was already present
            const bool pt_inserted = owners_found.insert(pt_n).second;
            if (pt_inserted)
              point_owners[rank_bbox.second].emplace_back(pt_n);
          }
        Assert(owners_found.size() > 0,
               ExcMessage("No owners found for the point " +
                          std::to_string(pt_n)));
        if (owners_found.size() == 1)
          map_owners_found[pt_n] = *owners_found.begin();
        else
          // Multiple owners
          std::copy(owners_found.begin(),
                    owners_found.end(),
                    std::back_inserter(map_owners_guessed[pt_n]));
      }

    return std::make_tuple(std::move(point_owners),
                           std::move(map_owners_found),
                           std::move(map_owners_guessed));
  }


  template <int dim, int spacedim>
  std::vector<
    std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
  vertex_to_cell_map(const Triangulation<dim, spacedim> &triangulation)
  {
    std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      vertex_to_cell_map(triangulation.n_vertices());
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (const unsigned int i : cell->vertex_indices())
        vertex_to_cell_map[cell->vertex_index(i)].insert(cell);

    // Take care of hanging nodes
    cell = triangulation.begin_active();
    for (; cell != endc; ++cell)
      {
        for (unsigned int i : cell->face_indices())
          {
            if ((cell->at_boundary(i) == false) &&
                (cell->neighbor(i)->is_active()))
              {
                typename Triangulation<dim, spacedim>::active_cell_iterator
                  adjacent_cell = cell->neighbor(i);
                for (unsigned int j = 0; j < cell->face(i)->n_vertices(); ++j)
                  vertex_to_cell_map[cell->face(i)->vertex_index(j)].insert(
                    adjacent_cell);
              }
          }

        // in 3d also loop over the edges
        if (dim == 3)
          {
            for (unsigned int i = 0; i < cell->n_lines(); ++i)
              if (cell->line(i)->has_children())
                // the only place where this vertex could have been
                // hiding is on the mid-edge point of the edge we
                // are looking at
                vertex_to_cell_map[cell->line(i)->child(0)->vertex_index(1)]
                  .insert(cell);
          }
      }

    return vertex_to_cell_map;
  }



  template <int dim, int spacedim>
  std::map<unsigned int, types::global_vertex_index>
  compute_local_to_global_vertex_index_map(
    const parallel::distributed::Triangulation<dim, spacedim> &triangulation)
  {
    std::map<unsigned int, types::global_vertex_index>
      local_to_global_vertex_index;

#ifndef DEAL_II_WITH_MPI

    // without MPI, this function doesn't make sense because on cannot
    // use parallel::distributed::Triangulation in any meaningful
    // way
    (void)triangulation;
    Assert(false,
           ExcMessage("This function does not make any sense "
                      "for parallel::distributed::Triangulation "
                      "objects if you do not have MPI enabled."));

#else

    using active_cell_iterator =
      typename Triangulation<dim, spacedim>::active_cell_iterator;
    const std::vector<std::set<active_cell_iterator>> vertex_to_cell =
      vertex_to_cell_map(triangulation);

    // Create a local index for the locally "owned" vertices
    types::global_vertex_index next_index      = 0;
    unsigned int               max_cellid_size = 0;
    std::set<std::pair<types::subdomain_id, types::global_vertex_index>>
                                                          vertices_added;
    std::map<types::subdomain_id, std::set<unsigned int>> vertices_to_recv;
    std::map<types::subdomain_id,
             std::vector<std::tuple<types::global_vertex_index,
                                    types::global_vertex_index,
                                    std::string>>>
                         vertices_to_send;
    active_cell_iterator cell = triangulation.begin_active(),
                         endc = triangulation.end();
    std::set<active_cell_iterator> missing_vert_cells;
    std::set<unsigned int>         used_vertex_index;
    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            for (const unsigned int i : cell->vertex_indices())
              {
                types::subdomain_id lowest_subdomain_id = cell->subdomain_id();
                typename std::set<active_cell_iterator>::iterator
                  adjacent_cell = vertex_to_cell[cell->vertex_index(i)].begin(),
                  end_adj_cell  = vertex_to_cell[cell->vertex_index(i)].end();
                for (; adjacent_cell != end_adj_cell; ++adjacent_cell)
                  lowest_subdomain_id =
                    std::min(lowest_subdomain_id,
                             (*adjacent_cell)->subdomain_id());

                // See if I "own" this vertex
                if (lowest_subdomain_id == cell->subdomain_id())
                  {
                    // Check that the vertex we are working on a vertex that has
                    // not be dealt with yet
                    if (used_vertex_index.find(cell->vertex_index(i)) ==
                        used_vertex_index.end())
                      {
                        // Set the local index
                        local_to_global_vertex_index[cell->vertex_index(i)] =
                          next_index++;

                        // Store the information that will be sent to the
                        // adjacent cells on other subdomains
                        adjacent_cell =
                          vertex_to_cell[cell->vertex_index(i)].begin();
                        for (; adjacent_cell != end_adj_cell; ++adjacent_cell)
                          if ((*adjacent_cell)->subdomain_id() !=
                              cell->subdomain_id())
                            {
                              std::pair<types::subdomain_id,
                                        types::global_vertex_index>
                                tmp((*adjacent_cell)->subdomain_id(),
                                    cell->vertex_index(i));
                              if (vertices_added.find(tmp) ==
                                  vertices_added.end())
                                {
                                  vertices_to_send[(*adjacent_cell)
                                                     ->subdomain_id()]
                                    .emplace_back(i,
                                                  cell->vertex_index(i),
                                                  cell->id().to_string());
                                  if (cell->id().to_string().size() >
                                      max_cellid_size)
                                    max_cellid_size =
                                      cell->id().to_string().size();
                                  vertices_added.insert(tmp);
                                }
                            }
                        used_vertex_index.insert(cell->vertex_index(i));
                      }
                  }
                else
                  {
                    // We don't own the vertex so we will receive its global
                    // index
                    vertices_to_recv[lowest_subdomain_id].insert(
                      cell->vertex_index(i));
                    missing_vert_cells.insert(cell);
                  }
              }
          }

        // Some hanging nodes are vertices of ghost cells. They need to be
        // received.
        if (cell->is_ghost())
          {
            for (unsigned int i : cell->face_indices())
              {
                if (cell->at_boundary(i) == false)
                  {
                    if (cell->neighbor(i)->is_active())
                      {
                        typename Triangulation<dim,
                                               spacedim>::active_cell_iterator
                          adjacent_cell = cell->neighbor(i);
                        if ((adjacent_cell->is_locally_owned()))
                          {
                            types::subdomain_id adj_subdomain_id =
                              adjacent_cell->subdomain_id();
                            if (cell->subdomain_id() < adj_subdomain_id)
                              for (unsigned int j = 0;
                                   j < cell->face(i)->n_vertices();
                                   ++j)
                                {
                                  vertices_to_recv[cell->subdomain_id()].insert(
                                    cell->face(i)->vertex_index(j));
                                  missing_vert_cells.insert(cell);
                                }
                          }
                      }
                  }
              }
          }
      }

    // Get the size of the largest CellID string
    max_cellid_size =
      Utilities::MPI::max(max_cellid_size, triangulation.get_communicator());

    // Make indices global by getting the number of vertices owned by each
    // processors and shifting the indices accordingly
    types::global_vertex_index shift = 0;
    int                        ierr  = MPI_Exscan(&next_index,
                          &shift,
                          1,
                          DEAL_II_VERTEX_INDEX_MPI_TYPE,
                          MPI_SUM,
                          triangulation.get_communicator());
    AssertThrowMPI(ierr);

    std::map<unsigned int, types::global_vertex_index>::iterator
      global_index_it  = local_to_global_vertex_index.begin(),
      global_index_end = local_to_global_vertex_index.end();
    for (; global_index_it != global_index_end; ++global_index_it)
      global_index_it->second += shift;


    const int mpi_tag = Utilities::MPI::internal::Tags::
      grid_tools_compute_local_to_global_vertex_index_map;
    const int mpi_tag2 = Utilities::MPI::internal::Tags::
      grid_tools_compute_local_to_global_vertex_index_map2;


    // In a first message, send the global ID of the vertices and the local
    // positions in the cells. In a second messages, send the cell ID as a
    // resize string. This is done in two messages so that types are not mixed

    // Send the first message
    std::vector<std::vector<types::global_vertex_index>> vertices_send_buffers(
      vertices_to_send.size());
    std::vector<MPI_Request> first_requests(vertices_to_send.size());
    typename std::map<types::subdomain_id,
                      std::vector<std::tuple<types::global_vertex_index,
                                             types::global_vertex_index,
                                             std::string>>>::iterator
      vert_to_send_it  = vertices_to_send.begin(),
      vert_to_send_end = vertices_to_send.end();
    for (unsigned int i = 0; vert_to_send_it != vert_to_send_end;
         ++vert_to_send_it, ++i)
      {
        int                destination = vert_to_send_it->first;
        const unsigned int n_vertices  = vert_to_send_it->second.size();
        const int          buffer_size = 2 * n_vertices;
        vertices_send_buffers[i].resize(buffer_size);

        // fill the buffer
        for (unsigned int j = 0; j < n_vertices; ++j)
          {
            vertices_send_buffers[i][2 * j] =
              std::get<0>(vert_to_send_it->second[j]);
            vertices_send_buffers[i][2 * j + 1] =
              local_to_global_vertex_index[std::get<1>(
                vert_to_send_it->second[j])];
          }

        // Send the message
        ierr = MPI_Isend(vertices_send_buffers[i].data(),
                         buffer_size,
                         DEAL_II_VERTEX_INDEX_MPI_TYPE,
                         destination,
                         mpi_tag,
                         triangulation.get_communicator(),
                         &first_requests[i]);
        AssertThrowMPI(ierr);
      }

    // Receive the first message
    std::vector<std::vector<types::global_vertex_index>> vertices_recv_buffers(
      vertices_to_recv.size());
    typename std::map<types::subdomain_id, std::set<unsigned int>>::iterator
      vert_to_recv_it  = vertices_to_recv.begin(),
      vert_to_recv_end = vertices_to_recv.end();
    for (unsigned int i = 0; vert_to_recv_it != vert_to_recv_end;
         ++vert_to_recv_it, ++i)
      {
        int                source      = vert_to_recv_it->first;
        const unsigned int n_vertices  = vert_to_recv_it->second.size();
        const int          buffer_size = 2 * n_vertices;
        vertices_recv_buffers[i].resize(buffer_size);

        // Receive the message
        ierr = MPI_Recv(vertices_recv_buffers[i].data(),
                        buffer_size,
                        DEAL_II_VERTEX_INDEX_MPI_TYPE,
                        source,
                        mpi_tag,
                        triangulation.get_communicator(),
                        MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
      }


    // Send second message
    std::vector<std::vector<char>> cellids_send_buffers(
      vertices_to_send.size());
    std::vector<MPI_Request> second_requests(vertices_to_send.size());
    vert_to_send_it = vertices_to_send.begin();
    for (unsigned int i = 0; vert_to_send_it != vert_to_send_end;
         ++vert_to_send_it, ++i)
      {
        int                destination = vert_to_send_it->first;
        const unsigned int n_vertices  = vert_to_send_it->second.size();
        const int          buffer_size = max_cellid_size * n_vertices;
        cellids_send_buffers[i].resize(buffer_size);

        // fill the buffer
        unsigned int pos = 0;
        for (unsigned int j = 0; j < n_vertices; ++j)
          {
            std::string cell_id = std::get<2>(vert_to_send_it->second[j]);
            for (unsigned int k = 0; k < max_cellid_size; ++k, ++pos)
              {
                if (k < cell_id.size())
                  cellids_send_buffers[i][pos] = cell_id[k];
                // if necessary fill up the reserved part of the buffer with an
                // invalid value
                else
                  cellids_send_buffers[i][pos] = '-';
              }
          }

        // Send the message
        ierr = MPI_Isend(cellids_send_buffers[i].data(),
                         buffer_size,
                         MPI_CHAR,
                         destination,
                         mpi_tag2,
                         triangulation.get_communicator(),
                         &second_requests[i]);
        AssertThrowMPI(ierr);
      }

    // Receive the second message
    std::vector<std::vector<char>> cellids_recv_buffers(
      vertices_to_recv.size());
    vert_to_recv_it = vertices_to_recv.begin();
    for (unsigned int i = 0; vert_to_recv_it != vert_to_recv_end;
         ++vert_to_recv_it, ++i)
      {
        int                source      = vert_to_recv_it->first;
        const unsigned int n_vertices  = vert_to_recv_it->second.size();
        const int          buffer_size = max_cellid_size * n_vertices;
        cellids_recv_buffers[i].resize(buffer_size);

        // Receive the message
        ierr = MPI_Recv(cellids_recv_buffers[i].data(),
                        buffer_size,
                        MPI_CHAR,
                        source,
                        mpi_tag2,
                        triangulation.get_communicator(),
                        MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
      }


    // Match the data received with the required vertices
    vert_to_recv_it = vertices_to_recv.begin();
    for (unsigned int i = 0; vert_to_recv_it != vert_to_recv_end;
         ++i, ++vert_to_recv_it)
      {
        for (unsigned int j = 0; j < vert_to_recv_it->second.size(); ++j)
          {
            const unsigned int local_pos_recv = vertices_recv_buffers[i][2 * j];
            const types::global_vertex_index global_id_recv =
              vertices_recv_buffers[i][2 * j + 1];
            const std::string cellid_recv(
              &cellids_recv_buffers[i][max_cellid_size * j],
              &cellids_recv_buffers[i][max_cellid_size * j] + max_cellid_size);
            bool found = false;
            typename std::set<active_cell_iterator>::iterator
              cell_set_it  = missing_vert_cells.begin(),
              end_cell_set = missing_vert_cells.end();
            for (; (found == false) && (cell_set_it != end_cell_set);
                 ++cell_set_it)
              {
                typename std::set<active_cell_iterator>::iterator
                  candidate_cell =
                    vertex_to_cell[(*cell_set_it)->vertex_index(i)].begin(),
                  end_cell =
                    vertex_to_cell[(*cell_set_it)->vertex_index(i)].end();
                for (; candidate_cell != end_cell; ++candidate_cell)
                  {
                    std::string current_cellid =
                      (*candidate_cell)->id().to_string();
                    current_cellid.resize(max_cellid_size, '-');
                    if (current_cellid.compare(cellid_recv) == 0)
                      {
                        local_to_global_vertex_index
                          [(*candidate_cell)->vertex_index(local_pos_recv)] =
                            global_id_recv;
                        found = true;

                        break;
                      }
                  }
              }
          }
      }
#endif

    return local_to_global_vertex_index;
  }



  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern &            cell_connectivity)
  {
    cell_connectivity.reinit(triangulation.n_active_cells(),
                             triangulation.n_active_cells());

    // loop over all cells and their neighbors to build the sparsity
    // pattern. note that it's a bit hard to enter all the connections when a
    // neighbor has children since we would need to find out which of its
    // children is adjacent to the current cell. this problem can be omitted
    // if we only do something if the neighbor has no children -- in that case
    // it is either on the same or a coarser level than we are. in return, we
    // have to add entries in both directions for both cells
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        const unsigned int index = cell->active_cell_index();
        cell_connectivity.add(index, index);
        for (auto f : cell->face_indices())
          if ((cell->at_boundary(f) == false) &&
              (cell->neighbor(f)->has_children() == false))
            {
              const unsigned int other_index =
                cell->neighbor(f)->active_cell_index();
              cell_connectivity.add(index, other_index);
              cell_connectivity.add(other_index, index);
            }
      }
  }



  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern &            cell_connectivity)
  {
    std::vector<std::vector<unsigned int>> vertex_to_cell(
      triangulation.n_vertices());
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        for (const unsigned int v : cell->vertex_indices())
          vertex_to_cell[cell->vertex_index(v)].push_back(
            cell->active_cell_index());
      }

    cell_connectivity.reinit(triangulation.n_active_cells(),
                             triangulation.n_active_cells());
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        for (const unsigned int v : cell->vertex_indices())
          for (unsigned int n = 0;
               n < vertex_to_cell[cell->vertex_index(v)].size();
               ++n)
            cell_connectivity.add(cell->active_cell_index(),
                                  vertex_to_cell[cell->vertex_index(v)][n]);
      }
  }


  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells_on_level(
    const Triangulation<dim, spacedim> &triangulation,
    const unsigned int                  level,
    DynamicSparsityPattern &            cell_connectivity)
  {
    std::vector<std::vector<unsigned int>> vertex_to_cell(
      triangulation.n_vertices());
    for (typename Triangulation<dim, spacedim>::cell_iterator cell =
           triangulation.begin(level);
         cell != triangulation.end(level);
         ++cell)
      {
        for (const unsigned int v : cell->vertex_indices())
          vertex_to_cell[cell->vertex_index(v)].push_back(cell->index());
      }

    cell_connectivity.reinit(triangulation.n_cells(level),
                             triangulation.n_cells(level));
    for (typename Triangulation<dim, spacedim>::cell_iterator cell =
           triangulation.begin(level);
         cell != triangulation.end(level);
         ++cell)
      {
        for (const unsigned int v : cell->vertex_indices())
          for (unsigned int n = 0;
               n < vertex_to_cell[cell->vertex_index(v)].size();
               ++n)
            cell_connectivity.add(cell->index(),
                                  vertex_to_cell[cell->vertex_index(v)][n]);
      }
  }



  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          Triangulation<dim, spacedim> &   triangulation,
                          const SparsityTools::Partitioner partitioner)
  {
    Assert((dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              &triangulation) == nullptr),
           ExcMessage("Objects of type parallel::distributed::Triangulation "
                      "are already partitioned implicitly and can not be "
                      "partitioned again explicitly."));

    std::vector<unsigned int> cell_weights;

    // Get cell weighting if a signal has been attached to the triangulation
    if (!triangulation.signals.cell_weight.empty())
      {
        cell_weights.resize(triangulation.n_active_cells(), 0U);

        // In a first step, obtain the weights of the locally owned
        // cells. For all others, the weight remains at the zero the
        // vector was initialized with above.
        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            cell_weights[cell->active_cell_index()] =
              triangulation.signals.cell_weight(
                cell, Triangulation<dim, spacedim>::CellStatus::CELL_PERSIST);

        // If this is a parallel triangulation, we then need to also
        // get the weights for all other cells. We have asserted above
        // that this function can't be used for
        // parallel::distribute::Triangulation objects, so the only
        // ones we have to worry about here are
        // parallel::shared::Triangulation
        if (const auto shared_tria =
              dynamic_cast<parallel::shared::Triangulation<dim, spacedim> *>(
                &triangulation))
          Utilities::MPI::sum(cell_weights,
                              shared_tria->get_communicator(),
                              cell_weights);
      }

    // Call the other more general function
    partition_triangulation(n_partitions,
                            cell_weights,
                            triangulation,
                            partitioner);
  }



  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          const std::vector<unsigned int> &cell_weights,
                          Triangulation<dim, spacedim> &   triangulation,
                          const SparsityTools::Partitioner partitioner)
  {
    Assert((dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              &triangulation) == nullptr),
           ExcMessage("Objects of type parallel::distributed::Triangulation "
                      "are already partitioned implicitly and can not be "
                      "partitioned again explicitly."));
    Assert(n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));

    // check for an easy return
    if (n_partitions == 1)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_subdomain_id(0);
        return;
      }

    // we decompose the domain by first
    // generating the connection graph of all
    // cells with their neighbors, and then
    // passing this graph off to METIS.
    // finally defer to the other function for
    // partitioning and assigning subdomain ids
    DynamicSparsityPattern cell_connectivity;
    get_face_connectivity_of_cells(triangulation, cell_connectivity);

    SparsityPattern sp_cell_connectivity;
    sp_cell_connectivity.copy_from(cell_connectivity);
    partition_triangulation(n_partitions,
                            cell_weights,
                            sp_cell_connectivity,
                            triangulation,
                            partitioner);
  }



  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int            n_partitions,
                          const SparsityPattern &       cell_connection_graph,
                          Triangulation<dim, spacedim> &triangulation,
                          const SparsityTools::Partitioner partitioner)
  {
    Assert((dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              &triangulation) == nullptr),
           ExcMessage("Objects of type parallel::distributed::Triangulation "
                      "are already partitioned implicitly and can not be "
                      "partitioned again explicitly."));

    std::vector<unsigned int> cell_weights;

    // Get cell weighting if a signal has been attached to the triangulation
    if (!triangulation.signals.cell_weight.empty())
      {
        cell_weights.resize(triangulation.n_active_cells(), 0U);

        // In a first step, obtain the weights of the locally owned
        // cells. For all others, the weight remains at the zero the
        // vector was initialized with above.
        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            cell_weights[cell->active_cell_index()] =
              triangulation.signals.cell_weight(
                cell, Triangulation<dim, spacedim>::CellStatus::CELL_PERSIST);

        // If this is a parallel triangulation, we then need to also
        // get the weights for all other cells. We have asserted above
        // that this function can't be used for
        // parallel::distribute::Triangulation objects, so the only
        // ones we have to worry about here are
        // parallel::shared::Triangulation
        if (const auto shared_tria =
              dynamic_cast<parallel::shared::Triangulation<dim, spacedim> *>(
                &triangulation))
          Utilities::MPI::sum(cell_weights,
                              shared_tria->get_communicator(),
                              cell_weights);
      }

    // Call the other more general function
    partition_triangulation(n_partitions,
                            cell_weights,
                            cell_connection_graph,
                            triangulation,
                            partitioner);
  }



  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          const std::vector<unsigned int> &cell_weights,
                          const SparsityPattern &       cell_connection_graph,
                          Triangulation<dim, spacedim> &triangulation,
                          const SparsityTools::Partitioner partitioner)
  {
    Assert((dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              &triangulation) == nullptr),
           ExcMessage("Objects of type parallel::distributed::Triangulation "
                      "are already partitioned implicitly and can not be "
                      "partitioned again explicitly."));
    Assert(n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));
    Assert(cell_connection_graph.n_rows() == triangulation.n_active_cells(),
           ExcMessage("Connectivity graph has wrong size"));
    Assert(cell_connection_graph.n_cols() == triangulation.n_active_cells(),
           ExcMessage("Connectivity graph has wrong size"));

    // signal that partitioning is going to happen
    triangulation.signals.pre_partition();

    // check for an easy return
    if (n_partitions == 1)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_subdomain_id(0);
        return;
      }

    // partition this connection graph and get
    // back a vector of indices, one per degree
    // of freedom (which is associated with a
    // cell)
    std::vector<unsigned int> partition_indices(triangulation.n_active_cells());
    SparsityTools::partition(cell_connection_graph,
                             cell_weights,
                             n_partitions,
                             partition_indices,
                             partitioner);

    // finally loop over all cells and set the subdomain ids
    for (const auto &cell : triangulation.active_cell_iterators())
      cell->set_subdomain_id(partition_indices[cell->active_cell_index()]);
  }


  namespace internal
  {
    /**
     * recursive helper function for partition_triangulation_zorder
     */
    template <class IT>
    void
    set_subdomain_id_in_zorder_recursively(IT                 cell,
                                           unsigned int &     current_proc_idx,
                                           unsigned int &     current_cell_idx,
                                           const unsigned int n_active_cells,
                                           const unsigned int n_partitions)
    {
      if (cell->is_active())
        {
          while (current_cell_idx >=
                 std::floor(static_cast<uint_least64_t>(n_active_cells) *
                            (current_proc_idx + 1) / n_partitions))
            ++current_proc_idx;
          cell->set_subdomain_id(current_proc_idx);
          ++current_cell_idx;
        }
      else
        {
          for (unsigned int n = 0; n < cell->n_children(); ++n)
            set_subdomain_id_in_zorder_recursively(cell->child(n),
                                                   current_proc_idx,
                                                   current_cell_idx,
                                                   n_active_cells,
                                                   n_partitions);
        }
    }
  } // namespace internal

  template <int dim, int spacedim>
  void
  partition_triangulation_zorder(const unsigned int            n_partitions,
                                 Triangulation<dim, spacedim> &triangulation,
                                 const bool                    group_siblings)
  {
    Assert((dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              &triangulation) == nullptr),
           ExcMessage("Objects of type parallel::distributed::Triangulation "
                      "are already partitioned implicitly and can not be "
                      "partitioned again explicitly."));
    Assert(n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));

    // signal that partitioning is going to happen
    triangulation.signals.pre_partition();

    // check for an easy return
    if (n_partitions == 1)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_subdomain_id(0);
        return;
      }

    // Duplicate the coarse cell reordoring
    // as done in p4est
    std::vector<types::global_dof_index> coarse_cell_to_p4est_tree_permutation;
    std::vector<types::global_dof_index> p4est_tree_to_coarse_cell_permutation;

    DynamicSparsityPattern cell_connectivity;
    GridTools::get_vertex_connectivity_of_cells_on_level(triangulation,
                                                         0,
                                                         cell_connectivity);
    coarse_cell_to_p4est_tree_permutation.resize(triangulation.n_cells(0));
    SparsityTools::reorder_hierarchical(cell_connectivity,
                                        coarse_cell_to_p4est_tree_permutation);

    p4est_tree_to_coarse_cell_permutation =
      Utilities::invert_permutation(coarse_cell_to_p4est_tree_permutation);

    unsigned int       current_proc_idx = 0;
    unsigned int       current_cell_idx = 0;
    const unsigned int n_active_cells   = triangulation.n_active_cells();

    // set subdomain id for active cell descendants
    // of each coarse cell in permuted order
    for (unsigned int idx = 0; idx < triangulation.n_cells(0); ++idx)
      {
        const unsigned int coarse_cell_idx =
          p4est_tree_to_coarse_cell_permutation[idx];
        typename Triangulation<dim, spacedim>::cell_iterator coarse_cell(
          &triangulation, 0, coarse_cell_idx);

        internal::set_subdomain_id_in_zorder_recursively(coarse_cell,
                                                         current_proc_idx,
                                                         current_cell_idx,
                                                         n_active_cells,
                                                         n_partitions);
      }

    // if all children of a cell are active (e.g. we
    // have a cell that is refined once and no part
    // is refined further), p4est places all of them
    // on the same processor. The new owner will be
    // the processor with the largest number of children
    // (ties are broken by picking the lower rank).
    // Duplicate this logic here.
    if (group_siblings)
      {
        typename Triangulation<dim, spacedim>::cell_iterator
          cell = triangulation.begin(),
          endc = triangulation.end();
        for (; cell != endc; ++cell)
          {
            if (cell->is_active())
              continue;
            bool                                 all_children_active = true;
            std::map<unsigned int, unsigned int> map_cpu_n_cells;
            for (unsigned int n = 0; n < cell->n_children(); ++n)
              if (!cell->child(n)->is_active())
                {
                  all_children_active = false;
                  break;
                }
              else
                ++map_cpu_n_cells[cell->child(n)->subdomain_id()];

            if (!all_children_active)
              continue;

            unsigned int new_owner = cell->child(0)->subdomain_id();
            for (std::map<unsigned int, unsigned int>::iterator it =
                   map_cpu_n_cells.begin();
                 it != map_cpu_n_cells.end();
                 ++it)
              if (it->second > map_cpu_n_cells[new_owner])
                new_owner = it->first;

            for (unsigned int n = 0; n < cell->n_children(); ++n)
              cell->child(n)->set_subdomain_id(new_owner);
          }
      }
  }


  template <int dim, int spacedim>
  void
  partition_multigrid_levels(Triangulation<dim, spacedim> &triangulation)
  {
    unsigned int n_levels = triangulation.n_levels();
    for (int lvl = n_levels - 1; lvl >= 0; --lvl)
      {
        typename Triangulation<dim, spacedim>::cell_iterator
          cell = triangulation.begin(lvl),
          endc = triangulation.end(lvl);
        for (; cell != endc; ++cell)
          {
            if (cell->is_active())
              cell->set_level_subdomain_id(cell->subdomain_id());
            else
              {
                Assert(cell->child(0)->level_subdomain_id() !=
                         numbers::artificial_subdomain_id,
                       ExcInternalError());
                cell->set_level_subdomain_id(
                  cell->child(0)->level_subdomain_id());
              }
          }
      }
  }



  template <int dim, int spacedim>
  std::vector<types::subdomain_id>
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            const std::vector<CellId> &         cell_ids)
  {
    std::vector<types::subdomain_id> subdomain_ids;
    subdomain_ids.reserve(cell_ids.size());

    if (dynamic_cast<
          const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
          &triangulation) != nullptr)
      {
        Assert(false, ExcNotImplemented());
      }
    else if (const parallel::distributed::Triangulation<dim, spacedim>
               *parallel_tria = dynamic_cast<
                 const parallel::distributed::Triangulation<dim, spacedim> *>(
                 &triangulation))
      {
#ifndef DEAL_II_WITH_P4EST
        Assert(
          false,
          ExcMessage(
            "You are attempting to use a functionality that is only available "
            "if deal.II was configured to use p4est, but cmake did not find a "
            "valid p4est library."));
#else
        // for parallel distributed triangulations, we will ask the p4est oracle
        // about the global partitioning of active cells since this information
        // is stored on every process
        for (const auto &cell_id : cell_ids)
          {
            // find descendent from coarse quadrant
            typename dealii::internal::p4est::types<dim>::quadrant p4est_cell,
              p4est_children[GeometryInfo<dim>::max_children_per_cell];

            dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_cell);
            for (const auto &child_index : cell_id.get_child_indices())
              {
                dealii::internal::p4est::init_quadrant_children<dim>(
                  p4est_cell, p4est_children);
                p4est_cell =
                  p4est_children[static_cast<unsigned int>(child_index)];
              }

            // find owning process, i.e., the subdomain id
            const int owner =
              dealii::internal::p4est::functions<dim>::comm_find_owner(
                const_cast<typename dealii::internal::p4est::types<dim>::forest
                             *>(parallel_tria->get_p4est()),
                cell_id.get_coarse_cell_id(),
                &p4est_cell,
                Utilities::MPI::this_mpi_process(
                  parallel_tria->get_communicator()));

            Assert(owner >= 0, ExcMessage("p4est should know the owner."));

            subdomain_ids.push_back(owner);
          }
#endif
      }
    else if (const parallel::shared::Triangulation<dim, spacedim> *shared_tria =
               dynamic_cast<const parallel::shared::Triangulation<dim, spacedim>
                              *>(&triangulation))
      {
        // for parallel shared triangulations, we need to access true subdomain
        // ids which are also valid for artificial cells
        const std::vector<types::subdomain_id> &true_subdomain_ids_of_cells =
          shared_tria->get_true_subdomain_ids_of_cells();

        for (const auto &cell_id : cell_ids)
          {
            const unsigned int active_cell_index =
              shared_tria->create_cell_iterator(cell_id)->active_cell_index();
            subdomain_ids.push_back(
              true_subdomain_ids_of_cells[active_cell_index]);
          }
      }
    else
      {
        // the most general type of triangulation is the serial one. here, all
        // subdomain information is directly available
        for (const auto &cell_id : cell_ids)
          {
            subdomain_ids.push_back(
              triangulation.create_cell_iterator(cell_id)->subdomain_id());
          }
      }

    return subdomain_ids;
  }



  template <int dim, int spacedim>
  void
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            std::vector<types::subdomain_id> &  subdomain)
  {
    Assert(subdomain.size() == triangulation.n_active_cells(),
           ExcDimensionMismatch(subdomain.size(),
                                triangulation.n_active_cells()));
    for (const auto &cell : triangulation.active_cell_iterators())
      subdomain[cell->active_cell_index()] = cell->subdomain_id();
  }



  template <int dim, int spacedim>
  unsigned int
  count_cells_with_subdomain_association(
    const Triangulation<dim, spacedim> &triangulation,
    const types::subdomain_id           subdomain)
  {
    unsigned int count = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->subdomain_id() == subdomain)
        ++count;

    return count;
  }



  template <int dim, int spacedim>
  std::vector<bool>
  get_locally_owned_vertices(const Triangulation<dim, spacedim> &triangulation)
  {
    // start with all vertices
    std::vector<bool> locally_owned_vertices =
      triangulation.get_used_vertices();

    // if the triangulation is distributed, eliminate those that
    // are owned by other processors -- either because the vertex is
    // on an artificial cell, or because it is on a ghost cell with
    // a smaller subdomain
    if (const auto *tr = dynamic_cast<
          const parallel::DistributedTriangulationBase<dim, spacedim> *>(
          &triangulation))
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->is_artificial() ||
            (cell->is_ghost() &&
             (cell->subdomain_id() < tr->locally_owned_subdomain())))
          for (const unsigned int v : cell->vertex_indices())
            locally_owned_vertices[cell->vertex_index(v)] = false;

    return locally_owned_vertices;
  }



  template <int dim, int spacedim>
  double
  minimal_cell_diameter(const Triangulation<dim, spacedim> &triangulation,
                        const Mapping<dim, spacedim> &      mapping)
  {
    double min_diameter = std::numeric_limits<double>::max();
    for (const auto &cell : triangulation.active_cell_iterators())
      if (!cell->is_artificial())
        min_diameter = std::min(min_diameter, cell->diameter(mapping));

    double global_min_diameter = 0;

#ifdef DEAL_II_WITH_MPI
    if (const parallel::TriangulationBase<dim, spacedim> *p_tria =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &triangulation))
      global_min_diameter =
        Utilities::MPI::min(min_diameter, p_tria->get_communicator());
    else
#endif
      global_min_diameter = min_diameter;

    return global_min_diameter;
  }



  template <int dim, int spacedim>
  double
  maximal_cell_diameter(const Triangulation<dim, spacedim> &triangulation,
                        const Mapping<dim, spacedim> &      mapping)
  {
    double max_diameter = 0.;
    for (const auto &cell : triangulation.active_cell_iterators())
      if (!cell->is_artificial())
        max_diameter = std::max(max_diameter, cell->diameter(mapping));

    double global_max_diameter = 0;

#ifdef DEAL_II_WITH_MPI
    if (const parallel::TriangulationBase<dim, spacedim> *p_tria =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &triangulation))
      global_max_diameter =
        Utilities::MPI::max(max_diameter, p_tria->get_communicator());
    else
#endif
      global_max_diameter = max_diameter;

    return global_max_diameter;
  }



  namespace internal
  {
    namespace FixUpDistortedChildCells
    {
      // compute the mean square
      // deviation of the alternating
      // forms of the children of the
      // given object from that of
      // the object itself. for
      // objects with
      // structdim==spacedim, the
      // alternating form is the
      // determinant of the jacobian,
      // whereas for faces with
      // structdim==spacedim-1, the
      // alternating form is the
      // (signed and scaled) normal
      // vector
      //
      // this average square
      // deviation is computed for an
      // object where the center node
      // has been replaced by the
      // second argument to this
      // function
      template <typename Iterator, int spacedim>
      double
      objective_function(const Iterator &       object,
                         const Point<spacedim> &object_mid_point)
      {
        const unsigned int structdim =
          Iterator::AccessorType::structure_dimension;
        Assert(spacedim == Iterator::AccessorType::dimension,
               ExcInternalError());

        // everything below is wrong
        // if not for the following
        // condition
        Assert(object->refinement_case() ==
                 RefinementCase<structdim>::isotropic_refinement,
               ExcNotImplemented());
        // first calculate the
        // average alternating form
        // for the parent cell/face
        Point<spacedim>
          parent_vertices[GeometryInfo<structdim>::vertices_per_cell];
        Tensor<spacedim - structdim, spacedim>
          parent_alternating_forms[GeometryInfo<structdim>::vertices_per_cell];

        for (const unsigned int i : object->vertex_indices())
          parent_vertices[i] = object->vertex(i);

        GeometryInfo<structdim>::alternating_form_at_vertices(
          parent_vertices, parent_alternating_forms);

        const Tensor<spacedim - structdim, spacedim>
          average_parent_alternating_form =
            std::accumulate(parent_alternating_forms,
                            parent_alternating_forms +
                              GeometryInfo<structdim>::vertices_per_cell,
                            Tensor<spacedim - structdim, spacedim>());

        // now do the same
        // computation for the
        // children where we use the
        // given location for the
        // object mid point instead of
        // the one the triangulation
        // currently reports
        Point<spacedim>
          child_vertices[GeometryInfo<structdim>::max_children_per_cell]
                        [GeometryInfo<structdim>::vertices_per_cell];
        Tensor<spacedim - structdim, spacedim> child_alternating_forms
          [GeometryInfo<structdim>::max_children_per_cell]
          [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int c = 0; c < object->n_children(); ++c)
          for (const unsigned int i : object->child(c)->vertex_indices())
            child_vertices[c][i] = object->child(c)->vertex(i);

        // replace mid-object
        // vertex. note that for
        // child i, the mid-object
        // vertex happens to have the
        // number
        // max_children_per_cell-i
        for (unsigned int c = 0; c < object->n_children(); ++c)
          child_vertices[c][GeometryInfo<structdim>::max_children_per_cell - c -
                            1] = object_mid_point;

        for (unsigned int c = 0; c < object->n_children(); ++c)
          GeometryInfo<structdim>::alternating_form_at_vertices(
            child_vertices[c], child_alternating_forms[c]);

        // on a uniformly refined
        // hypercube object, the child
        // alternating forms should
        // all be smaller by a factor
        // of 2^structdim than the
        // ones of the parent. as a
        // consequence, we'll use the
        // squared deviation from
        // this ideal value as an
        // objective function
        double objective = 0;
        for (unsigned int c = 0; c < object->n_children(); ++c)
          for (const unsigned int i : object->child(c)->vertex_indices())
            objective +=
              (child_alternating_forms[c][i] -
               average_parent_alternating_form / std::pow(2., 1. * structdim))
                .norm_square();

        return objective;
      }


      /**
       * Return the location of the midpoint
       * of the 'f'th face (vertex) of this 1d
       * object.
       */
      template <typename Iterator>
      Point<Iterator::AccessorType::space_dimension>
      get_face_midpoint(const Iterator &   object,
                        const unsigned int f,
                        std::integral_constant<int, 1>)
      {
        return object->vertex(f);
      }



      /**
       * Return the location of the midpoint
       * of the 'f'th face (line) of this 2d
       * object.
       */
      template <typename Iterator>
      Point<Iterator::AccessorType::space_dimension>
      get_face_midpoint(const Iterator &   object,
                        const unsigned int f,
                        std::integral_constant<int, 2>)
      {
        return object->line(f)->center();
      }



      /**
       * Return the location of the midpoint
       * of the 'f'th face (quad) of this 3d
       * object.
       */
      template <typename Iterator>
      Point<Iterator::AccessorType::space_dimension>
      get_face_midpoint(const Iterator &   object,
                        const unsigned int f,
                        std::integral_constant<int, 3>)
      {
        return object->face(f)->center();
      }



      /**
       * Compute the minimal diameter of an
       * object by looking for the minimal
       * distance between the mid-points of
       * its faces. This minimal diameter is
       * used to determine the step length
       * for our grid cell improvement
       * algorithm, and it should be small
       * enough that the point moves around
       * within the cell even if it is highly
       * elongated -- thus, the diameter of
       * the object is not a good measure,
       * while the minimal diameter is. Note
       * that the algorithm below works for
       * both cells that are long rectangles
       * with parallel sides where the
       * nearest distance is between opposite
       * edges as well as highly slanted
       * parallelograms where the shortest
       * distance is between neighboring
       * edges.
       */
      template <typename Iterator>
      double
      minimal_diameter(const Iterator &object)
      {
        const unsigned int structdim =
          Iterator::AccessorType::structure_dimension;

        double diameter = object->diameter();
        for (const unsigned int f : object->face_indices())
          for (unsigned int e = f + 1; e < object->n_faces(); ++e)
            diameter = std::min(
              diameter,
              get_face_midpoint(object,
                                f,
                                std::integral_constant<int, structdim>())
                .distance(get_face_midpoint(
                  object, e, std::integral_constant<int, structdim>())));

        return diameter;
      }



      /**
       * Try to fix up a single cell by moving around its midpoint. Return
       * whether we succeeded with this.
       */
      template <typename Iterator>
      bool
      fix_up_object(const Iterator &object)
      {
        const unsigned int structdim =
          Iterator::AccessorType::structure_dimension;
        const unsigned int spacedim = Iterator::AccessorType::space_dimension;

        // right now we can only deal with cells that have been refined
        // isotropically because that is the only case where we have a cell
        // mid-point that can be moved around without having to consider
        // boundary information
        Assert(object->has_children(), ExcInternalError());
        Assert(object->refinement_case() ==
                 RefinementCase<structdim>::isotropic_refinement,
               ExcNotImplemented());

        // get the current location of the object mid-vertex:
        Point<spacedim> object_mid_point = object->child(0)->vertex(
          GeometryInfo<structdim>::max_children_per_cell - 1);

        // now do a few steepest descent steps to reduce the objective
        // function. compute the diameter in the helper function above
        unsigned int iteration = 0;
        const double diameter  = minimal_diameter(object);

        // current value of objective function and initial delta
        double current_value = objective_function(object, object_mid_point);
        double initial_delta = 0;

        do
          {
            // choose a step length that is initially 1/4 of the child
            // objects' diameter, and a sequence whose sum does not converge
            // (to avoid premature termination of the iteration)
            const double step_length = diameter / 4 / (iteration + 1);

            // compute the objective function's derivative using a two-sided
            // difference formula with eps=step_length/10
            Tensor<1, spacedim> gradient;
            for (unsigned int d = 0; d < spacedim; ++d)
              {
                const double eps = step_length / 10;

                Tensor<1, spacedim> h;
                h[d] = eps / 2;

                gradient[d] =
                  (objective_function(
                     object, project_to_object(object, object_mid_point + h)) -
                   objective_function(
                     object, project_to_object(object, object_mid_point - h))) /
                  eps;
              }

            // there is nowhere to go
            if (gradient.norm() == 0)
              break;

            // We need to go in direction -gradient. the optimal value of the
            // objective function is zero, so assuming that the model is
            // quadratic we would have to go -2*val/||gradient|| in this
            // direction, make sure we go at most step_length into this
            // direction
            object_mid_point -=
              std::min(2 * current_value / (gradient * gradient),
                       step_length / gradient.norm()) *
              gradient;
            object_mid_point = project_to_object(object, object_mid_point);

            // compute current value of the objective function
            const double previous_value = current_value;
            current_value = objective_function(object, object_mid_point);

            if (iteration == 0)
              initial_delta = (previous_value - current_value);

            // stop if we aren't moving much any more
            if ((iteration >= 1) &&
                ((previous_value - current_value < 0) ||
                 (std::fabs(previous_value - current_value) <
                  0.001 * initial_delta)))
              break;

            ++iteration;
          }
        while (iteration < 20);

        // verify that the new
        // location is indeed better
        // than the one before. check
        // this by comparing whether
        // the minimum value of the
        // products of parent and
        // child alternating forms is
        // positive. for cells this
        // means that the
        // determinants have the same
        // sign, for faces that the
        // face normals of parent and
        // children point in the same
        // general direction
        double old_min_product, new_min_product;

        Point<spacedim>
          parent_vertices[GeometryInfo<structdim>::vertices_per_cell];
        for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
          parent_vertices[i] = object->vertex(i);

        Tensor<spacedim - structdim, spacedim>
          parent_alternating_forms[GeometryInfo<structdim>::vertices_per_cell];
        GeometryInfo<structdim>::alternating_form_at_vertices(
          parent_vertices, parent_alternating_forms);

        Point<spacedim>
          child_vertices[GeometryInfo<structdim>::max_children_per_cell]
                        [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int c = 0; c < object->n_children(); ++c)
          for (const unsigned int i : object->child(c)->vertex_indices())
            child_vertices[c][i] = object->child(c)->vertex(i);

        Tensor<spacedim - structdim, spacedim> child_alternating_forms
          [GeometryInfo<structdim>::max_children_per_cell]
          [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int c = 0; c < object->n_children(); ++c)
          GeometryInfo<structdim>::alternating_form_at_vertices(
            child_vertices[c], child_alternating_forms[c]);

        old_min_product =
          child_alternating_forms[0][0] * parent_alternating_forms[0];
        for (unsigned int c = 0; c < object->n_children(); ++c)
          for (const unsigned int i : object->child(c)->vertex_indices())
            for (const unsigned int j : object->vertex_indices())
              old_min_product = std::min<double>(old_min_product,
                                                 child_alternating_forms[c][i] *
                                                   parent_alternating_forms[j]);

        // for the new minimum value,
        // replace mid-object
        // vertex. note that for child
        // i, the mid-object vertex
        // happens to have the number
        // max_children_per_cell-i
        for (unsigned int c = 0; c < object->n_children(); ++c)
          child_vertices[c][GeometryInfo<structdim>::max_children_per_cell - c -
                            1] = object_mid_point;

        for (unsigned int c = 0; c < object->n_children(); ++c)
          GeometryInfo<structdim>::alternating_form_at_vertices(
            child_vertices[c], child_alternating_forms[c]);

        new_min_product =
          child_alternating_forms[0][0] * parent_alternating_forms[0];
        for (unsigned int c = 0; c < object->n_children(); ++c)
          for (const unsigned int i : object->child(c)->vertex_indices())
            for (const unsigned int j : object->vertex_indices())
              new_min_product = std::min<double>(new_min_product,
                                                 child_alternating_forms[c][i] *
                                                   parent_alternating_forms[j]);

        // if new minimum value is
        // better than before, then set the
        // new mid point. otherwise
        // return this object as one of
        // those that can't apparently
        // be fixed
        if (new_min_product >= old_min_product)
          object->child(0)->vertex(
            GeometryInfo<structdim>::max_children_per_cell - 1) =
            object_mid_point;

        // return whether after this
        // operation we have an object that
        // is well oriented
        return (std::max(new_min_product, old_min_product) > 0);
      }



      // possibly fix up the faces of a cell by moving around its mid-points
      template <int dim, int spacedim>
      void
      fix_up_faces(
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &cell,
        std::integral_constant<int, dim>,
        std::integral_constant<int, spacedim>)
      {
        // see if we first can fix up some of the faces of this object. We can
        // mess with faces if and only if the neighboring cell is not even
        // more refined than we are (since in that case the sub-faces have
        // themselves children that we can't move around any more). however,
        // the latter case shouldn't happen anyway: if the current face is
        // distorted but the neighbor is even more refined, then the face had
        // been deformed before already, and had been ignored at the time; we
        // should then also be able to ignore it this time as well
        for (auto f : cell->face_indices())
          {
            Assert(cell->face(f)->has_children(), ExcInternalError());
            Assert(cell->face(f)->refinement_case() ==
                     RefinementCase<dim - 1>::isotropic_refinement,
                   ExcInternalError());

            bool subface_is_more_refined = false;
            for (unsigned int g = 0;
                 g < GeometryInfo<dim>::max_children_per_face;
                 ++g)
              if (cell->face(f)->child(g)->has_children())
                {
                  subface_is_more_refined = true;
                  break;
                }

            if (subface_is_more_refined == true)
              continue;

            // we finally know that we can do something about this face
            fix_up_object(cell->face(f));
          }
      }
    } /* namespace FixUpDistortedChildCells */
  }   /* namespace internal */


  template <int dim, int spacedim>
  typename Triangulation<dim, spacedim>::DistortedCellList
  fix_up_distorted_child_cells(
    const typename Triangulation<dim, spacedim>::DistortedCellList
      &distorted_cells,
    Triangulation<dim, spacedim> & /*triangulation*/)
  {
    static_assert(
      dim != 1 && spacedim != 1,
      "This function is only valid when dim != 1 or spacedim != 1.");
    typename Triangulation<dim, spacedim>::DistortedCellList unfixable_subset;

    // loop over all cells that we have to fix up
    for (typename std::list<
           typename Triangulation<dim, spacedim>::cell_iterator>::const_iterator
           cell_ptr = distorted_cells.distorted_cells.begin();
         cell_ptr != distorted_cells.distorted_cells.end();
         ++cell_ptr)
      {
        const typename Triangulation<dim, spacedim>::cell_iterator cell =
          *cell_ptr;

        Assert(!cell->is_active(),
               ExcMessage(
                 "This function is only valid for a list of cells that "
                 "have children (i.e., no cell in the list may be active)."));

        internal::FixUpDistortedChildCells::fix_up_faces(
          cell,
          std::integral_constant<int, dim>(),
          std::integral_constant<int, spacedim>());

        // If possible, fix up the object.
        if (!internal::FixUpDistortedChildCells::fix_up_object(cell))
          unfixable_subset.distorted_cells.push_back(cell);
      }

    return unfixable_subset;
  }



  template <int dim, int spacedim>
  void
  copy_boundary_to_manifold_id(Triangulation<dim, spacedim> &tria,
                               const bool                    reset_boundary_ids)
  {
    const auto                      src_boundary_ids = tria.get_boundary_ids();
    std::vector<types::manifold_id> dst_manifold_ids(src_boundary_ids.size());
    auto                            m_it = dst_manifold_ids.begin();
    for (const auto b : src_boundary_ids)
      {
        *m_it = static_cast<types::manifold_id>(b);
        ++m_it;
      }
    const std::vector<types::boundary_id> reset_boundary_id =
      reset_boundary_ids ?
        std::vector<types::boundary_id>(src_boundary_ids.size(), 0) :
        src_boundary_ids;
    map_boundary_to_manifold_ids(src_boundary_ids,
                                 dst_manifold_ids,
                                 tria,
                                 reset_boundary_id);
  }



  template <int dim, int spacedim>
  void
  map_boundary_to_manifold_ids(
    const std::vector<types::boundary_id> &src_boundary_ids,
    const std::vector<types::manifold_id> &dst_manifold_ids,
    Triangulation<dim, spacedim> &         tria,
    const std::vector<types::boundary_id> &reset_boundary_ids_)
  {
    AssertDimension(src_boundary_ids.size(), dst_manifold_ids.size());
    const auto reset_boundary_ids =
      reset_boundary_ids_.size() ? reset_boundary_ids_ : src_boundary_ids;
    AssertDimension(reset_boundary_ids.size(), src_boundary_ids.size());

    // in 3d, we not only have to copy boundary ids of faces, but also of edges
    // because we see them twice (once from each adjacent boundary face),
    // we cannot immediately reset their boundary ids. thus, copy first
    // and reset later
    if (dim >= 3)
      for (const auto &cell : tria.active_cell_iterators())
        for (auto f : cell->face_indices())
          if (cell->face(f)->at_boundary())
            for (unsigned int e = 0; e < cell->face(f)->n_lines(); ++e)
              {
                const auto         bid = cell->face(f)->line(e)->boundary_id();
                const unsigned int ind = std::find(src_boundary_ids.begin(),
                                                   src_boundary_ids.end(),
                                                   bid) -
                                         src_boundary_ids.begin();
                if (ind < src_boundary_ids.size())
                  cell->face(f)->line(e)->set_manifold_id(
                    dst_manifold_ids[ind]);
              }

    // now do cells
    for (const auto &cell : tria.active_cell_iterators())
      for (auto f : cell->face_indices())
        if (cell->face(f)->at_boundary())
          {
            const auto         bid = cell->face(f)->boundary_id();
            const unsigned int ind =
              std::find(src_boundary_ids.begin(), src_boundary_ids.end(), bid) -
              src_boundary_ids.begin();

            if (ind < src_boundary_ids.size())
              {
                // assign the manifold id
                cell->face(f)->set_manifold_id(dst_manifold_ids[ind]);
                // then reset boundary id
                cell->face(f)->set_boundary_id(reset_boundary_ids[ind]);
              }

            if (dim >= 3)
              for (unsigned int e = 0; e < cell->face(f)->n_lines(); ++e)
                {
                  const auto bid = cell->face(f)->line(e)->boundary_id();
                  const unsigned int ind = std::find(src_boundary_ids.begin(),
                                                     src_boundary_ids.end(),
                                                     bid) -
                                           src_boundary_ids.begin();
                  if (ind < src_boundary_ids.size())
                    cell->face(f)->line(e)->set_boundary_id(
                      reset_boundary_ids[ind]);
                }
          }
  }


  template <int dim, int spacedim>
  void
  copy_material_to_manifold_id(Triangulation<dim, spacedim> &tria,
                               const bool                    compute_face_ids)
  {
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = tria.begin_active(),
      endc = tria.end();

    for (; cell != endc; ++cell)
      {
        cell->set_manifold_id(cell->material_id());
        if (compute_face_ids == true)
          {
            for (auto f : cell->face_indices())
              {
                if (cell->at_boundary(f) == false)
                  cell->face(f)->set_manifold_id(
                    std::min(cell->material_id(),
                             cell->neighbor(f)->material_id()));
                else
                  cell->face(f)->set_manifold_id(cell->material_id());
              }
          }
      }
  }


  template <int dim, int spacedim>
  void
  assign_co_dimensional_manifold_indicators(
    Triangulation<dim, spacedim> &            tria,
    const std::function<types::manifold_id(
      const std::set<types::manifold_id> &)> &disambiguation_function,
    bool                                      overwrite_only_flat_manifold_ids)
  {
    // Easy case first:
    if (dim == 1)
      return;
    const unsigned int n_subobjects =
      dim == 2 ? tria.n_lines() : tria.n_lines() + tria.n_quads();

    // If user index is zero, then it has not been set.
    std::vector<std::set<types::manifold_id>> manifold_ids(n_subobjects + 1);
    std::vector<unsigned int>                 backup;
    tria.save_user_indices(backup);
    tria.clear_user_data();

    unsigned next_index = 1;
    for (auto &cell : tria.active_cell_iterators())
      {
        if (dim > 1)
          for (unsigned int l = 0; l < cell->n_lines(); ++l)
            {
              if (cell->line(l)->user_index() == 0)
                {
                  AssertIndexRange(next_index, n_subobjects + 1);
                  manifold_ids[next_index].insert(cell->manifold_id());
                  cell->line(l)->set_user_index(next_index++);
                }
              else
                manifold_ids[cell->line(l)->user_index()].insert(
                  cell->manifold_id());
            }
        if (dim > 2)
          for (unsigned int l = 0; l < cell->n_faces(); ++l)
            {
              if (cell->quad(l)->user_index() == 0)
                {
                  AssertIndexRange(next_index, n_subobjects + 1);
                  manifold_ids[next_index].insert(cell->manifold_id());
                  cell->quad(l)->set_user_index(next_index++);
                }
              else
                manifold_ids[cell->quad(l)->user_index()].insert(
                  cell->manifold_id());
            }
      }
    for (auto &cell : tria.active_cell_iterators())
      {
        if (dim > 1)
          for (unsigned int l = 0; l < cell->n_lines(); ++l)
            {
              const auto id = cell->line(l)->user_index();
              // Make sure we change the manifold indicator only once
              if (id != 0)
                {
                  if (cell->line(l)->manifold_id() ==
                        numbers::flat_manifold_id ||
                      overwrite_only_flat_manifold_ids == false)
                    cell->line(l)->set_manifold_id(
                      disambiguation_function(manifold_ids[id]));
                  cell->line(l)->set_user_index(0);
                }
            }
        if (dim > 2)
          for (unsigned int l = 0; l < cell->n_faces(); ++l)
            {
              const auto id = cell->quad(l)->user_index();
              // Make sure we change the manifold indicator only once
              if (id != 0)
                {
                  if (cell->quad(l)->manifold_id() ==
                        numbers::flat_manifold_id ||
                      overwrite_only_flat_manifold_ids == false)
                    cell->quad(l)->set_manifold_id(
                      disambiguation_function(manifold_ids[id]));
                  cell->quad(l)->set_user_index(0);
                }
            }
      }
    tria.load_user_indices(backup);
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


  template <int dim, int spacedim>
  void
  remove_hanging_nodes(Triangulation<dim, spacedim> &tria,
                       const bool                    isotropic,
                       const unsigned int            max_iterations)
  {
    unsigned int iter                = 0;
    bool         continue_refinement = true;

    while (continue_refinement && (iter < max_iterations))
      {
        if (max_iterations != numbers::invalid_unsigned_int)
          iter++;
        continue_refinement = false;

        for (const auto &cell : tria.active_cell_iterators())
          for (const unsigned int j : cell->face_indices())
            if (cell->at_boundary(j) == false &&
                cell->neighbor(j)->has_children())
              {
                if (isotropic)
                  {
                    cell->set_refine_flag();
                    continue_refinement = true;
                  }
                else
                  continue_refinement |= cell->flag_for_face_refinement(j);
              }

        tria.execute_coarsening_and_refinement();
      }
  }

  template <int dim, int spacedim>
  void
  remove_anisotropy(Triangulation<dim, spacedim> &tria,
                    const double                  max_ratio,
                    const unsigned int            max_iterations)
  {
    unsigned int iter                = 0;
    bool         continue_refinement = true;

    while (continue_refinement && (iter < max_iterations))
      {
        iter++;
        continue_refinement = false;
        for (const auto &cell : tria.active_cell_iterators())
          {
            std::pair<unsigned int, double> info =
              GridTools::get_longest_direction<dim, spacedim>(cell);
            if (info.second > max_ratio)
              {
                cell->set_refine_flag(
                  RefinementCase<dim>::cut_axis(info.first));
                continue_refinement = true;
              }
          }
        tria.execute_coarsening_and_refinement();
      }
  }


  template <int dim, int spacedim>
  void
  regularize_corner_cells(Triangulation<dim, spacedim> &tria,
                          const double                  limit_angle_fraction)
  {
    if (dim == 1)
      return; // Nothing to do

    // Check that we don't have hanging nodes
    AssertThrow(!tria.has_hanging_nodes(),
                ExcMessage("The input Triangulation cannot "
                           "have hanging nodes."));


    bool has_cells_with_more_than_dim_faces_on_boundary = true;
    bool has_cells_with_dim_faces_on_boundary           = false;

    unsigned int refinement_cycles = 0;

    while (has_cells_with_more_than_dim_faces_on_boundary)
      {
        has_cells_with_more_than_dim_faces_on_boundary = false;

        for (const auto &cell : tria.active_cell_iterators())
          {
            unsigned int boundary_face_counter = 0;
            for (auto f : cell->face_indices())
              if (cell->face(f)->at_boundary())
                boundary_face_counter++;
            if (boundary_face_counter > dim)
              {
                has_cells_with_more_than_dim_faces_on_boundary = true;
                break;
              }
            else if (boundary_face_counter == dim)
              has_cells_with_dim_faces_on_boundary = true;
          }
        if (has_cells_with_more_than_dim_faces_on_boundary)
          {
            tria.refine_global(1);
            refinement_cycles++;
          }
      }

    if (has_cells_with_dim_faces_on_boundary)
      {
        tria.refine_global(1);
        refinement_cycles++;
      }
    else
      {
        while (refinement_cycles > 0)
          {
            for (const auto &cell : tria.active_cell_iterators())
              cell->set_coarsen_flag();
            tria.execute_coarsening_and_refinement();
            refinement_cycles--;
          }
        return;
      }

    std::vector<bool>            cells_to_remove(tria.n_active_cells(), false);
    std::vector<Point<spacedim>> vertices = tria.get_vertices();

    std::vector<bool> faces_to_remove(tria.n_raw_faces(), false);

    std::vector<CellData<dim>> cells_to_add;
    SubCellData                subcelldata_to_add;

    // Trick compiler for dimension independent things
    const unsigned int v0 = 0, v1 = 1, v2 = (dim > 1 ? 2 : 0),
                       v3 = (dim > 1 ? 3 : 0);

    for (const auto &cell : tria.active_cell_iterators())
      {
        double       angle_fraction   = 0;
        unsigned int vertex_at_corner = numbers::invalid_unsigned_int;

        if (dim == 2)
          {
            Tensor<1, spacedim> p0;
            p0[spacedim > 1 ? 1 : 0] = 1;
            Tensor<1, spacedim> p1;
            p1[0] = 1;

            if (cell->face(v0)->at_boundary() && cell->face(v3)->at_boundary())
              {
                p0               = cell->vertex(v0) - cell->vertex(v2);
                p1               = cell->vertex(v3) - cell->vertex(v2);
                vertex_at_corner = v2;
              }
            else if (cell->face(v3)->at_boundary() &&
                     cell->face(v1)->at_boundary())
              {
                p0               = cell->vertex(v2) - cell->vertex(v3);
                p1               = cell->vertex(v1) - cell->vertex(v3);
                vertex_at_corner = v3;
              }
            else if (cell->face(1)->at_boundary() &&
                     cell->face(2)->at_boundary())
              {
                p0               = cell->vertex(v0) - cell->vertex(v1);
                p1               = cell->vertex(v3) - cell->vertex(v1);
                vertex_at_corner = v1;
              }
            else if (cell->face(2)->at_boundary() &&
                     cell->face(0)->at_boundary())
              {
                p0               = cell->vertex(v2) - cell->vertex(v0);
                p1               = cell->vertex(v1) - cell->vertex(v0);
                vertex_at_corner = v0;
              }
            p0 /= p0.norm();
            p1 /= p1.norm();
            angle_fraction = std::acos(p0 * p1) / numbers::PI;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }

        if (angle_fraction > limit_angle_fraction)
          {
            auto flags_removal = [&](unsigned int f1,
                                     unsigned int f2,
                                     unsigned int n1,
                                     unsigned int n2) -> void {
              cells_to_remove[cell->active_cell_index()]               = true;
              cells_to_remove[cell->neighbor(n1)->active_cell_index()] = true;
              cells_to_remove[cell->neighbor(n2)->active_cell_index()] = true;

              faces_to_remove[cell->face(f1)->index()] = true;
              faces_to_remove[cell->face(f2)->index()] = true;

              faces_to_remove[cell->neighbor(n1)->face(f1)->index()] = true;
              faces_to_remove[cell->neighbor(n2)->face(f2)->index()] = true;
            };

            auto cell_creation = [&](const unsigned int vv0,
                                     const unsigned int vv1,
                                     const unsigned int f0,
                                     const unsigned int f1,

                                     const unsigned int n0,
                                     const unsigned int v0n0,
                                     const unsigned int v1n0,

                                     const unsigned int n1,
                                     const unsigned int v0n1,
                                     const unsigned int v1n1) {
              CellData<dim> c1, c2;
              CellData<1>   l1, l2;

              c1.vertices[v0] = cell->vertex_index(vv0);
              c1.vertices[v1] = cell->vertex_index(vv1);
              c1.vertices[v2] = cell->neighbor(n0)->vertex_index(v0n0);
              c1.vertices[v3] = cell->neighbor(n0)->vertex_index(v1n0);

              c1.manifold_id = cell->manifold_id();
              c1.material_id = cell->material_id();

              c2.vertices[v0] = cell->vertex_index(vv0);
              c2.vertices[v1] = cell->neighbor(n1)->vertex_index(v0n1);
              c2.vertices[v2] = cell->vertex_index(vv1);
              c2.vertices[v3] = cell->neighbor(n1)->vertex_index(v1n1);

              c2.manifold_id = cell->manifold_id();
              c2.material_id = cell->material_id();

              l1.vertices[0] = cell->vertex_index(vv0);
              l1.vertices[1] = cell->neighbor(n0)->vertex_index(v0n0);

              l1.boundary_id = cell->line(f0)->boundary_id();
              l1.manifold_id = cell->line(f0)->manifold_id();
              subcelldata_to_add.boundary_lines.push_back(l1);

              l2.vertices[0] = cell->vertex_index(vv0);
              l2.vertices[1] = cell->neighbor(n1)->vertex_index(v0n1);

              l2.boundary_id = cell->line(f1)->boundary_id();
              l2.manifold_id = cell->line(f1)->manifold_id();
              subcelldata_to_add.boundary_lines.push_back(l2);

              cells_to_add.push_back(c1);
              cells_to_add.push_back(c2);
            };

            if (dim == 2)
              {
                switch (vertex_at_corner)
                  {
                    case 0:
                      flags_removal(0, 2, 3, 1);
                      cell_creation(0, 3, 0, 2, 3, 2, 3, 1, 1, 3);
                      break;
                    case 1:
                      flags_removal(1, 2, 3, 0);
                      cell_creation(1, 2, 2, 1, 0, 0, 2, 3, 3, 2);
                      break;
                    case 2:
                      flags_removal(3, 0, 1, 2);
                      cell_creation(2, 1, 3, 0, 1, 3, 1, 2, 0, 1);
                      break;
                    case 3:
                      flags_removal(3, 1, 0, 2);
                      cell_creation(3, 0, 1, 3, 2, 1, 0, 0, 2, 0);
                      break;
                  }
              }
            else
              {
                Assert(false, ExcNotImplemented());
              }
          }
      }

    // if no cells need to be added, then no regularization is necessary.
    // Restore things as they were before this function was called.
    if (cells_to_add.size() == 0)
      {
        while (refinement_cycles > 0)
          {
            for (const auto &cell : tria.active_cell_iterators())
              cell->set_coarsen_flag();
            tria.execute_coarsening_and_refinement();
            refinement_cycles--;
          }
        return;
      }

    // add the cells that were not marked as skipped
    for (const auto &cell : tria.active_cell_iterators())
      {
        if (cells_to_remove[cell->active_cell_index()] == false)
          {
            CellData<dim> c;
            for (const unsigned int v : cell->vertex_indices())
              c.vertices[v] = cell->vertex_index(v);
            c.manifold_id = cell->manifold_id();
            c.material_id = cell->material_id();
            cells_to_add.push_back(c);
          }
      }

    // Face counter for both dim == 2 and dim == 3
    typename Triangulation<dim, spacedim>::active_face_iterator
      face = tria.begin_active_face(),
      endf = tria.end_face();
    for (; face != endf; ++face)
      if ((face->at_boundary() ||
           face->manifold_id() != numbers::flat_manifold_id) &&
          faces_to_remove[face->index()] == false)
        {
          for (unsigned int l = 0; l < face->n_lines(); ++l)
            {
              CellData<1> line;
              if (dim == 2)
                {
                  for (const unsigned int v : face->vertex_indices())
                    line.vertices[v] = face->vertex_index(v);
                  line.boundary_id = face->boundary_id();
                  line.manifold_id = face->manifold_id();
                }
              else
                {
                  for (const unsigned int v : face->line(l)->vertex_indices())
                    line.vertices[v] = face->line(l)->vertex_index(v);
                  line.boundary_id = face->line(l)->boundary_id();
                  line.manifold_id = face->line(l)->manifold_id();
                }
              subcelldata_to_add.boundary_lines.push_back(line);
            }
          if (dim == 3)
            {
              CellData<2> quad;
              for (const unsigned int v : face->vertex_indices())
                quad.vertices[v] = face->vertex_index(v);
              quad.boundary_id = face->boundary_id();
              quad.manifold_id = face->manifold_id();
              subcelldata_to_add.boundary_quads.push_back(quad);
            }
        }
    GridTools::delete_unused_vertices(vertices,
                                      cells_to_add,
                                      subcelldata_to_add);
    GridTools::consistently_order_cells(cells_to_add);

    // Save manifolds
    auto manifold_ids = tria.get_manifold_ids();
    std::map<types::manifold_id, std::unique_ptr<Manifold<dim, spacedim>>>
      manifolds;
    // Set manifolds in new Triangulation
    for (const auto manifold_id : manifold_ids)
      if (manifold_id != numbers::flat_manifold_id)
        manifolds[manifold_id] = tria.get_manifold(manifold_id).clone();

    tria.clear();

    tria.create_triangulation(vertices, cells_to_add, subcelldata_to_add);

    // Restore manifolds
    for (const auto manifold_id : manifold_ids)
      if (manifold_id != numbers::flat_manifold_id)
        tria.set_manifold(manifold_id, *manifolds[manifold_id]);
  }



  template <int dim, int spacedim>
#ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>>
#else
  return_type
#endif
  compute_point_locations(
    const Cache<dim, spacedim> &        cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint)
  {
    const auto cqmp = compute_point_locations_try_all(cache, points, cell_hint);
    // Splitting the tuple's components
    auto &cells          = std::get<0>(cqmp);
    auto &qpoints        = std::get<1>(cqmp);
    auto &maps           = std::get<2>(cqmp);
    auto &missing_points = std::get<3>(cqmp);
    // If a point was not found, throwing an error, as the old
    // implementation of compute_point_locations would have done
    AssertThrow(std::get<3>(cqmp).size() == 0,
                ExcPointNotFound<spacedim>(points[missing_points[0]]));

    (void)missing_points;

    return std::make_tuple(std::move(cells),
                           std::move(qpoints),
                           std::move(maps));
  }



  template <int dim, int spacedim>
#ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>,
    std::vector<unsigned int>>
#else
  return_type
#endif
  compute_point_locations_try_all(
    const Cache<dim, spacedim> &        cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint)
  {
    // How many points are here?
    const unsigned int np = points.size();

    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>
                                           cells_out;
    std::vector<std::vector<Point<dim>>>   qpoints_out;
    std::vector<std::vector<unsigned int>> maps_out;
    std::vector<unsigned int>              missing_points_out;

    // Now the easy case.
    if (np == 0)
      return std::make_tuple(std::move(cells_out),
                             std::move(qpoints_out),
                             std::move(maps_out),
                             std::move(missing_points_out));

    // For the search we shall use the following tree
    const auto &b_tree = cache.get_cell_bounding_boxes_rtree();

    // We begin by finding the cell/transform of the first point
    std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
              Point<dim>>
      my_pair;

    bool         found          = false;
    unsigned int points_checked = 0;

    // If a hint cell was given, use it
    if (cell_hint.state() == IteratorState::valid)
      {
        try
          {
            my_pair = GridTools::find_active_cell_around_point(cache,
                                                               points[0],
                                                               cell_hint);
            found   = true;
          }
        catch (const GridTools::ExcPointNotFound<dim> &)
          {
            missing_points_out.emplace_back(0);
          }
        ++points_checked;
      }

    // The tree search returns
    // - a bounding box covering the cell
    // - the active cell iterator
    std::vector<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>>
      box_cell;

    // This is used as an index for box_cell
    int cell_candidate_idx = -1;
    // If any of the cells in box_cell is a ghost cell,
    // an artificial cell or at the boundary,
    // we want to use try/catch
    bool use_try = false;

    while (!found && points_checked < np)
      {
        box_cell.clear();
        b_tree.query(boost::geometry::index::intersects(points[points_checked]),
                     std::back_inserter(box_cell));

        // Checking box_cell result for a suitable candidate
        cell_candidate_idx = -1;
        for (unsigned int i = 0; i < box_cell.size(); ++i)
          {
            // As a candidate we don't want artificial cells
            if (!box_cell[i].second->is_artificial())
              cell_candidate_idx = i;

            // If the cell is not locally owned or at boundary
            // we check for exceptions
            if (cell_candidate_idx != -1 &&
                (!box_cell[i].second->is_locally_owned() ||
                 box_cell[i].second->at_boundary()))
              use_try = true;


            if (cell_candidate_idx != -1)
              break;
          }

        // If a suitable cell was found, use it as hint
        if (cell_candidate_idx != -1)
          {
            if (use_try)
              {
                try
                  {
                    my_pair = GridTools::find_active_cell_around_point(
                      cache,
                      points[points_checked],
                      box_cell[cell_candidate_idx].second);
                    found = true;
                  }
                catch (const GridTools::ExcPointNotFound<dim> &)
                  {
                    missing_points_out.emplace_back(points_checked);
                  }
              }
            else
              {
                my_pair = GridTools::find_active_cell_around_point(
                  cache,
                  points[points_checked],
                  box_cell[cell_candidate_idx].second);
                found = true;
              }
          }
        else
          {
            try
              {
                my_pair = GridTools::find_active_cell_around_point(
                  cache, points[points_checked]);
                // If we arrive here the cell was not among
                // the candidates returned by the tree, so we're adding it
                // by hand
                found              = true;
                cell_candidate_idx = box_cell.size();
                box_cell.push_back(
                  std::make_pair(my_pair.first->bounding_box(), my_pair.first));
              }
            catch (const GridTools::ExcPointNotFound<dim> &)
              {
                missing_points_out.emplace_back(points_checked);
              }
          }

        // Updating the position of the analyzed points
        ++points_checked;
      }

    // If the point has been found in a cell, adding it
    if (found)
      {
        cells_out.emplace_back(my_pair.first);
        qpoints_out.emplace_back(1, my_pair.second);
        maps_out.emplace_back(1, points_checked - 1);
      }

    // Now the second easy case.
    if (np == qpoints_out.size())
      return std::make_tuple(std::move(cells_out),
                             std::move(qpoints_out),
                             std::move(maps_out),
                             std::move(missing_points_out));

    // Cycle over all points left
    for (unsigned int p = points_checked; p < np; ++p)
      {
        // We assume the last used cell contains the point: checking it
        if (cell_candidate_idx != -1)
          if (!box_cell[cell_candidate_idx].first.point_inside(points[p]))
            // Point outside candidate cell: we have no candidate
            cell_candidate_idx = -1;

        // If there's no candidate, run a tree search
        if (cell_candidate_idx == -1)
          {
            // Using the b_tree to find new candidates
            box_cell.clear();
            b_tree.query(boost::geometry::index::intersects(points[p]),
                         std::back_inserter(box_cell));
            // Checking the returned bounding boxes/cells
            use_try            = false;
            cell_candidate_idx = -1;
            for (unsigned int i = 0; i < box_cell.size(); ++i)
              {
                // As a candidate we don't want artificial cells
                if (!box_cell[i].second->is_artificial())
                  cell_candidate_idx = i;

                // If the cell is not locally owned or at boundary
                // we check for exceptions
                if (cell_candidate_idx != -1 &&
                    (!box_cell[i].second->is_locally_owned() ||
                     box_cell[i].second->at_boundary()))
                  use_try = true;

                // If a cell candidate was found we can stop
                if (cell_candidate_idx != -1)
                  break;
              }
          }

        if (cell_candidate_idx == -1)
          {
            // No candidate cell, but the cell might
            // still be inside the mesh, this is our final check:
            try
              {
                my_pair =
                  GridTools::find_active_cell_around_point(cache, points[p]);
                // If we arrive here the cell was not among
                // the candidates returned by the tree, so we're adding it
                // by hand
                cell_candidate_idx = box_cell.size();
                box_cell.push_back(
                  std::make_pair(my_pair.first->bounding_box(), my_pair.first));
              }
            catch (const GridTools::ExcPointNotFound<dim> &)
              {
                missing_points_out.emplace_back(p);
                continue;
              }
          }
        else
          {
            // We have a candidate cell
            if (use_try)
              {
                try
                  {
                    my_pair = GridTools::find_active_cell_around_point(
                      cache, points[p], box_cell[cell_candidate_idx].second);
                  }
                catch (const GridTools::ExcPointNotFound<dim> &)
                  {
                    missing_points_out.push_back(p);
                    continue;
                  }
              }
            else
              {
                my_pair = GridTools::find_active_cell_around_point(
                  cache, points[p], box_cell[cell_candidate_idx].second);
              }

            // If the point was found in another cell,
            // updating cell_candidate_idx
            if (my_pair.first != box_cell[cell_candidate_idx].second)
              {
                for (unsigned int i = 0; i < box_cell.size(); ++i)
                  {
                    if (my_pair.first == box_cell[i].second)
                      {
                        cell_candidate_idx = i;
                        break;
                      }
                  }

                if (my_pair.first != box_cell[cell_candidate_idx].second)
                  {
                    // The cell was not among the candidates returned by the
                    // tree
                    cell_candidate_idx = box_cell.size();
                    box_cell.push_back(
                      std::make_pair(my_pair.first->bounding_box(),
                                     my_pair.first));
                  }
              }
          }


        // Assuming the point is more likely to be in the last
        // used cell
        if (my_pair.first == cells_out.back())
          {
            // Found in the last cell: adding the data
            qpoints_out.back().emplace_back(my_pair.second);
            maps_out.back().emplace_back(p);
          }
        else
          {
            // Check if it is in another cell already found
            typename std::vector<typename Triangulation<dim, spacedim>::
                                   active_cell_iterator>::iterator cells_it =
              std::find(cells_out.begin(), cells_out.end() - 1, my_pair.first);

            if (cells_it == cells_out.end() - 1)
              {
                // Cell not found: adding a new cell
                cells_out.emplace_back(my_pair.first);
                qpoints_out.emplace_back(1, my_pair.second);
                maps_out.emplace_back(1, p);
              }
            else
              {
                // Cell found: just adding the point index and qpoint to the
                // list
                unsigned int current_cell = cells_it - cells_out.begin();
                qpoints_out[current_cell].emplace_back(my_pair.second);
                maps_out[current_cell].emplace_back(p);
              }
          }
      }

    // Debug Checking
    Assert(cells_out.size() == maps_out.size(),
           ExcDimensionMismatch(cells_out.size(), maps_out.size()));

    Assert(cells_out.size() == qpoints_out.size(),
           ExcDimensionMismatch(cells_out.size(), qpoints_out.size()));

#ifdef DEBUG
    unsigned int c   = cells_out.size();
    unsigned int qps = 0;
    // The number of points in all
    // the cells must be the same as
    // the number of points we
    // started off from,
    // plus the points which were ignored
    for (unsigned int n = 0; n < c; ++n)
      {
        Assert(qpoints_out[n].size() == maps_out[n].size(),
               ExcDimensionMismatch(qpoints_out[n].size(), maps_out[n].size()));
        qps += qpoints_out[n].size();
      }

    Assert(qps + missing_points_out.size() == np,
           ExcDimensionMismatch(qps + missing_points_out.size(), np));
#endif

    return std::make_tuple(std::move(cells_out),
                           std::move(qpoints_out),
                           std::move(maps_out),
                           std::move(missing_points_out));
  }



  template <int dim, int spacedim>
#ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>,
    std::vector<std::vector<Point<spacedim>>>,
    std::vector<std::vector<unsigned int>>>
#else
  return_type
#endif
  distributed_compute_point_locations(
    const GridTools::Cache<dim, spacedim> &                cache,
    const std::vector<Point<spacedim>> &                   points,
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const double                                           tolerance)
  {
    // run internal function ...
    const auto all = internal::distributed_compute_point_locations(
                       cache, points, global_bboxes, tolerance, false)
                       .send_components;

    // ... and reshuffle the data
    std::tuple<
      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
      std::vector<std::vector<Point<dim>>>,
      std::vector<std::vector<unsigned int>>,
      std::vector<std::vector<Point<spacedim>>>,
      std::vector<std::vector<unsigned int>>>
      result;

    std::pair<int, int> dummy{-1, -1};

    for (unsigned int i = 0; i < all.size(); ++i)
      {
        if (dummy != std::get<0>(all[i]))
          {
            std::get<0>(result).push_back(
              typename Triangulation<dim, spacedim>::active_cell_iterator{
                &cache.get_triangulation(),
                std::get<0>(all[i]).first,
                std::get<0>(all[i]).second});

            const unsigned int new_size = std::get<0>(result).size();

            std::get<1>(result).resize(new_size);
            std::get<2>(result).resize(new_size);
            std::get<3>(result).resize(new_size);
            std::get<4>(result).resize(new_size);

            dummy = std::get<0>(all[i]);
          }

        std::get<1>(result).back().push_back(
          std::get<3>(all[i])); // reference point
        std::get<2>(result).back().push_back(std::get<2>(all[i])); // index
        std::get<3>(result).back().push_back(std::get<4>(all[i])); // real point
        std::get<4>(result).back().push_back(std::get<1>(all[i])); // rank
      }

    return result;
  }



  namespace internal
  {
    template <int spacedim>
    std::tuple<std::vector<unsigned int>,
               std::vector<unsigned int>,
               std::vector<unsigned int>>
    guess_point_owner(
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const std::vector<Point<spacedim>> &                   points)
    {
      std::vector<std::pair<unsigned int, unsigned int>> ranks_and_indices;
      ranks_and_indices.reserve(points.size());

      for (unsigned int i = 0; i < points.size(); ++i)
        {
          const auto &point = points[i];
          for (unsigned rank = 0; rank < global_bboxes.size(); ++rank)
            for (const auto &box : global_bboxes[rank])
              if (box.point_inside(point))
                {
                  ranks_and_indices.emplace_back(rank, i);
                  break;
                }
        }

      // convert to CRS
      std::sort(ranks_and_indices.begin(), ranks_and_indices.end());

      std::vector<unsigned int> ranks;
      std::vector<unsigned int> ptr;
      std::vector<unsigned int> indices;

      unsigned int dummy_rank = numbers::invalid_unsigned_int;

      for (const auto &i : ranks_and_indices)
        {
          if (dummy_rank != i.first)
            {
              dummy_rank = i.first;
              ranks.push_back(dummy_rank);
              ptr.push_back(indices.size());
            }

          indices.push_back(i.second);
        }
      ptr.push_back(indices.size());

      return std::make_tuple(std::move(ranks),
                             std::move(ptr),
                             std::move(indices));
    }



    template <int dim, int spacedim>
    std::vector<
      std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                Point<dim>>>
    find_all_locally_owned_active_cells_around_point(
      const Cache<dim, spacedim> &                                 cache,
      const Point<spacedim> &                                      point,
      typename Triangulation<dim, spacedim>::active_cell_iterator &cell_hint,
      const std::vector<bool> &marked_vertices,
      const double             tolerance)
    {
      std::vector<
        std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                  Point<dim>>>
        locally_owned_active_cells_around_point;

      try
        {
          const auto first_cell = GridTools::find_active_cell_around_point(
            cache, point, cell_hint, marked_vertices, tolerance);

          cell_hint = first_cell.first;

          const auto active_cells_around_point =
            GridTools::find_all_active_cells_around_point(
              cache.get_mapping(),
              cache.get_triangulation(),
              point,
              tolerance,
              first_cell);

          locally_owned_active_cells_around_point.reserve(
            active_cells_around_point.size());

          for (const auto &cell : active_cells_around_point)
            if (cell.first->is_locally_owned())
              locally_owned_active_cells_around_point.push_back(cell);
        }
      catch (...)
        {}

      return locally_owned_active_cells_around_point;
    }



    template <int dim, int spacedim>
    DistributedComputePointLocationsInternal<dim, spacedim>
    distributed_compute_point_locations(
      const GridTools::Cache<dim, spacedim> &                cache,
      const std::vector<Point<spacedim>> &                   points,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const double                                           tolerance,
      const bool                                             perform_handshake)
    {
      DistributedComputePointLocationsInternal<dim, spacedim> result;

      auto &send_components = result.send_components;
      auto &send_ranks      = result.send_ranks;
      auto &send_ptrs       = result.send_ptrs;
      auto &recv_components = result.recv_components;
      auto &recv_ranks      = result.recv_ranks;
      auto &recv_ptrs       = result.recv_ptrs;

      const auto potential_owners =
        internal::guess_point_owner(global_bboxes, points);

      const auto &potential_owners_ranks   = std::get<0>(potential_owners);
      const auto &potential_owners_ptrs    = std::get<1>(potential_owners);
      const auto &potential_owners_indices = std::get<2>(potential_owners);

      const std::vector<bool> marked_vertices;
      auto cell_hint = cache.get_triangulation().begin_active();

      const auto translate = [&](const unsigned int other_rank) {
        const auto ptr = std::find(potential_owners_ranks.begin(),
                                   potential_owners_ranks.end(),
                                   other_rank);

        Assert(ptr != potential_owners_ranks.end(), ExcInternalError());

        const auto other_rank_index =
          std::distance(potential_owners_ranks.begin(), ptr);

        return other_rank_index;
      };

      Utilities::MPI::ConsensusAlgorithms::AnonymousProcess<char, char> process(
        [&]() { return potential_owners_ranks; },
        [&](const unsigned int other_rank, std::vector<char> &send_buffer) {
          const auto other_rank_index = translate(other_rank);

          std::vector<std::pair<unsigned int, Point<spacedim>>> temp;
          temp.reserve(potential_owners_ptrs[other_rank_index + 1] -
                       potential_owners_ptrs[other_rank_index]);

          for (unsigned int i = potential_owners_ptrs[other_rank_index];
               i < potential_owners_ptrs[other_rank_index + 1];
               ++i)
            temp.emplace_back(potential_owners_indices[i],
                              points[potential_owners_indices[i]]);

          send_buffer = Utilities::pack(temp, false);
        },
        [&](const unsigned int &     other_rank,
            const std::vector<char> &recv_buffer,
            std::vector<char> &      request_buffer) {
          const auto recv_buffer_unpacked = Utilities::unpack<
            std::vector<std::pair<unsigned int, Point<spacedim>>>>(recv_buffer,
                                                                   false);

          std::vector<unsigned int> request_buffer_temp(
            recv_buffer_unpacked.size(), 0);

          cell_hint = cache.get_triangulation().begin_active();

          for (unsigned int i = 0; i < recv_buffer_unpacked.size(); ++i)
            {
              const auto &index_and_point = recv_buffer_unpacked[i];

              const auto cells_and_reference_positions =
                find_all_locally_owned_active_cells_around_point(
                  cache,
                  index_and_point.second,
                  cell_hint,
                  marked_vertices,
                  tolerance);

              for (const auto &cell_and_reference_position :
                   cells_and_reference_positions)
                {
                  send_components.emplace_back(
                    std::pair<int, int>(
                      cell_and_reference_position.first->level(),
                      cell_and_reference_position.first->index()),
                    other_rank,
                    index_and_point.first,
                    cell_and_reference_position.second,
                    index_and_point.second,
                    numbers::invalid_unsigned_int);
                }

              if (perform_handshake)
                request_buffer_temp[i] = cells_and_reference_positions.size();
            }

          if (perform_handshake)
            request_buffer = Utilities::pack(request_buffer_temp, false);
        },
        [&](const unsigned int other_rank, std::vector<char> &recv_buffer) {
          if (perform_handshake)
            {
              const auto other_rank_index = translate(other_rank);

              recv_buffer =
                Utilities::pack(std::vector<unsigned int>(
                                  potential_owners_ptrs[other_rank_index + 1] -
                                  potential_owners_ptrs[other_rank_index]),
                                false);
            }
        },
        [&](const unsigned int       other_rank,
            const std::vector<char> &recv_buffer) {
          if (perform_handshake)
            {
              const auto recv_buffer_unpacked =
                Utilities::unpack<std::vector<unsigned int>>(recv_buffer,
                                                             false);

              const auto other_rank_index = translate(other_rank);

              for (unsigned int i = 0; i < recv_buffer_unpacked.size(); ++i)
                for (unsigned int j = 0; j < recv_buffer_unpacked[i]; ++j)
                  recv_components.emplace_back(
                    other_rank,
                    potential_owners_indices
                      [i + potential_owners_ptrs[other_rank_index]],
                    numbers::invalid_unsigned_int);
            }
        });

      Utilities::MPI::ConsensusAlgorithms::Selector<char, char>(
        process, cache.get_triangulation().get_communicator())
        .run();

      if (true)
        {
          // sort according to rank (and point index and cell) -> make
          // deterministic
          std::sort(send_components.begin(),
                    send_components.end(),
                    [&](const auto &a, const auto &b) {
                      if (std::get<1>(a) != std::get<1>(b)) // rank
                        return std::get<1>(a) < std::get<1>(b);

                      if (std::get<2>(a) != std::get<2>(b)) // point index
                        return std::get<2>(a) < std::get<2>(b);

                      return std::get<0>(a) < std::get<0>(b); // cell
                    });

          // perform enumeration and extract rank information
          for (unsigned int i = 0, dummy = numbers::invalid_unsigned_int;
               i < send_components.size();
               ++i)
            {
              std::get<5>(send_components[i]) = i;

              if (dummy != std::get<1>(send_components[i]))
                {
                  dummy = std::get<1>(send_components[i]);
                  send_ranks.push_back(dummy);
                  send_ptrs.push_back(i);
                }
            }
          send_ptrs.push_back(send_components.size());

          // sort according to cell, rank, point index (while keeping
          // partial ordering)
          std::sort(send_components.begin(),
                    send_components.end(),
                    [&](const auto &a, const auto &b) {
                      if (std::get<0>(a) != std::get<0>(b))
                        return std::get<0>(a) < std::get<0>(b); // cell

                      if (std::get<1>(a) != std::get<1>(b))
                        return std::get<1>(a) < std::get<1>(b); // rank

                      if (std::get<2>(a) != std::get<2>(b))
                        return std::get<2>(a) < std::get<2>(b); // point index

                      return std::get<5>(a) < std::get<5>(b); // enumeration
                    });
        }

      if (perform_handshake)
        {
          // sort according to rank (and point index) -> make deterministic
          std::sort(recv_components.begin(),
                    recv_components.end(),
                    [&](const auto &a, const auto &b) {
                      if (std::get<0>(a) != std::get<0>(b))
                        return std::get<0>(a) < std::get<0>(b); // rank

                      return std::get<1>(a) < std::get<1>(b); // point index
                    });

          // perform enumeration and extract rank information
          for (unsigned int i = 0, dummy = numbers::invalid_unsigned_int;
               i < recv_components.size();
               ++i)
            {
              std::get<2>(recv_components[i]) = i;

              if (dummy != std::get<0>(recv_components[i]))
                {
                  dummy = std::get<0>(recv_components[i]);
                  recv_ranks.push_back(dummy);
                  recv_ptrs.push_back(i);
                }
            }
          recv_ptrs.push_back(recv_components.size());

          // sort according to point index and rank (while keeping partial
          // ordering)
          std::sort(recv_components.begin(),
                    recv_components.end(),
                    [&](const auto &a, const auto &b) {
                      if (std::get<1>(a) != std::get<1>(b))
                        return std::get<1>(a) < std::get<1>(b); // point index

                      if (std::get<0>(a) != std::get<0>(b))
                        return std::get<0>(a) < std::get<0>(b); // rank

                      return std::get<2>(a) < std::get<2>(b); // enumeration
                    });
        }

      return result;
    }
  } // namespace internal



  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  extract_used_vertices(const Triangulation<dim, spacedim> &container,
                        const Mapping<dim, spacedim> &      mapping)
  {
    std::map<unsigned int, Point<spacedim>> result;
    for (const auto &cell : container.active_cell_iterators())
      {
        if (!cell->is_artificial())
          {
            const auto vs = mapping.get_vertices(cell);
            for (unsigned int i = 0; i < vs.size(); ++i)
              result[cell->vertex_index(i)] = vs[i];
          }
      }
    return result;
  }


  template <int spacedim>
  unsigned int
  find_closest_vertex(const std::map<unsigned int, Point<spacedim>> &vertices,
                      const Point<spacedim> &                        p)
  {
    auto id_and_v = std::min_element(
      vertices.begin(),
      vertices.end(),
      [&](const std::pair<const unsigned int, Point<spacedim>> &p1,
          const std::pair<const unsigned int, Point<spacedim>> &p2) -> bool {
        return p1.second.distance(p) < p2.second.distance(p);
      });
    return id_and_v->first;
  }


  template <int dim, int spacedim>
  std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
            Point<dim>>
  find_active_cell_around_point(
    const Cache<dim, spacedim> &cache,
    const Point<spacedim> &     p,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &                      cell_hint,
    const std::vector<bool> &marked_vertices,
    const double             tolerance)
  {
    const auto &mesh            = cache.get_triangulation();
    const auto &mapping         = cache.get_mapping();
    const auto &vertex_to_cells = cache.get_vertex_to_cell_map();
    const auto &vertex_to_cell_centers =
      cache.get_vertex_to_cell_centers_directions();
    const auto &used_vertices_rtree = cache.get_used_vertices_rtree();

    return find_active_cell_around_point(mapping,
                                         mesh,
                                         p,
                                         vertex_to_cells,
                                         vertex_to_cell_centers,
                                         cell_hint,
                                         marked_vertices,
                                         used_vertices_rtree,
                                         tolerance);
  }

  template <int spacedim>
  std::vector<std::vector<BoundingBox<spacedim>>>
  exchange_local_bounding_boxes(
    const std::vector<BoundingBox<spacedim>> &local_bboxes,
    const MPI_Comm &                          mpi_communicator)
  {
#ifndef DEAL_II_WITH_MPI
    (void)local_bboxes;
    (void)mpi_communicator;
    Assert(false,
           ExcMessage(
             "GridTools::exchange_local_bounding_boxes() requires MPI."));
    return {};
#else
    // Step 1: preparing data to be sent
    unsigned int n_bboxes = local_bboxes.size();
    // Dimension of the array to be exchanged (number of double)
    int n_local_data = 2 * spacedim * n_bboxes;
    // data array stores each entry of each point describing the bounding boxes
    std::vector<double> loc_data_array(n_local_data);
    for (unsigned int i = 0; i < n_bboxes; ++i)
      for (unsigned int d = 0; d < spacedim; ++d)
        {
          // Extracting the coordinates of each boundary point
          loc_data_array[2 * i * spacedim + d] =
            local_bboxes[i].get_boundary_points().first[d];
          loc_data_array[2 * i * spacedim + spacedim + d] =
            local_bboxes[i].get_boundary_points().second[d];
        }

    // Step 2: exchanging the size of local data
    unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);

    // Vector to store the size of loc_data_array for every process
    std::vector<int> size_all_data(n_procs);

    // Exchanging the number of bboxes
    int ierr = MPI_Allgather(&n_local_data,
                             1,
                             MPI_INT,
                             size_all_data.data(),
                             1,
                             MPI_INT,
                             mpi_communicator);
    AssertThrowMPI(ierr);

    // Now computing the the displacement, relative to recvbuf,
    // at which to store the incoming data
    std::vector<int> rdispls(n_procs);
    rdispls[0] = 0;
    for (unsigned int i = 1; i < n_procs; ++i)
      rdispls[i] = rdispls[i - 1] + size_all_data[i - 1];

    // Step 3: exchange the data and bounding boxes:
    // Allocating a vector to contain all the received data
    std::vector<double> data_array(rdispls.back() + size_all_data.back());

    ierr = MPI_Allgatherv(loc_data_array.data(),
                          n_local_data,
                          MPI_DOUBLE,
                          data_array.data(),
                          size_all_data.data(),
                          rdispls.data(),
                          MPI_DOUBLE,
                          mpi_communicator);
    AssertThrowMPI(ierr);

    // Step 4: create the array of bboxes for output
    std::vector<std::vector<BoundingBox<spacedim>>> global_bboxes(n_procs);
    unsigned int                                    begin_idx = 0;
    for (unsigned int i = 0; i < n_procs; ++i)
      {
        // Number of local bounding boxes
        unsigned int n_bbox_i = size_all_data[i] / (spacedim * 2);
        global_bboxes[i].resize(n_bbox_i);
        for (unsigned int bbox = 0; bbox < n_bbox_i; ++bbox)
          {
            Point<spacedim> p1, p2; // boundary points for bbox
            for (unsigned int d = 0; d < spacedim; ++d)
              {
                p1[d] = data_array[begin_idx + 2 * bbox * spacedim + d];
                p2[d] =
                  data_array[begin_idx + 2 * bbox * spacedim + spacedim + d];
              }
            BoundingBox<spacedim> loc_bbox(std::make_pair(p1, p2));
            global_bboxes[i][bbox] = loc_bbox;
          }
        // Shifting the first index to the start of the next vector
        begin_idx += size_all_data[i];
      }
    return global_bboxes;
#endif // DEAL_II_WITH_MPI
  }



  template <int spacedim>
  RTree<std::pair<BoundingBox<spacedim>, unsigned int>>
  build_global_description_tree(
    const std::vector<BoundingBox<spacedim>> &local_description,
    const MPI_Comm &                          mpi_communicator)
  {
#ifndef DEAL_II_WITH_MPI
    (void)mpi_communicator;
    // Building a tree with the only boxes available without MPI
    std::vector<std::pair<BoundingBox<spacedim>, unsigned int>> boxes_index(
      local_description.size());
    // Adding to each box the rank of the process owning it
    for (unsigned int i = 0; i < local_description.size(); ++i)
      boxes_index[i] = std::make_pair(local_description[i], 0u);
    return pack_rtree(boxes_index);
#else
    // Exchanging local bounding boxes
    const std::vector<std::vector<BoundingBox<spacedim>>> global_bboxes =
      Utilities::MPI::all_gather(mpi_communicator, local_description);

    // Preparing to flatten the vector
    const unsigned int n_procs =
      Utilities::MPI::n_mpi_processes(mpi_communicator);
    // The i'th element of the following vector contains the index of the first
    // local bounding box from the process of rank i
    std::vector<unsigned int> bboxes_position(n_procs);

    unsigned int tot_bboxes = 0;
    for (const auto &process_bboxes : global_bboxes)
      tot_bboxes += process_bboxes.size();

    // Now flattening the vector
    std::vector<std::pair<BoundingBox<spacedim>, unsigned int>>
      flat_global_bboxes;
    flat_global_bboxes.reserve(tot_bboxes);
    unsigned int process_index = 0;
    for (const auto &process_bboxes : global_bboxes)
      {
        // Initialize a vector containing bounding boxes and rank of a process
        std::vector<std::pair<BoundingBox<spacedim>, unsigned int>>
          boxes_and_indices(process_bboxes.size());

        // Adding to each box the rank of the process owning it
        for (unsigned int i = 0; i < process_bboxes.size(); ++i)
          boxes_and_indices[i] =
            std::make_pair(process_bboxes[i], process_index);

        flat_global_bboxes.insert(flat_global_bboxes.end(),
                                  boxes_and_indices.begin(),
                                  boxes_and_indices.end());

        ++process_index;
      }

    // Build a tree out of the bounding boxes.  We avoid using the
    // insert method so that boost uses the packing algorithm
    return RTree<std::pair<BoundingBox<spacedim>, unsigned int>>(
      flat_global_bboxes.begin(), flat_global_bboxes.end());
#endif // DEAL_II_WITH_MPI
  }



  template <int dim, int spacedim>
  void
  collect_coinciding_vertices(
    const Triangulation<dim, spacedim> &               tria,
    std::map<unsigned int, std::vector<unsigned int>> &coinciding_vertex_groups,
    std::map<unsigned int, unsigned int> &vertex_to_coinciding_vertex_group)
  {
    // 1) determine for each vertex a vertex it concides with and
    //    put it into a map
    {
      static const int lookup_table_2d[2][2] =
        //           flip:
        {
          {0, 1}, // false
          {1, 0}  // true
        };

      static const int lookup_table_3d[2][2][2][4] =
        //                   orientation flip  rotation
        {{{
            {0, 2, 1, 3}, // false       false false
            {2, 3, 0, 1}  // false       false true
          },
          {
            {3, 1, 2, 0}, // false       true  false
            {1, 0, 3, 2}  // false       true  true
          }},
         {{
            {0, 1, 2, 3}, // true        false false
            {1, 3, 0, 2}  // true        false true
          },
          {
            {3, 2, 1, 0}, // true        true  false
            {2, 0, 3, 1}  // true        true  true
          }}};

      // loop over all periodic face pairs
      for (const auto &pair : tria.get_periodic_face_map())
        {
          if (pair.first.first->level() != pair.second.first.first->level())
            continue;

          const auto face_a = pair.first.first->face(pair.first.second);
          const auto face_b =
            pair.second.first.first->face(pair.second.first.second);
          const auto mask = pair.second.second;

          AssertDimension(face_a->n_vertices(), face_b->n_vertices());

          // loop over all vertices on face
          for (unsigned int i = 0; i < face_a->n_vertices(); ++i)
            {
              const bool face_orientation = mask[0];
              const bool face_flip        = mask[1];
              const bool face_rotation    = mask[2];

              // find the right local vertex index for the second face
              unsigned int j = 0;
              switch (dim)
                {
                  case 1:
                    j = i;
                    break;
                  case 2:
                    j = lookup_table_2d[face_flip][i];
                    break;
                  case 3:
                    j = lookup_table_3d[face_orientation][face_flip]
                                       [face_rotation][i];
                    break;
                  default:
                    AssertThrow(false, ExcNotImplemented());
                }

              // get vertex indices and store in map
              const auto   vertex_a = face_a->vertex_index(i);
              const auto   vertex_b = face_b->vertex_index(j);
              unsigned int temp     = std::min(vertex_a, vertex_b);

              auto it_a = vertex_to_coinciding_vertex_group.find(vertex_a);
              if (it_a != vertex_to_coinciding_vertex_group.end())
                temp = std::min(temp, it_a->second);

              auto it_b = vertex_to_coinciding_vertex_group.find(vertex_b);
              if (it_b != vertex_to_coinciding_vertex_group.end())
                temp = std::min(temp, it_b->second);

              if (it_a != vertex_to_coinciding_vertex_group.end())
                it_a->second = temp;
              else
                vertex_to_coinciding_vertex_group[vertex_a] = temp;

              if (it_b != vertex_to_coinciding_vertex_group.end())
                it_b->second = temp;
              else
                vertex_to_coinciding_vertex_group[vertex_b] = temp;
            }
        }

      // 2) compress map: let vertices point to the coinciding vertex with
      //    the smallest index
      for (auto &p : vertex_to_coinciding_vertex_group)
        {
          if (p.first == p.second)
            continue;
          unsigned int temp = p.second;
          while (temp != vertex_to_coinciding_vertex_group[temp])
            temp = vertex_to_coinciding_vertex_group[temp];
          p.second = temp;
        }

      // 3) create a map: smallest index of coinciding index -> all
      //    coinciding indices
      for (auto p : vertex_to_coinciding_vertex_group)
        coinciding_vertex_groups[p.second] = {};

      for (auto p : vertex_to_coinciding_vertex_group)
        coinciding_vertex_groups[p.second].push_back(p.first);
    }
  }



  template <int dim, int spacedim>
  std::map<unsigned int, std::set<dealii::types::subdomain_id>>
  compute_vertices_with_ghost_neighbors(
    const Triangulation<dim, spacedim> &tria)
  {
    if (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &tria) == nullptr) // nothing to do for a serial triangulation
      return {};

    // 1) collect for each vertex on periodic faces all vertices it coincides
    //    with
    std::map<unsigned int, std::vector<unsigned int>> coinciding_vertex_groups;
    std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;

    GridTools::collect_coinciding_vertices(tria,
                                           coinciding_vertex_groups,
                                           vertex_to_coinciding_vertex_group);

    // 2) collect vertices belonging to local cells
    std::vector<bool> vertex_of_own_cell(tria.n_vertices(), false);
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const unsigned int v : cell->vertex_indices())
          vertex_of_own_cell[cell->vertex_index(v)] = true;

    // 3) for each vertex belonging to a locally owned cell all ghost
    //    neighbors (including the periodic own)
    std::map<unsigned int, std::set<types::subdomain_id>> result;

    // loop over all active ghost cells
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->is_ghost())
        {
          const types::subdomain_id owner = cell->subdomain_id();

          // loop over all its vertices
          for (const unsigned int v : cell->vertex_indices())
            {
              // set owner if vertex belongs to a local cell
              if (vertex_of_own_cell[cell->vertex_index(v)])
                result[cell->vertex_index(v)].insert(owner);

              // mark also nodes coinciding due to periodicity
              auto coinciding_vertex_group =
                vertex_to_coinciding_vertex_group.find(cell->vertex_index(v));
              if (coinciding_vertex_group !=
                  vertex_to_coinciding_vertex_group.end())
                for (auto coinciding_vertex :
                     coinciding_vertex_groups[coinciding_vertex_group->second])
                  if (vertex_of_own_cell[coinciding_vertex])
                    result[coinciding_vertex].insert(owner);
            }
        }

    return result;
  }

} /* namespace GridTools */


// explicit instantiations
#include "grid_tools.inst"

DEAL_II_NAMESPACE_CLOSE

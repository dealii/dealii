// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#ifdef DEAL_II_WITH_ARBORX
#  include <deal.II/arborx/access_traits.h>
#  include <deal.II/arborx/distributed_tree.h>
#endif

#ifdef DEAL_II_WITH_CGAL
#  include <deal.II/cgal/intersections.h>
#  include <deal.II/cgal/utilities.h>
#endif

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
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/constrained_linear_operator.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>

#include <deal.II/physics/transformations.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_map>

DEAL_II_NAMESPACE_OPEN

#ifndef DEAL_II_WITH_ARBORX

// If we configured without ArborX, we still need to have a couple of
// dummy types that we can reference in code below. They do not
// actually do anything useful.
template <int dim, typename Number>
class BoundingBox;

namespace ArborXWrappers
{
  class DistributedTree
  {
  public:
    template <int dim, typename Number>
    DistributedTree(const MPI_Comm &,
                    const std::vector<BoundingBox<dim, Number>> &);

    template <typename QueryType>
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>>
    query(const QueryType &queries);
  };

  class BoundingBoxIntersectPredicate
  {};
} // namespace ArborXWrappers
#endif


namespace GridTools
{
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


    // Transformation to rotate around one of the cartesian z-axis in 2d.
    class Rotate2d
    {
    public:
      explicit Rotate2d(const double angle)
        : rotation_matrix(
            Physics::Transformations::Rotations::rotation_matrix_2d(angle))
      {}
      Point<2>
      operator()(const Point<2> &p) const
      {
        return static_cast<Point<2>>(rotation_matrix * p);
      }

    private:
      const Tensor<2, 2, double> rotation_matrix;
    };


    // Transformation to rotate around one of the cartesian axes.
    class Rotate3d
    {
    public:
      Rotate3d(const Tensor<1, 3, double> &axis, const double angle)
        : rotation_matrix(
            Physics::Transformations::Rotations::rotation_matrix_3d(axis,
                                                                    angle))
      {}

      Point<3>
      operator()(const Point<3> &p) const
      {
        return static_cast<Point<3>>(rotation_matrix * p);
      }

    private:
      const Tensor<2, 3, double> rotation_matrix;
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
  shift(const Tensor<1, spacedim>    &shift_vector,
        Triangulation<dim, spacedim> &triangulation)
  {
    transform(internal::Shift<spacedim>(shift_vector), triangulation);
  }



  template <int dim, int spacedim>
  void
  rotate(const double /*angle*/,
         Triangulation<dim, spacedim> & /*triangulation*/)
  {
    AssertThrow(false,
                ExcMessage(
                  "GridTools::rotate() is only available for spacedim = 2."));
  }



  template <>
  void
  rotate(const double angle, Triangulation<1, 2> &triangulation)
  {
    transform(internal::Rotate2d(angle), triangulation);
  }



  template <>
  void
  rotate(const double angle, Triangulation<2, 2> &triangulation)
  {
    transform(internal::Rotate2d(angle), triangulation);
  }


  template <int dim>
  void
  rotate(const Tensor<1, 3, double> &axis,
         const double                angle,
         Triangulation<dim, 3>      &triangulation)
  {
    transform(internal::Rotate3d(axis, angle), triangulation);
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
    laplace_solve(const SparseMatrix<double>      &S,
                  const AffineConstraints<double> &constraints,
                  Vector<double>                  &u)
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
                    Triangulation<dim>                       &triangulation,
                    const Function<dim>                      *coefficient,
                    const bool solve_for_absolute_positions)
  {
    if (dim == 1)
      DEAL_II_NOT_IMPLEMENTED();

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

    const QGauss<dim> quadrature(4);

    Assert(triangulation.all_reference_cells_are_hyper_cube(),
           ExcNotImplemented());
    const auto reference_cell = ReferenceCells::get_hypercube<dim>();
    MatrixCreator::create_laplace_matrix(
      reference_cell.template get_default_linear_mapping<dim, dim>(),
      dof_handler,
      quadrature,
      S,
      coefficient);

    // set up the boundary values for the laplace problem
    std::array<AffineConstraints<double>, dim>                  constraints;
    typename std::map<unsigned int, Point<dim>>::const_iterator map_end =
      new_points.end();

    // Fill these maps using the data given by new_points
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        // Loop over all vertices of the cell and see if it is listed in the map
        // given as first argument of the function. We visit vertices multiple
        // times, so also check that if we have already added a constraint, we
        // don't do it a second time again.
        for (const unsigned int vertex_no : cell->vertex_indices())
          {
            const unsigned int vertex_index = cell->vertex_index(vertex_no);
            const Point<dim>  &vertex_point = cell->vertex(vertex_no);

            const typename std::map<unsigned int, Point<dim>>::const_iterator
              map_iter = new_points.find(vertex_index);

            if (map_iter != map_end)
              for (unsigned int i = 0; i < dim; ++i)
                if (constraints[i].is_constrained(
                      cell->vertex_dof_index(vertex_no, 0)) == false)
                  {
                    constraints[i].add_constraint(
                      cell->vertex_dof_index(vertex_no, 0),
                      {},
                      (solve_for_absolute_positions ?
                         map_iter->second[i] :
                         map_iter->second[i] - vertex_point[i]));
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
                v[i] = us[i](dof_index);
              else
                v[i] += us[i](dof_index);

            vertex_touched[cell->vertex_index(vertex_no)] = true;
          }
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
                    shift_vector[d] = uniform_distribution(rng);

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
                  shift_vector[d] = uniform_distribution(rng);

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
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  unsigned int find_closest_vertex(const MeshType<dim, spacedim> &mesh,
                                   const Point<spacedim>         &p,
                                   const std::vector<bool> &marked_vertices)
  {
    // first get the underlying triangulation from the mesh and determine
    // vertices and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const std::vector<Point<spacedim>> &vertices = tria.get_vertices();

    Assert(tria.get_vertices().size() == marked_vertices.size() ||
             marked_vertices.empty(),
           ExcDimensionMismatch(tria.get_vertices().size(),
                                marked_vertices.size()));

    // marked_vertices is expected to be a subset of used_vertices. Thus,
    // comparing the range marked_vertices.begin() to marked_vertices.end() with
    // the range used_vertices.begin() to used_vertices.end() the element in the
    // second range must be valid if the element in the first range is valid.
    Assert(
      marked_vertices.empty() ||
        std::equal(marked_vertices.begin(),
                   marked_vertices.end(),
                   tria.get_used_vertices().begin(),
                   [](bool p, bool q) { return !p || q; }),
      ExcMessage(
        "marked_vertices should be a subset of used vertices in the triangulation "
        "but marked_vertices contains one or more vertices that are not used vertices!"));

    // If marked_indices is empty, consider all used_vertices for finding the
    // closest vertex to the point. Otherwise, marked_indices is used.
    const std::vector<bool> &vertices_to_use =
      (marked_vertices.empty()) ? tria.get_used_vertices() : marked_vertices;

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
    for (unsigned int j = best_vertex + 1; j < vertices.size(); ++j)
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
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  unsigned int find_closest_vertex(const Mapping<dim, spacedim>  &mapping,
                                   const MeshType<dim, spacedim> &mesh,
                                   const Point<spacedim>         &p,
                                   const std::vector<bool> &marked_vertices)
  {
    // Take a shortcut in the simple case.
    if (mapping.preserves_vertex_locations() == true)
      return find_closest_vertex(mesh, p, marked_vertices);

    // first get the underlying triangulation from the mesh and determine
    // vertices and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    auto vertices = extract_used_vertices(tria, mapping);

    Assert(tria.get_vertices().size() == marked_vertices.size() ||
             marked_vertices.empty(),
           ExcDimensionMismatch(tria.get_vertices().size(),
                                marked_vertices.size()));

    // marked_vertices is expected to be a subset of used_vertices. Thus,
    // comparing the range marked_vertices.begin() to marked_vertices.end()
    // with the range used_vertices.begin() to used_vertices.end() the element
    // in the second range must be valid if the element in the first range is
    // valid.
    Assert(
      marked_vertices.empty() ||
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
      const Tensor<1, spacedim>              &point_direction,
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
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#endif
    find_active_cell_around_point(
      const Mapping<dim, spacedim>  &mapping,
      const MeshType<dim, spacedim> &mesh,
      const Point<spacedim>         &p,
      const std::vector<
        std::set<typename MeshType<dim, spacedim>::active_cell_iterator>>
        &vertex_to_cells,
      const std::vector<std::vector<Tensor<1, spacedim>>>
        &vertex_to_cell_centers,
      const typename MeshType<dim, spacedim>::active_cell_iterator &cell_hint,
      const std::vector<bool> &marked_vertices,
      const RTree<std::pair<Point<spacedim>, unsigned int>>
                  &used_vertices_rtree,
      const double tolerance,
      const RTree<
        std::pair<BoundingBox<spacedim>,
                  typename Triangulation<dim, spacedim>::active_cell_iterator>>
        *relevant_cell_bounding_boxes_rtree)
  {
    std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
              Point<dim>>
      cell_and_position;
    cell_and_position.first = mesh.end();

    // To handle points at the border we keep track of points which are close to
    // the unit cell:
    std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
              Point<dim>>
      cell_and_position_approx;

    if (relevant_cell_bounding_boxes_rtree != nullptr &&
        !relevant_cell_bounding_boxes_rtree->empty())
      {
        // create a bounding box around point p with 2*tolerance as side length.
        const auto bb = BoundingBox<spacedim>(p).create_extended(tolerance);

        if (relevant_cell_bounding_boxes_rtree->qbegin(
              boost::geometry::index::intersects(bb)) ==
            relevant_cell_bounding_boxes_rtree->qend())
          return cell_and_position;
      }

    bool found_cell  = false;
    bool approx_cell = false;

    unsigned int closest_vertex_index = 0;
    // ensure closest vertex index is a marked one, otherwise cell (with vertex
    // 0) might be found even though it is not marked. This is only relevant if
    // searching with rtree, using find_closest_vertex already can manage not
    // finding points
    if (marked_vertices.size() && !used_vertices_rtree.empty())
      {
        const auto itr =
          std::find(marked_vertices.begin(), marked_vertices.end(), true);
        Assert(itr != marked_vertices.end(),
               dealii::ExcMessage("No vertex has been marked!"));
        closest_vertex_index = std::distance(marked_vertices.begin(), itr);
      }

    Tensor<1, spacedim> vertex_to_point;
    auto                current_cell = cell_hint;

    // check whether cell has at least one marked vertex
    const auto cell_marked = [&mesh, &marked_vertices](const auto &cell) {
      if (marked_vertices.empty())
        return true;

      if (cell != mesh.active_cell_iterators().end())
        for (unsigned int i = 0; i < cell->n_vertices(); ++i)
          if (marked_vertices[cell->vertex_index(i)])
            return true;

      return false;
    };

    // check whether any cell in collection is marked
    const auto any_cell_marked = [&cell_marked](const auto &cells) {
      return std::any_of(cells.begin(),
                         cells.end(),
                         [&cell_marked](const auto &cell) {
                           return cell_marked(cell);
                         });
    };
    (void)any_cell_marked;

    while (found_cell == false)
      {
        // First look at the vertices of the cell cell_hint. If it's an
        // invalid cell, then query for the closest global vertex
        if (current_cell.state() == IteratorState::valid &&
            cell_marked(cell_hint))
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
            // For some clang-based compilers and boost versions the call to
            // RTree::query doesn't compile. Since using an rtree here is just a
            // performance improvement disabling this branch is OK.
            // This is fixed in boost in
            // https://github.com/boostorg/numeric_conversion/commit/50a1eae942effb0a9b90724323ef8f2a67e7984a
#if defined(DEAL_II_WITH_BOOST_BUNDLED) ||                \
  !(defined(__clang_major__) && __clang_major__ >= 16) || \
  BOOST_VERSION >= 108100
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

                // Searching for a point which is located outside the
                // triangulation results in res.size() = 0
                Assert(res.size() < 2,
                       dealii::ExcMessage("There can not be multiple results"));

                if (res.size() > 0)
                  if (any_cell_marked(vertex_to_cells[res[0].second]))
                    closest_vertex_index = res[0].second;
              }
            else
#endif
              {
                closest_vertex_index = GridTools::find_closest_vertex(
                  mapping, mesh, p, marked_vertices);
              }
            vertex_to_point = p - mesh.get_vertices()[closest_vertex_index];
          }

        if constexpr (running_in_debug_mode())
          {
            {
              // Double-check if found index is at marked cell
              Assert(any_cell_marked(vertex_to_cells[closest_vertex_index]),
                     dealii::ExcMessage("Found non-marked vertex"));
            }
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
        // hint. Most likely we will find the point in them.
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
                    if ((*cell)->reference_cell().contains_point(p_unit,
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
                    const double dist = p_unit.distance(
                      (*cell)->reference_cell().closest_point(p_unit));
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
    const Point<spacedim>                                             &position,
    const Mapping<dim, spacedim>                                      &mapping)
  {
    const auto         vertices         = mapping.get_vertices(cell);
    double             minimum_distance = position.distance_square(vertices[0]);
    unsigned int       closest_vertex   = 0;
    const unsigned int n_vertices       = cell->n_vertices();

    for (unsigned int v = 1; v < n_vertices; ++v)
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
      template <typename MeshType>
      DEAL_II_CXX20_REQUIRES(
        concepts::is_triangulation_or_dof_handler<MeshType>)
      std::tuple<
        BoundingBox<MeshType::space_dimension>,
        bool> compute_cell_predicate_bounding_box(const typename MeshType::
                                                    cell_iterator &parent_cell,
                                                  const std::function<bool(
                                                    const typename MeshType::
                                                      active_cell_iterator &)>
                                                    &predicate)
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
        if (active_cells.empty() || i == active_cells.size())
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



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<BoundingBox<MeshType::space_dimension>> compute_mesh_predicate_bounding_box(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::active_cell_iterator &)>
                        &predicate,
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
            volumes.reserve(merged_b_boxes.size());
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
    const std::vector<Point<spacedim>>                    &points)
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
    const std::vector<Point<spacedim>>                          &points)
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
  std::map<unsigned int, types::global_vertex_index>
  compute_local_to_global_vertex_index_map(
    const Triangulation<dim, spacedim> &triangulation)
  {
    std::map<unsigned int, types::global_vertex_index>
      local_to_global_vertex_index;

#ifndef DEAL_II_WITH_MPI

    // If we don't have MPI then all vertices are local
    for (unsigned int i = 0; i < triangulation.n_vertices(); ++i)
      local_to_global_vertex_index[i] = i;

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
    std::set<active_cell_iterator> missing_vert_cells;
    std::set<unsigned int>         used_vertex_index;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            for (const unsigned int i : cell->vertex_indices())
              {
                types::subdomain_id lowest_subdomain_id = cell->subdomain_id();
                for (const auto &adjacent_cell :
                     vertex_to_cell[cell->vertex_index(i)])
                  lowest_subdomain_id = std::min(lowest_subdomain_id,
                                                 adjacent_cell->subdomain_id());

                // See if this process "owns" this vertex
                if (lowest_subdomain_id == cell->subdomain_id())
                  {
                    // Check that the vertex we are working on is a vertex that
                    // has not been dealt with yet
                    if (used_vertex_index.find(cell->vertex_index(i)) ==
                        used_vertex_index.end())
                      {
                        // Set the local index
                        local_to_global_vertex_index[cell->vertex_index(i)] =
                          next_index++;

                        // Store the information that will be sent to the
                        // adjacent cells on other subdomains
                        for (const auto &adjacent_cell :
                             vertex_to_cell[cell->vertex_index(i)])
                          if (adjacent_cell->subdomain_id() !=
                              cell->subdomain_id())
                            {
                              std::pair<types::subdomain_id,
                                        types::global_vertex_index>
                                tmp(adjacent_cell->subdomain_id(),
                                    cell->vertex_index(i));
                              if (vertices_added.find(tmp) ==
                                  vertices_added.end())
                                {
                                  vertices_to_send[adjacent_cell
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
            for (const unsigned int i : cell->face_indices())
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
    max_cellid_size = Utilities::MPI::max(max_cellid_size,
                                          triangulation.get_mpi_communicator());

    // Make indices global by getting the number of vertices owned by each
    // processors and shifting the indices accordingly
    types::global_vertex_index shift = 0;
    int                        ierr  = MPI_Exscan(
      &next_index,
      &shift,
      1,
      Utilities::MPI::mpi_type_id_for_type<types::global_vertex_index>,
      MPI_SUM,
      triangulation.get_mpi_communicator());
    AssertThrowMPI(ierr);

    for (auto &global_index_it : local_to_global_vertex_index)
      global_index_it.second += shift;


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
        ierr = MPI_Isend(
          vertices_send_buffers[i].data(),
          buffer_size,
          Utilities::MPI::mpi_type_id_for_type<types::global_vertex_index>,
          destination,
          mpi_tag,
          triangulation.get_mpi_communicator(),
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
        ierr = MPI_Recv(
          vertices_recv_buffers[i].data(),
          buffer_size,
          Utilities::MPI::mpi_type_id_for_type<types::global_vertex_index>,
          source,
          mpi_tag,
          triangulation.get_mpi_communicator(),
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
      }

    // At this point, wait for all of the isend operations to finish:
    MPI_Waitall(first_requests.size(),
                first_requests.data(),
                MPI_STATUSES_IGNORE);


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
                         triangulation.get_mpi_communicator(),
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
                        triangulation.get_mpi_communicator(),
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

    // At this point, wait for all of the isend operations of the second round
    // to finish:
    MPI_Waitall(second_requests.size(),
                second_requests.data(),
                MPI_STATUSES_IGNORE);
#endif

    return local_to_global_vertex_index;
  }



  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          Triangulation<dim, spacedim>    &triangulation,
                          const SparsityTools::Partitioner partitioner)
  {
    Assert((dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              &triangulation) == nullptr),
           ExcMessage("Objects of type parallel::distributed::Triangulation "
                      "are already partitioned implicitly and can not be "
                      "partitioned again explicitly."));

    std::vector<unsigned int> cell_weights;

    // Get cell weighting if a signal has been attached to the triangulation
    if (!triangulation.signals.weight.empty())
      {
        cell_weights.resize(triangulation.n_active_cells(), 0U);

        // In a first step, obtain the weights of the locally owned
        // cells. For all others, the weight remains at the zero the
        // vector was initialized with above.
        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            cell_weights[cell->active_cell_index()] =
              triangulation.signals.weight(cell, CellStatus::cell_will_persist);

        // If this is a parallel triangulation, we then need to also
        // get the weights for all other cells. We have asserted above
        // that this function can't be used for
        // parallel::distributed::Triangulation objects, so the only
        // ones we have to worry about here are
        // parallel::shared::Triangulation
        if (const auto shared_tria =
              dynamic_cast<parallel::shared::Triangulation<dim, spacedim> *>(
                &triangulation))
          Utilities::MPI::sum(cell_weights,
                              shared_tria->get_mpi_communicator(),
                              cell_weights);

        // verify that the global sum of weights is larger than 0
        Assert(std::accumulate(cell_weights.begin(),
                               cell_weights.end(),
                               std::uint64_t(0)) > 0,
               ExcMessage("The global sum of weights over all active cells "
                          "is zero. Please verify how you generate weights."));
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
                          Triangulation<dim, spacedim>    &triangulation,
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
                          const SparsityPattern        &cell_connection_graph,
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
    if (!triangulation.signals.weight.empty())
      {
        cell_weights.resize(triangulation.n_active_cells(), 0U);

        // In a first step, obtain the weights of the locally owned
        // cells. For all others, the weight remains at the zero the
        // vector was initialized with above.
        for (const auto &cell : triangulation.active_cell_iterators() |
                                  IteratorFilters::LocallyOwnedCell())
          cell_weights[cell->active_cell_index()] =
            triangulation.signals.weight(cell, CellStatus::cell_will_persist);

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
                              shared_tria->get_mpi_communicator(),
                              cell_weights);

        // verify that the global sum of weights is larger than 0
        Assert(std::accumulate(cell_weights.begin(),
                               cell_weights.end(),
                               std::uint64_t(0)) > 0,
               ExcMessage("The global sum of weights over all active cells "
                          "is zero. Please verify how you generate weights."));
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
                          const SparsityPattern        &cell_connection_graph,
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
                                           unsigned int      &current_proc_idx,
                                           unsigned int      &current_cell_idx,
                                           const unsigned int n_active_cells,
                                           const unsigned int n_partitions)
    {
      if (cell->is_active())
        {
          while (current_cell_idx >=
                 std::floor(static_cast<std::uint_least64_t>(n_active_cells) *
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
    Assert(triangulation.signals.weight.empty(), ExcNotImplemented());

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
        for (const auto &cell : triangulation.cell_iterators_on_level(lvl))
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

  namespace internal
  {
    namespace
    {
      // Split get_subdomain_association() for p::d::T since we want to compile
      // it in 1d but none of the p4est stuff is available in 1d.
      template <int dim, int spacedim>
      void
      get_subdomain_association(
        const parallel::distributed::Triangulation<dim, spacedim>
                                                          &triangulation,
        const std::vector<CellId>                         &cell_ids,
        [[maybe_unused]] std::vector<types::subdomain_id> &subdomain_ids)
      {
#ifndef DEAL_II_WITH_P4EST
        (void)triangulation;
        (void)cell_ids;
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
                             *>(triangulation.get_p4est()),
                cell_id.get_coarse_cell_id(),
                &p4est_cell,
                Utilities::MPI::this_mpi_process(
                  triangulation.get_mpi_communicator()));

            Assert(owner >= 0, ExcMessage("p4est should know the owner."));

            subdomain_ids.push_back(owner);
          }
#endif
      }



      template <int spacedim>
      void
      get_subdomain_association(
        const parallel::distributed::Triangulation<1, spacedim> &,
        const std::vector<CellId> &,
        std::vector<types::subdomain_id> &)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }
    } // anonymous namespace
  }   // namespace internal



  template <int dim, int spacedim>
  std::vector<types::subdomain_id>
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            const std::vector<CellId>          &cell_ids)
  {
    std::vector<types::subdomain_id> subdomain_ids;
    subdomain_ids.reserve(cell_ids.size());

    if (dynamic_cast<
          const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
          &triangulation) != nullptr)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }
    else if (const parallel::distributed::Triangulation<dim, spacedim>
               *parallel_tria = dynamic_cast<
                 const parallel::distributed::Triangulation<dim, spacedim> *>(
                 &triangulation))
      {
        internal::get_subdomain_association(*parallel_tria,
                                            cell_ids,
                                            subdomain_ids);
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
                            std::vector<types::subdomain_id>   &subdomain)
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
      objective_function(const Iterator        &object,
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
            objective += (child_alternating_forms[c][i] -
                          average_parent_alternating_form /
                            Utilities::fixed_power<structdim>(2))
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
      get_face_midpoint(const Iterator    &object,
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
      get_face_midpoint(const Iterator    &object,
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
      get_face_midpoint(const Iterator    &object,
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
        const typename Triangulation<dim, spacedim>::cell_iterator &cell =
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
    Triangulation<dim, spacedim>          &tria,
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
    Triangulation<dim, spacedim>             &tria,
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

    AssertThrow(tria.all_reference_cells_are_hyper_cube(), ExcNotImplemented());

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
                ++boundary_face_counter;
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
            ++refinement_cycles;
          }
      }

    if (has_cells_with_dim_faces_on_boundary)
      {
        tria.refine_global(1);
        ++refinement_cycles;
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
            DEAL_II_NOT_IMPLEMENTED();
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
                DEAL_II_NOT_IMPLEMENTED();
              }
          }
      }

    // if no cells need to be added, then no regularization is necessary.
    // Restore things as they were before this function was called.
    if (cells_to_add.empty())
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
            CellData<dim> c(cell->n_vertices());
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
              CellData<2> quad(face->n_vertices());
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
    const Cache<dim, spacedim>         &cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint)
  {
    const auto cqmp = compute_point_locations_try_all(cache, points, cell_hint);
    // Splitting the tuple's components
    auto &cells   = std::get<0>(cqmp);
    auto &qpoints = std::get<1>(cqmp);
    auto &maps    = std::get<2>(cqmp);

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
    const Cache<dim, spacedim>         &cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint)
  {
    Assert((dim == spacedim),
           ExcMessage("Only implemented for dim==spacedim."));

    // Alias
    namespace bgi = boost::geometry::index;

    // Get the mapping
    const auto &mapping = cache.get_mapping();

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

    // Now make a tree of indices for the points
    // [TODO] This would work better with pack_rtree_of_indices, but
    // windows does not like it. Build a tree with pairs of point and id
    std::vector<std::pair<Point<spacedim>, unsigned int>> points_and_ids(np);
    for (unsigned int i = 0; i < np; ++i)
      points_and_ids[i] = std::make_pair(points[i], i);
    const auto p_tree = pack_rtree(points_and_ids);

    // Keep track of all found points
    std::vector<bool> found_points(points.size(), false);

    // Check if a point was found
    const auto already_found = [&found_points](const auto &id) {
      AssertIndexRange(id.second, found_points.size());
      return found_points[id.second];
    };

    // check if the given cell was already in the vector of cells before. If so,
    // insert in the corresponding vectors the reference point and the id.
    // Otherwise append a new entry to all vectors.
    const auto store_cell_point_and_id =
      [&](
        const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
        const Point<dim>   &ref_point,
        const unsigned int &id) {
        const auto it = std::find(cells_out.rbegin(), cells_out.rend(), cell);
        if (it != cells_out.rend())
          {
            const auto cell_id =
              (cells_out.size() - 1 - (it - cells_out.rbegin()));
            qpoints_out[cell_id].emplace_back(ref_point);
            maps_out[cell_id].emplace_back(id);
          }
        else
          {
            cells_out.emplace_back(cell);
            qpoints_out.emplace_back(std::vector<Point<dim>>({ref_point}));
            maps_out.emplace_back(std::vector<unsigned int>({id}));
          }
      };

    // Check all points within a given pair of box and cell
    const auto check_all_points_within_box = [&](const auto &leaf) {
      const double                relative_tolerance = 1e-12;
      const BoundingBox<spacedim> box =
        leaf.first.create_extended_relative(relative_tolerance);
      const auto &cell_hint = leaf.second;

      for (const auto &point_and_id :
           p_tree | bgi::adaptors::queried(!bgi::satisfies(already_found) &&
                                           bgi::intersects(box)))
        {
          const auto id = point_and_id.second;
          const auto cell_and_ref =
            GridTools::find_active_cell_around_point(cache,
                                                     points[id],
                                                     cell_hint);
          const auto &cell      = cell_and_ref.first;
          const auto &ref_point = cell_and_ref.second;

          if (cell.state() == IteratorState::valid)
            store_cell_point_and_id(cell, ref_point, id);
          else
            missing_points_out.emplace_back(id);

          // Don't look anymore for this point
          found_points[id] = true;
        }
    };

    // If a hint cell was given, use it
    if (cell_hint.state() == IteratorState::valid)
      check_all_points_within_box(
        std::make_pair(mapping.get_bounding_box(cell_hint), cell_hint));

    // Now loop over all points that have not been found yet
    for (unsigned int i = 0; i < np; ++i)
      if (found_points[i] == false)
        {
          // Get the closest cell to this point
          const auto leaf = b_tree.qbegin(bgi::nearest(points[i], 1));
          // Now checks all points that fall within this box
          if (leaf != b_tree.qend())
            check_all_points_within_box(*leaf);
          else
            {
              // We should not get here. Throw an error.
              DEAL_II_ASSERT_UNREACHABLE();
            }
        }
    // Now make sure we send out the rest of the points that we did not find.
    for (unsigned int i = 0; i < np; ++i)
      if (found_points[i] == false)
        missing_points_out.emplace_back(i);

    // Debug Checking
    AssertDimension(cells_out.size(), maps_out.size());
    AssertDimension(cells_out.size(), qpoints_out.size());

    if constexpr (running_in_debug_mode())
      {
        unsigned int c   = cells_out.size();
        unsigned int qps = 0;
        // The number of points in all
        // the cells must be the same as
        // the number of points we
        // started off from,
        // plus the points which were ignored
        for (unsigned int n = 0; n < c; ++n)
          {
            AssertDimension(qpoints_out[n].size(), maps_out[n].size());
            qps += qpoints_out[n].size();
          }

        Assert(qps + missing_points_out.size() == np,
               ExcDimensionMismatch(qps + missing_points_out.size(), np));
      }

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
    const GridTools::Cache<dim, spacedim>                 &cache,
    const std::vector<Point<spacedim>>                    &points,
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const double                                           tolerance,
    const std::vector<bool>                               &marked_vertices,
    const bool enforce_unique_mapping)
  {
    // run internal function ...
    const auto all =
      internal::distributed_compute_point_locations(cache,
                                                    points,
                                                    global_bboxes,
                                                    marked_vertices,
                                                    tolerance,
                                                    false,
                                                    enforce_unique_mapping)
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
    /**
     * Determine for each rank which entry of @p entities it
     * might own. The first entry of the returned tuple is a list of
     * ranks and the second and third entry give CRS data
     * structure (pointers within a list of indices).
     */
    template <int spacedim, typename T>
    std::tuple<std::vector<unsigned int>,
               std::vector<unsigned int>,
               std::vector<unsigned int>>
    guess_owners_of_entities(
      const MPI_Comm                                         comm,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const std::vector<T>                                  &entities,
      const double                                           tolerance)
    {
      std::vector<std::pair<unsigned int, unsigned int>> ranks_and_indices;
      ranks_and_indices.reserve(entities.size());

#if defined(DEAL_II_WITH_ARBORX)
      static constexpr bool use_arborx = true;
#else
      static constexpr bool use_arborx = false;
#endif
      // Lambda to process bboxes if global_bboxes.size()>1 or ArborX not
      // available
      const auto process_bboxes = [&]() -> void {
        std::vector<std::vector<BoundingBox<spacedim>>> global_bboxes_temp;
        auto *global_bboxes_to_be_used = &global_bboxes;

        if (global_bboxes.size() == 1 && use_arborx == false)
          {
            global_bboxes_temp =
              Utilities::MPI::all_gather(comm, global_bboxes[0]);
            global_bboxes_to_be_used = &global_bboxes_temp;
          }

        // helper function to determine if a bounding box is valid
        const auto is_valid = [](const auto &bb) {
          for (unsigned int i = 0; i < spacedim; ++i)
            if (bb.get_boundary_points().first[i] >
                bb.get_boundary_points().second[i])
              return false;

          return true;
        };

        // linearize vector of vectors
        std::vector<std::pair<BoundingBox<spacedim>, unsigned int>>
          boxes_and_ranks;

        for (unsigned rank = 0; rank < global_bboxes_to_be_used->size(); ++rank)
          for (const auto &box : (*global_bboxes_to_be_used)[rank])
            if (is_valid(box))
              boxes_and_ranks.emplace_back(box, rank);

        // pack boxes into r-tree
        const auto tree = pack_rtree(boxes_and_ranks);

        // loop over all entities
        for (unsigned int i = 0; i < entities.size(); ++i)
          {
            // create a bounding box with tolerance
            const auto bb =
              BoundingBox<spacedim>(entities[i]).create_extended(tolerance);

            // determine ranks potentially owning point/bounding box
            std::set<unsigned int> my_ranks;

            for (const auto &box_and_rank :
                 tree | boost::geometry::index::adaptors::queried(
                          boost::geometry::index::intersects(bb)))
              my_ranks.insert(box_and_rank.second);

            for (const auto rank : my_ranks)
              ranks_and_indices.emplace_back(rank, i);
          }
      };

      if constexpr (use_arborx)
        {
          if (global_bboxes.size() == 1)
            {
              ArborXWrappers::DistributedTree distributed_tree(
                comm, global_bboxes[0]);
              std::vector<BoundingBox<spacedim>> query_bounding_boxes;
              query_bounding_boxes.reserve(entities.size());
              for (const auto &entity : entities)
                query_bounding_boxes.emplace_back(
                  BoundingBox<spacedim>(entity).create_extended(tolerance));

              ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
                query_bounding_boxes);
              const auto &[indices_ranks, offsets] =
                distributed_tree.query(bb_intersect);

              for (unsigned long int i = 0; i < offsets.size() - 1; ++i)
                {
                  std::set<unsigned int> my_ranks;
                  for (int j = offsets[i]; j < offsets[i + 1]; ++j)
                    my_ranks.insert(indices_ranks[j].second);

                  for (const auto rank : my_ranks)
                    ranks_and_indices.emplace_back(rank, i);
                }
            }
          else
            {
              // global_bboxes.size()>1
              process_bboxes();
            }
        }
      else
        {
          // No ArborX
          process_bboxes();
        }


      // convert to CRS
      std::sort(ranks_and_indices.begin(), ranks_and_indices.end());

      std::vector<unsigned int> ranks;
      std::vector<unsigned int> ptr;
      std::vector<unsigned int> indices;

      unsigned int current_rank = numbers::invalid_unsigned_int;

      for (const std::pair<unsigned int, unsigned int> &i : ranks_and_indices)
        {
          if (current_rank != i.first)
            {
              current_rank = i.first;
              ranks.push_back(current_rank);
              ptr.push_back(indices.size());
            }

          indices.push_back(i.second);
        }
      ptr.push_back(indices.size());

      return {std::move(ranks), std::move(ptr), std::move(indices)};
    }



    template <int dim, int spacedim>
    std::vector<
      std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                Point<dim>>>
    find_all_locally_owned_active_cells_around_point(
      const Cache<dim, spacedim>                                  &cache,
      const Point<spacedim>                                       &point,
      typename Triangulation<dim, spacedim>::active_cell_iterator &cell_hint,
      const std::vector<bool> &marked_vertices,
      const double             tolerance,
      const bool               enforce_unique_mapping)
    {
      std::vector<
        std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
                  Point<dim>>>
        locally_owned_active_cells_around_point;

      const auto first_cell = GridTools::find_active_cell_around_point(
        cache.get_mapping(),
        cache.get_triangulation(),
        point,
        cache.get_vertex_to_cell_map(),
        cache.get_vertex_to_cell_centers_directions(),
        cell_hint,
        marked_vertices,
        cache.get_used_vertices_rtree(),
        tolerance,
        &cache.get_locally_owned_cell_bounding_boxes_rtree());

      const unsigned int my_rank = Utilities::MPI::this_mpi_process(
        cache.get_triangulation().get_mpi_communicator());

      cell_hint = first_cell.first;
      if (cell_hint.state() == IteratorState::valid)
        {
          const auto active_cells_around_point =
            GridTools::find_all_active_cells_around_point(
              cache.get_mapping(),
              cache.get_triangulation(),
              point,
              tolerance,
              first_cell,
              &cache.get_vertex_to_cell_map());

          if (enforce_unique_mapping)
            {
              // check if the rank of this process is the lowest of all cells
              // if not, the other process will handle this cell and we don't
              // have to do here anything in the case of unique mapping
              unsigned int lowes_rank = numbers::invalid_unsigned_int;

              for (const auto &cell : active_cells_around_point)
                lowes_rank = std::min(lowes_rank, cell.first->subdomain_id());

              if (lowes_rank != my_rank)
                return {};
            }

          locally_owned_active_cells_around_point.reserve(
            active_cells_around_point.size());

          for (const auto &cell : active_cells_around_point)
            if (cell.first->is_locally_owned())
              locally_owned_active_cells_around_point.push_back(cell);
        }

      std::sort(locally_owned_active_cells_around_point.begin(),
                locally_owned_active_cells_around_point.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });

      if (enforce_unique_mapping &&
          locally_owned_active_cells_around_point.size() > 1)
        // in the case of unique mapping, we only need a single cell
        return {locally_owned_active_cells_around_point.front()};
      else
        return locally_owned_active_cells_around_point;
    }

    template <int dim, int spacedim>
    DistributedComputePointLocationsInternal<dim, spacedim>::
      DistributedComputePointLocationsInternal()
      : n_searched_points(numbers::invalid_unsigned_int)
    {}

    template <int dim, int spacedim>
    void
    DistributedComputePointLocationsInternal<dim, spacedim>::finalize_setup()
    {
      // before reshuffeling the data check if data.recv_components and
      // n_searched_points are in a valid state.
      Assert(n_searched_points != numbers::invalid_unsigned_int,
             ExcInternalError());
      Assert(recv_components.empty() ||
               std::get<1>(*std::max_element(recv_components.begin(),
                                             recv_components.end(),
                                             [](const auto &a, const auto &b) {
                                               return std::get<1>(a) <
                                                      std::get<1>(b);
                                             })) < n_searched_points,
             ExcInternalError());

      send_ranks.clear();
      recv_ranks.clear();
      send_ptrs.clear();
      recv_ptrs.clear();

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

      if (recv_components.size() > 0)
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
    }



    template <int dim, int spacedim>
    DistributedComputePointLocationsInternal<dim, spacedim>
    distributed_compute_point_locations(
      const GridTools::Cache<dim, spacedim>                 &cache,
      const std::vector<Point<spacedim>>                    &points,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const std::vector<bool>                               &marked_vertices,
      const double                                           tolerance,
      const bool                                             perform_handshake,
      const bool enforce_unique_mapping)
    {
      DistributedComputePointLocationsInternal<dim, spacedim> result;
      result.n_searched_points = points.size();

      auto &send_components = result.send_components;
      auto &recv_components = result.recv_components;

      const auto comm = cache.get_triangulation().get_mpi_communicator();

      const auto potential_owners = internal::guess_owners_of_entities(
        comm, global_bboxes, points, tolerance);

      const auto &potential_owners_ranks   = std::get<0>(potential_owners);
      const auto &potential_owners_ptrs    = std::get<1>(potential_owners);
      const auto &potential_owners_indices = std::get<2>(potential_owners);

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

      Assert(
        (marked_vertices.empty()) ||
          (marked_vertices.size() == cache.get_triangulation().n_vertices()),
        ExcMessage(
          "The marked_vertices vector has to be either empty or its size has "
          "to equal the number of vertices of the triangulation."));

      using RequestType = std::vector<std::pair<unsigned int, Point<spacedim>>>;
      using AnswerType  = std::vector<unsigned int>;

      // In the case that a marked_vertices vector has been given and none
      // of its entries is true, we know that this process does not own
      // any of the incoming points (and it will not send any data) so
      // that we can take a short cut.
      const bool has_relevant_vertices =
        (marked_vertices.empty()) ||
        (std::find(marked_vertices.begin(), marked_vertices.end(), true) !=
         marked_vertices.end());

      const auto create_request = [&](const unsigned int other_rank) {
        const auto other_rank_index = translate(other_rank);

        RequestType request;
        request.reserve(potential_owners_ptrs[other_rank_index + 1] -
                        potential_owners_ptrs[other_rank_index]);

        for (unsigned int i = potential_owners_ptrs[other_rank_index];
             i < potential_owners_ptrs[other_rank_index + 1];
             ++i)
          request.emplace_back(potential_owners_indices[i],
                               points[potential_owners_indices[i]]);

        return request;
      };

      const auto answer_request =
        [&](const unsigned int &other_rank,
            const RequestType  &request) -> AnswerType {
        AnswerType answer(request.size(), 0);

        if (has_relevant_vertices)
          {
            cell_hint = cache.get_triangulation().begin_active();

            for (unsigned int i = 0; i < request.size(); ++i)
              {
                const auto &index_and_point = request[i];

                const auto cells_and_reference_positions =
                  find_all_locally_owned_active_cells_around_point(
                    cache,
                    index_and_point.second,
                    cell_hint,
                    marked_vertices,
                    tolerance,
                    enforce_unique_mapping);

                if (cell_hint.state() != IteratorState::valid)
                  cell_hint = cache.get_triangulation().begin_active();

                for (const auto &cell_and_reference_position :
                     cells_and_reference_positions)
                  {
                    const auto cell = cell_and_reference_position.first;
                    auto       reference_position =
                      cell_and_reference_position.second;

                    reference_position =
                      cell->reference_cell().closest_point(reference_position);

                    send_components.emplace_back(
                      std::pair<int, int>(cell->level(), cell->index()),
                      other_rank,
                      index_and_point.first,
                      reference_position,
                      index_and_point.second,
                      numbers::invalid_unsigned_int);
                  }

                answer[i] = cells_and_reference_positions.size();
              }
          }

        if (perform_handshake)
          return answer;
        else
          return {};
      };

      const auto process_answer = [&](const unsigned int other_rank,
                                      const AnswerType  &answer) {
        if (perform_handshake)
          {
            const auto other_rank_index = translate(other_rank);

            for (unsigned int i = 0; i < answer.size(); ++i)
              for (unsigned int j = 0; j < answer[i]; ++j)
                recv_components.emplace_back(
                  other_rank,
                  potential_owners_indices
                    [i + potential_owners_ptrs[other_rank_index]],
                  numbers::invalid_unsigned_int);
          }
      };

      Utilities::MPI::ConsensusAlgorithms::selector<RequestType, AnswerType>(
        potential_owners_ranks,
        create_request,
        answer_request,
        process_answer,
        comm);

      result.finalize_setup();

      return result;
    }



    template <int structdim, int spacedim>
    template <int dim>
    DistributedComputePointLocationsInternal<dim, spacedim>
    DistributedComputeIntersectionLocationsInternal<structdim, spacedim>::
      convert_to_distributed_compute_point_locations_internal(
        const unsigned int                  n_points_1D,
        const Triangulation<dim, spacedim> &tria,
        const Mapping<dim, spacedim>       &mapping,
        std::vector<Quadrature<spacedim>>  *mapped_quadratures_recv_comp,
        const bool consistent_numbering_of_sender_and_receiver) const
    {
      using CellIterator =
        typename Triangulation<dim, spacedim>::active_cell_iterator;

      if (mapped_quadratures_recv_comp != nullptr)
        {
          AssertDimension(mapped_quadratures_recv_comp->size(), 0);
          mapped_quadratures_recv_comp->reserve(recv_components.size());
        }

      GridTools::internal::DistributedComputePointLocationsInternal<dim,
                                                                    spacedim>
        result;

      // We need quadrature rules for the intersections. We are using a
      // QGaussSimplex quadrature rule since CGAL always returns simplices
      // as intersections.
      const QGaussSimplex<structdim> quadrature(n_points_1D);

      // Resulting quadrature points get different indices. In the case the
      // requested intersections are unique also the resulting quadrature
      // points are unique and we can simply number the points in an
      // ascending way.
      for (const auto &recv_component : recv_components)
        {
          // dependent on the size of the intersection an empty quadrature
          // is returned. Therefore, we have to compute the quadrature also
          // here.
          const Quadrature<spacedim> &quad =
            quadrature.compute_affine_transformation(
              std::get<2>(recv_component));

          for (unsigned int i = 0; i < quad.size(); ++i)
            {
              // the third component of result.recv_components is not needed
              // before finalize_setup.
              result.recv_components.emplace_back(
                std::get<0>(recv_component),
                result.recv_components.size(), // number of point
                numbers::invalid_unsigned_int);
            }

          // append quadrature
          if (mapped_quadratures_recv_comp != nullptr)
            mapped_quadratures_recv_comp->push_back(quad);
        }

      // since empty quadratures might be present we have to set the number
      // of searched points after inserting the point indices into
      // recv_components
      result.n_searched_points = result.recv_components.size();

      // send_ranks, counter, and indices_of_rank is only needed if
      // consistent_numbering_of_sender_and_receiver==true
      // indices_of_rank is always empty if deal.II is compiled without MPI
      std::map<unsigned int, std::vector<unsigned int>> indices_of_rank;
      std::map<unsigned int, unsigned int>              counter;
      std::set<unsigned int>                            send_ranks;
      if (consistent_numbering_of_sender_and_receiver)
        {
          for (const auto &sc : send_components)
            send_ranks.insert(std::get<1>(sc));

          for (const auto rank : send_ranks)
            counter[rank] = 0;

          // indices assigned at recv side needed to fill send_components
          indices_of_rank = communicate_indices(result.recv_components,
                                                tria.get_mpi_communicator());
        }

      for (const auto &send_component : send_components)
        {
          const CellIterator cell(&tria,
                                  std::get<0>(send_component).first,
                                  std::get<0>(send_component).second);

          const Quadrature<spacedim> &quad =
            quadrature.compute_affine_transformation(
              std::get<3>(send_component));

          const auto rank = std::get<1>(send_component);

          for (unsigned int q = 0; q < quad.size(); ++q)
            {
              // the fifth component of result.send_components is filled
              // during sorting the data and initializing the CRS structures
              result.send_components.emplace_back(std::make_tuple(
                std::get<0>(send_component),
                rank,
                indices_of_rank.empty() ?
                  result.send_components.size() :
                  indices_of_rank.at(rank)[counter.at(rank)],
                mapping.transform_real_to_unit_cell(cell, quad.point(q)),
                quad.point(q),
                numbers::invalid_unsigned_int));

              if (!indices_of_rank.empty())
                ++counter[rank];
            }
        }

      result.finalize_setup();

      return result;
    }



    template <int structdim, int spacedim>
    std::map<unsigned int, std::vector<unsigned int>>
    DistributedComputeIntersectionLocationsInternal<structdim, spacedim>::
      communicate_indices(
        [[maybe_unused]] const std::vector<
          std::tuple<unsigned int, unsigned int, unsigned int>>
                                       &point_recv_components,
        [[maybe_unused]] const MPI_Comm comm) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      return {};
#else
      // since we are converting to DistributedComputePointLocationsInternal
      // we use the RPE tag
      const auto mpi_tag =
        Utilities::MPI::internal::Tags::remote_point_evaluation;

      const unsigned int my_rank = Utilities::MPI::this_mpi_process(comm);

      std::set<unsigned int> send_ranks;
      for (const auto &sc : send_components)
        send_ranks.insert(std::get<1>(sc));
      std::set<unsigned int> recv_ranks;
      for (const auto &rc : recv_components)
        recv_ranks.insert(std::get<0>(rc));

      std::vector<MPI_Request> requests;
      requests.reserve(send_ranks.size());

      // rank to used indices on the rank needed on sending side
      std::map<unsigned int, std::vector<unsigned int>> indices_of_rank;
      indices_of_rank[my_rank] = std::vector<unsigned int>();

      // rank to used indices on the rank known on recv side
      std::map<unsigned int, std::vector<unsigned int>> send_indices_of_rank;
      for (const auto rank : recv_ranks)
        if (rank != my_rank)
          send_indices_of_rank[rank] = std::vector<unsigned int>();

      // fill the maps
      for (const auto &point_recv_component : point_recv_components)
        {
          const auto rank = std::get<0>(point_recv_component);
          const auto idx  = std::get<1>(point_recv_component);

          if (rank == my_rank)
            indices_of_rank[rank].emplace_back(idx);
          else
            send_indices_of_rank[rank].emplace_back(idx);
        }

      // send indices to the ranks we normally receive from
      for (const auto rank : recv_ranks)
        {
          if (rank == my_rank)
            continue;

          auto buffer = Utilities::pack(send_indices_of_rank[rank], false);

          requests.push_back(MPI_Request());

          const int ierr = MPI_Isend(buffer.data(),
                                     buffer.size(),
                                     MPI_CHAR,
                                     rank,
                                     mpi_tag,
                                     comm,
                                     &requests.back());
          AssertThrowMPI(ierr);
        }

      // receive indices at the ranks we normally send from
      for (const auto rank : send_ranks)
        {
          if (rank == my_rank)
            continue;

          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
          AssertThrowMPI(ierr);

          int message_length;
          ierr = MPI_Get_count(&status, MPI_CHAR, &message_length);
          AssertThrowMPI(ierr);

          std::vector<char> buffer(message_length);

          ierr = MPI_Recv(buffer.data(),
                          buffer.size(),
                          MPI_CHAR,
                          status.MPI_SOURCE,
                          mpi_tag,
                          comm,
                          MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          indices_of_rank[status.MPI_SOURCE] =
            Utilities::unpack<std::vector<unsigned int>>(buffer, false);
        }

      // make sure all messages have been sent
      const int ierr =
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      return indices_of_rank;
#endif
    }



    template <int structdim, int dim, int spacedim>
    DistributedComputeIntersectionLocationsInternal<structdim, spacedim>
    distributed_compute_intersection_locations(
      const Cache<dim, spacedim>                      &cache,
      const std::vector<std::vector<Point<spacedim>>> &intersection_requests,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const std::vector<bool>                               &marked_vertices,
      const double                                           tolerance)
    {
      using IntersectionRequest = std::vector<Point<spacedim>>;

      using IntersectionAnswer =
        typename DistributedComputeIntersectionLocationsInternal<
          structdim,
          spacedim>::IntersectionType;

      const auto comm = cache.get_triangulation().get_mpi_communicator();

      DistributedComputeIntersectionLocationsInternal<structdim, spacedim>
        result;

      auto &send_components = result.send_components;
      auto &recv_components = result.recv_components;
      auto &recv_ptrs       = result.recv_ptrs;

      // search for potential owners
      const auto potential_owners = internal::guess_owners_of_entities(
        comm, global_bboxes, intersection_requests, tolerance);

      const auto &potential_owners_ranks   = std::get<0>(potential_owners);
      const auto &potential_owners_ptrs    = std::get<1>(potential_owners);
      const auto &potential_owners_indices = std::get<2>(potential_owners);

      const auto translate = [&](const unsigned int other_rank) {
        const auto ptr = std::find(potential_owners_ranks.begin(),
                                   potential_owners_ranks.end(),
                                   other_rank);

        Assert(ptr != potential_owners_ranks.end(), ExcInternalError());

        const auto other_rank_index =
          std::distance(potential_owners_ranks.begin(), ptr);

        return other_rank_index;
      };

      Assert(
        (marked_vertices.empty()) ||
          (marked_vertices.size() == cache.get_triangulation().n_vertices()),
        ExcMessage(
          "The marked_vertices vector has to be either empty or its size has "
          "to equal the number of vertices of the triangulation."));

      // In the case that a marked_vertices vector has been given and none
      // of its entries is true, we know that this process does not own
      // any of the incoming points (and it will not send any data) so
      // that we can take a short cut.
      const bool has_relevant_vertices =
        (marked_vertices.empty()) ||
        (std::find(marked_vertices.begin(), marked_vertices.end(), true) !=
         marked_vertices.end());

      // intersection between two cells:
      // One rank requests all intersections of owning cell:
      // owning cell index, cgal vertices of cell
      using RequestType =
        std::vector<std::pair<unsigned int, IntersectionRequest>>;
      // Other ranks send back all found intersections for requesting cell:
      // requesting cell index, cgal vertices of found intersections
      using AnswerType =
        std::vector<std::pair<unsigned int, IntersectionAnswer>>;

      const auto create_request = [&](const unsigned int other_rank) {
        const auto other_rank_index = translate(other_rank);

        RequestType request;
        request.reserve(potential_owners_ptrs[other_rank_index + 1] -
                        potential_owners_ptrs[other_rank_index]);

        for (unsigned int i = potential_owners_ptrs[other_rank_index];
             i < potential_owners_ptrs[other_rank_index + 1];
             ++i)
          request.emplace_back(
            potential_owners_indices[i],
            intersection_requests[potential_owners_indices[i]]);

        return request;
      };


      // TODO: this is potentially useful in many cases and it would be nice to
      // have cache.get_locally_owned_cell_bounding_boxes_rtree(marked_vertices)
      const auto construct_locally_owned_cell_bounding_boxes_rtree =
        [&cache](const std::vector<bool> &marked_verts) {
          const auto cell_marked = [&marked_verts](const auto &cell) {
            for (const unsigned int v : cell->vertex_indices())
              if (marked_verts[cell->vertex_index(v)])
                return true;
            return false;
          };

          const auto &boxes_and_cells =
            cache.get_locally_owned_cell_bounding_boxes_rtree();

          if (marked_verts.empty())
            return boxes_and_cells;

          std::vector<std::pair<
            BoundingBox<spacedim>,
            typename Triangulation<dim, spacedim>::active_cell_iterator>>
            potential_boxes_and_cells;

          for (const auto &box_and_cell : boxes_and_cells)
            if (cell_marked(box_and_cell.second))
              potential_boxes_and_cells.emplace_back(box_and_cell);

          return pack_rtree(potential_boxes_and_cells);
        };


      RTree<
        std::pair<BoundingBox<spacedim>,
                  typename Triangulation<dim, spacedim>::active_cell_iterator>>
        marked_cell_tree;

      const auto answer_request =
        [&]([[maybe_unused]] const unsigned int &other_rank,
            const RequestType                   &request) -> AnswerType {
        AnswerType answer;

        if (has_relevant_vertices)
          {
            if (marked_cell_tree.empty())
              {
                marked_cell_tree =
                  construct_locally_owned_cell_bounding_boxes_rtree(
                    marked_vertices);
              }

            // process requests
            for (unsigned int i = 0; i < request.size(); ++i)
              {
                // create a bounding box with tolerance
                const auto bb = BoundingBox<spacedim>(request[i].second)
                                  .create_extended(tolerance);

                for ([[maybe_unused]] const auto &box_cell :
                     marked_cell_tree |
                       boost::geometry::index::adaptors::queried(
                         boost::geometry::index::intersects(bb)))
                  {
#ifdef DEAL_II_WITH_CGAL
                    const auto &cell                   = box_cell.second;
                    const auto &request_index          = request[i].first;
                    auto        requested_intersection = request[i].second;
                    CGALWrappers::resort_dealii_vertices_to_cgal_order(
                      structdim, requested_intersection);

                    const auto &try_intersection =
                      CGALWrappers::get_vertices_in_cgal_order(
                        cell, cache.get_mapping());

                    const auto &found_intersections = CGALWrappers::
                      compute_intersection_of_cells<dim, structdim, spacedim>(
                        try_intersection, requested_intersection, tolerance);

                    if (found_intersections.size() > 0)
                      {
                        for (const auto &found_intersection :
                             found_intersections)
                          {
                            answer.emplace_back(request_index,
                                                found_intersection);

                            send_components.emplace_back(
                              std::make_pair(cell->level(), cell->index()),
                              other_rank,
                              request_index,
                              found_intersection);
                          }
                      }
#else
                    Assert(false, ExcNeedsCGAL());
#endif
                  }
              }
          }

        return answer;
      };

      const auto process_answer = [&](const unsigned int other_rank,
                                      const AnswerType  &answer) {
        for (unsigned int i = 0; i < answer.size(); ++i)
          recv_components.emplace_back(other_rank,
                                       answer[i].first,
                                       answer[i].second);
      };

      Utilities::MPI::ConsensusAlgorithms::selector<RequestType, AnswerType>(
        potential_owners_ranks,
        create_request,
        answer_request,
        process_answer,
        comm);

      // sort according to 1) intersection index and 2) rank (keeping the order
      // of recv components with same indices and ranks)
      std::stable_sort(recv_components.begin(),
                       recv_components.end(),
                       [&](const auto &a, const auto &b) {
                         // intersection index
                         if (std::get<1>(a) != std::get<1>(b))
                           return std::get<1>(a) < std::get<1>(b);

                         // rank
                         return std::get<0>(a) < std::get<0>(b);
                       });

      // sort according to 1) rank and 2) intersection index (keeping the
      // order of recv components with same indices and ranks)
      std::stable_sort(send_components.begin(),
                       send_components.end(),
                       [&](const auto &a, const auto &b) {
                         // rank
                         if (std::get<1>(a) != std::get<1>(b))
                           return std::get<1>(a) < std::get<1>(b);

                         // intersection idx
                         return std::get<2>(a) < std::get<2>(b);
                       });

      // construct recv_ptrs
      recv_ptrs.assign(intersection_requests.size() + 1, 0);
      for (const auto &rc : recv_components)
        ++recv_ptrs[std::get<1>(rc) + 1];
      for (unsigned int i = 0; i < intersection_requests.size(); ++i)
        recv_ptrs[i + 1] += recv_ptrs[i];

      return result;
    }

  } // namespace internal



  template <int spacedim>
  unsigned int
  find_closest_vertex(const std::map<unsigned int, Point<spacedim>> &vertices,
                      const Point<spacedim>                         &p)
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
    const Point<spacedim>      &p,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
                            &cell_hint,
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
    [[maybe_unused]] const std::vector<BoundingBox<spacedim>> &local_bboxes,
    [[maybe_unused]] const MPI_Comm                            mpi_communicator)
  {
#ifndef DEAL_II_WITH_MPI
    Assert(false,
           ExcMessage(
             "GridTools::exchange_local_bounding_boxes() requires MPI."));
    return {};
#else
    // Step 1: preparing data to be sent
    unsigned int n_bboxes = local_bboxes.size();
    // Dimension of the array to be exchanged (number of double)
    int n_local_data = 2 * spacedim * n_bboxes;
    // data array stores each entry of each point describing the bounding
    // boxes
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

    // Now computing the displacement, relative to recvbuf,
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
    [[maybe_unused]] const MPI_Comm           mpi_communicator)
  {
#ifndef DEAL_II_WITH_MPI
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
    // The i'th element of the following vector contains the index of the
    // first local bounding box from the process of rank i
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
    const Triangulation<dim, spacedim>                &tria,
    std::map<unsigned int, std::vector<unsigned int>> &coinciding_vertex_groups,
    std::map<unsigned int, unsigned int> &vertex_to_coinciding_vertex_group)
  {
    // 1) determine for each vertex a vertex it coincides with and
    //    put it into a map
    {
      // loop over all periodic face pairs
      for (const auto &pair : tria.get_periodic_face_map())
        {
          if (pair.first.first->level() != pair.second.first.first->level())
            continue;

          const auto face_a = pair.first.first->face(pair.first.second);
          const auto face_b =
            pair.second.first.first->face(pair.second.first.second);
          const auto reference_cell       = pair.first.first->reference_cell();
          const auto face_reference_cell  = face_a->reference_cell();
          const auto combined_orientation = pair.second.second;
          const auto inverse_combined_orientation =
            face_reference_cell.get_inverse_combined_orientation(
              combined_orientation);

          AssertDimension(face_a->n_vertices(), face_b->n_vertices());

          // loop over all vertices on face
          for (unsigned int i = 0; i < face_a->n_vertices(); ++i)
            {
              // find the right local vertex index for the second face
              const unsigned int j =
                reference_cell.standard_to_real_face_vertex(
                  i, pair.first.second, inverse_combined_orientation);

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
    for (const auto &cell :
         tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      for (const unsigned int v : cell->vertex_indices())
        vertex_of_own_cell[cell->vertex_index(v)] = true;

    // 3) for each vertex belonging to a locally owned cell, find all ghost
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



  namespace internal
  {
    template <int          dim,
              unsigned int n_vertices,
              unsigned int n_sub_vertices,
              unsigned int n_configurations,
              unsigned int n_lines,
              unsigned int n_cols,
              typename value_type>
    void
    process_sub_cell(
      const std::array<unsigned int, n_configurations>      &cut_line_table,
      const ndarray<unsigned int, n_configurations, n_cols> &new_line_table,
      const ndarray<unsigned int, n_lines, 2>       &line_to_vertex_table,
      const std::vector<value_type>                 &ls_values,
      const std::vector<Point<dim>>                 &points,
      const std::vector<unsigned int>               &mask,
      const double                                   iso_level,
      const double                                   tolerance,
      std::vector<Point<dim>>                       &vertices,
      std::vector<CellData<dim == 1 ? 1 : dim - 1>> &cells,
      const bool                                     write_back_cell_data)
    {
      // inspired by https://graphics.stanford.edu/~mdfisher/MarchingCubes.html

      constexpr unsigned int X = static_cast<unsigned int>(-1);

      // determine configuration
      unsigned int configuration = 0;
      for (unsigned int v = 0; v < n_vertices; ++v)
        if (ls_values[mask[v]] < iso_level)
          configuration |= (1 << v);

      // cell is not cut (nothing to do)
      if (cut_line_table[configuration] == 0)
        return;

      // helper function to determine where an edge (between index i and j) is
      // cut - see also: http://paulbourke.net/geometry/polygonise/
      const auto interpolate = [&](const unsigned int i, const unsigned int j) {
        if (std::abs(iso_level - ls_values[mask[i]]) < tolerance)
          return points[mask[i]];
        if (std::abs(iso_level - ls_values[mask[j]]) < tolerance)
          return points[mask[j]];
        if (std::abs(ls_values[mask[i]] - ls_values[mask[j]]) < tolerance)
          return points[mask[i]];

        const double mu = (iso_level - ls_values[mask[i]]) /
                          (ls_values[mask[j]] - ls_values[mask[i]]);

        return Point<dim>(points[mask[i]] +
                          mu * (points[mask[j]] - points[mask[i]]));
      };

      // determine the position where edges are cut (if they are cut)
      std::array<Point<dim>, n_lines> vertex_list_all;
      for (unsigned int l = 0; l < n_lines; ++l)
        if (cut_line_table[configuration] & (1 << l))
          vertex_list_all[l] =
            interpolate(line_to_vertex_table[l][0], line_to_vertex_table[l][1]);

      // merge duplicate vertices if possible
      unsigned int                      local_vertex_count = 0;
      std::array<Point<dim>, n_lines>   vertex_list_reduced;
      std::array<unsigned int, n_lines> local_remap;
      std::fill(local_remap.begin(), local_remap.end(), X);
      for (int i = 0; new_line_table[configuration][i] != X; ++i)
        if (local_remap[new_line_table[configuration][i]] == X)
          {
            vertex_list_reduced[local_vertex_count] =
              vertex_list_all[new_line_table[configuration][i]];
            local_remap[new_line_table[configuration][i]] = local_vertex_count;
            ++local_vertex_count;
          }

      // write back vertices
      const unsigned int n_vertices_old = vertices.size();
      for (unsigned int i = 0; i < local_vertex_count; ++i)
        vertices.push_back(vertex_list_reduced[i]);

      // write back cells
      if (write_back_cell_data && dim > 1)
        {
          for (unsigned int i = 0; new_line_table[configuration][i] != X;
               i += n_sub_vertices)
            {
              cells.resize(cells.size() + 1);
              cells.back().vertices.resize(n_sub_vertices);

              for (unsigned int v = 0; v < n_sub_vertices; ++v)
                cells.back().vertices[v] =
                  local_remap[new_line_table[configuration][i + v]] +
                  n_vertices_old;
            }
        }
    }
  } // namespace internal



  template <int dim, typename VectorType>
  MarchingCubeAlgorithm<dim, VectorType>::MarchingCubeAlgorithm(
    const Mapping<dim, dim>       &mapping,
    const FiniteElement<dim, dim> &fe,
    const unsigned int             n_subdivisions,
    const double                   tolerance)
    : n_subdivisions(n_subdivisions)
    , tolerance(tolerance)
    , fe_values(mapping,
                fe,
                create_quadrature_rule(n_subdivisions),
                update_values | update_quadrature_points)
  {}



  template <int dim, typename VectorType>
  Quadrature<dim>
  MarchingCubeAlgorithm<dim, VectorType>::create_quadrature_rule(
    const unsigned int n_subdivisions)
  {
    std::vector<Point<dim>> quadrature_points;

    if (dim == 1)
      {
        for (unsigned int i = 0; i <= n_subdivisions; ++i)
          quadrature_points.emplace_back(1.0 / n_subdivisions * i);
      }
    else if (dim == 2)
      {
        for (unsigned int j = 0; j <= n_subdivisions; ++j)
          for (unsigned int i = 0; i <= n_subdivisions; ++i)
            quadrature_points.emplace_back(1.0 / n_subdivisions * i,
                                           1.0 / n_subdivisions * j);
      }
    else
      {
        for (unsigned int k = 0; k <= n_subdivisions; ++k)
          for (unsigned int j = 0; j <= n_subdivisions; ++j)
            for (unsigned int i = 0; i <= n_subdivisions; ++i)
              quadrature_points.emplace_back(1.0 / n_subdivisions * i,
                                             1.0 / n_subdivisions * j,
                                             1.0 / n_subdivisions * k);
      }


    return {quadrature_points};
  }



  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process(
    const DoFHandler<dim>                         &background_dof_handler,
    const VectorType                              &ls_vector,
    const double                                   iso_level,
    std::vector<Point<dim>>                       &vertices,
    std::vector<CellData<dim == 1 ? 1 : dim - 1>> &cells) const
  {
    AssertThrow(
      dim > 1,
      ExcMessage(
        "Not implemented for dim==1. Use the alternative process()-function "
        "not returning a vector of CellData objects."));

    for (const auto &cell : background_dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      process_cell(cell, ls_vector, iso_level, vertices, cells);
  }

  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process(
    const DoFHandler<dim>   &background_dof_handler,
    const VectorType        &ls_vector,
    const double             iso_level,
    std::vector<Point<dim>> &vertices) const
  {
    for (const auto &cell : background_dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      process_cell(cell, ls_vector, iso_level, vertices);

    delete_duplicated_vertices(vertices, 1e-10 /*tol*/);
  }


  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &ls_vector,
    const double                                          iso_level,
    std::vector<Point<dim>>                              &vertices,
    std::vector<CellData<dim == 1 ? 1 : dim - 1>>        &cells) const
  {
    AssertThrow(
      dim > 1,
      ExcMessage(
        "Not implemented for dim==1. Use the alternative process_cell()-function "
        "not returning a vector of CellData objects."));

    std::vector<value_type> ls_values;

    fe_values.reinit(cell);
    ls_values.resize(fe_values.n_quadrature_points);
    fe_values.get_function_values(ls_vector, ls_values);
    process_cell(
      ls_values, fe_values.get_quadrature_points(), iso_level, vertices, cells);
  }

  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &ls_vector,
    const double                                          iso_level,
    std::vector<Point<dim>>                              &vertices) const
  {
    // This vector is just a placeholder to reuse the process_cell function.
    std::vector<CellData<dim == 1 ? 1 : dim - 1>> dummy_cells;

    std::vector<value_type> ls_values;

    fe_values.reinit(cell);
    ls_values.resize(fe_values.n_quadrature_points);
    fe_values.get_function_values(ls_vector, ls_values);

    process_cell(ls_values,
                 fe_values.get_quadrature_points(),
                 iso_level,
                 vertices,
                 dummy_cells,
                 false /*don't write back cell data*/);
  }


  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process_cell(
    std::vector<value_type>                       &ls_values,
    const std::vector<Point<dim>>                 &points,
    const double                                   iso_level,
    std::vector<Point<dim>>                       &vertices,
    std::vector<CellData<dim == 1 ? 1 : dim - 1>> &cells,
    const bool                                     write_back_cell_data) const
  {
    const unsigned p = n_subdivisions + 1;

    if (dim == 1)
      {
        for (unsigned int i = 0; i < n_subdivisions; ++i)
          {
            std::vector<unsigned int> mask{i + 0, i + 1};

            // check if a corner node is cut
            if (std::abs(iso_level - ls_values[mask[0]]) < tolerance)
              vertices.emplace_back(points[mask[0]]);
            else if (std::abs(iso_level - ls_values[mask[1]]) < tolerance)
              {
                if (i + 1 == n_subdivisions)
                  vertices.emplace_back(points[mask[1]]);
              }
            // check if the edge is cut
            else if (((ls_values[mask[0]] > iso_level) &&
                      (ls_values[mask[1]] < iso_level)) ||
                     ((ls_values[mask[0]] < iso_level) &&
                      (ls_values[mask[1]] > iso_level)))
              {
                // determine the interpolation weight (0<mu<1)
                const double mu = (iso_level - ls_values[mask[0]]) /
                                  (ls_values[mask[1]] - ls_values[mask[0]]);

                // interpolate
                vertices.emplace_back(points[mask[0]] +
                                      mu * (points[mask[1]] - points[mask[0]]));
              }
          }
      }
    else if (dim == 2)
      {
        for (unsigned int j = 0; j < n_subdivisions; ++j)
          for (unsigned int i = 0; i < n_subdivisions; ++i)
            {
              std::vector<unsigned int> mask{p * (j + 0) + (i + 0),
                                             p * (j + 0) + (i + 1),
                                             p * (j + 1) + (i + 1),
                                             p * (j + 1) + (i + 0)};

              process_sub_cell(ls_values,
                               points,
                               mask,
                               iso_level,
                               vertices,
                               cells,
                               write_back_cell_data);
            }
      }
    else if (dim == 3)
      {
        for (unsigned int k = 0; k < n_subdivisions; ++k)
          for (unsigned int j = 0; j < n_subdivisions; ++j)
            for (unsigned int i = 0; i < n_subdivisions; ++i)
              {
                std::vector<unsigned int> mask{
                  p * p * (k + 0) + p * (j + 0) + (i + 0),
                  p * p * (k + 0) + p * (j + 0) + (i + 1),
                  p * p * (k + 0) + p * (j + 1) + (i + 1),
                  p * p * (k + 0) + p * (j + 1) + (i + 0),
                  p * p * (k + 1) + p * (j + 0) + (i + 0),
                  p * p * (k + 1) + p * (j + 0) + (i + 1),
                  p * p * (k + 1) + p * (j + 1) + (i + 1),
                  p * p * (k + 1) + p * (j + 1) + (i + 0)};

                process_sub_cell(ls_values,
                                 points,
                                 mask,
                                 iso_level,
                                 vertices,
                                 cells,
                                 write_back_cell_data);
              }
      }
  }



  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process_sub_cell(
    const std::vector<value_type>   &ls_values,
    const std::vector<Point<2>>     &points,
    const std::vector<unsigned int> &mask,
    const double                     iso_level,
    std::vector<Point<2>>           &vertices,
    std::vector<CellData<1>>        &cells,
    const bool                       write_back_cell_data) const
  {
    // set up dimension-dependent sizes and tables
    constexpr unsigned int n_vertices       = 4;
    constexpr unsigned int n_sub_vertices   = 2;
    constexpr unsigned int n_lines          = 4;
    constexpr unsigned int n_configurations = Utilities::pow(2, n_vertices);
    constexpr unsigned int X                = static_cast<unsigned int>(-1);

    // table that indicates if an edge is cut (if the i-th bit is set the i-th
    // line is cut)
    constexpr std::array<unsigned int, n_configurations> cut_line_table = {
      {0b0000,
       0b0101,
       0b0110,
       0b0011,
       0b1010,
       0b0000,
       0b1100,
       0b1001,
       0b1001,
       0b1100,
       0b0000,
       0b1010,
       0b0011,
       0b0110,
       0b0101,
       0b0000}};

    // list of the definition of the newly created lines (each line is defined
    // by two edges it cuts)
    constexpr ndarray<unsigned int, n_configurations, 5> new_line_table = {
      {{{X, X, X, X, X}},
       {{0, 2, X, X, X}},
       {{1, 2, X, X, X}},
       {{0, 1, X, X, X}},
       {{1, 3, X, X, X}},
       {{X, X, X, X, X}},
       {{2, 3, X, X, X}},
       {{0, 3, X, X, X}},
       {{0, 3, X, X, X}},
       {{2, 3, X, X, X}},
       {{X, X, X, X, X}},
       {{1, 3, X, X, X}},
       {{0, 1, X, X, X}},
       {{2, 1, X, X, X}},
       {{0, 2, X, X, X}},
       {{X, X, X, X, X}}}};

    // vertices of each line
    constexpr ndarray<unsigned int, n_lines, 2> line_to_vertex_table = {
      {{{0, 3}}, {{1, 2}}, {{0, 1}}, {{3, 2}}}};

    // run dimension-independent code
    internal::process_sub_cell<2,
                               n_vertices,
                               n_sub_vertices,
                               n_configurations,
                               n_lines,
                               5>(cut_line_table,
                                  new_line_table,
                                  line_to_vertex_table,
                                  ls_values,
                                  points,
                                  mask,
                                  iso_level,
                                  tolerance,
                                  vertices,
                                  cells,
                                  write_back_cell_data);
  }



  template <int dim, typename VectorType>
  void
  MarchingCubeAlgorithm<dim, VectorType>::process_sub_cell(
    const std::vector<value_type>   &ls_values,
    const std::vector<Point<3>>     &points,
    const std::vector<unsigned int> &mask,
    const double                     iso_level,
    std::vector<Point<3>>           &vertices,
    std::vector<CellData<2>>        &cells,
    const bool                       write_back_cell_data) const
  {
    // set up dimension-dependent sizes and tables
    constexpr unsigned int n_vertices       = 8;
    constexpr unsigned int n_sub_vertices   = 3;
    constexpr unsigned int n_lines          = 12;
    constexpr unsigned int n_configurations = Utilities::pow(2, n_vertices);
    constexpr unsigned int X                = static_cast<unsigned int>(-1);

    // clang-format off
    // table that indicates if an edge is cut (if the i-th bit is set the i-th
    // line is cut)
    constexpr std::array<unsigned int, n_configurations> cut_line_table = {{
      0x0,   0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905,
      0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99,  0x393, 0x29a,
      0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93,
      0xf99, 0xe90, 0x230, 0x339, 0x33,  0x13a, 0x636, 0x73f, 0x435, 0x53c,
      0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9,
      0x1a3, 0xaa,  0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6,
      0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569, 0x663, 0x76a, 0x66,  0x16f,
      0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
      0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,  0x3f5, 0x2fc, 0xdfc, 0xcf5,
      0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a,
      0x256, 0x35f, 0x55,  0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
      0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
      0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9,
      0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,  0x1c5, 0x2cf, 0x3c6,
      0x4ca, 0x5c3, 0x6c9, 0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f,
      0xf55, 0xe5c, 0x15c, 0x55,  0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
      0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5,
      0xff,  0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a,
      0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66,  0x76a, 0x663,
      0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
      0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa,  0x1a3, 0x2a9, 0x3a0, 0xd30, 0xc39,
      0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636,
      0x13a, 0x33,  0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f,
      0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99,  0x190,
      0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605,
      0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0}};
    // clang-format on

    // list of the definition of the newly created triangles (each triangles is
    // defined by two edges it cuts)
    constexpr ndarray<unsigned int, n_configurations, 16> new_line_table = {
      {{{X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{0, 8, 3, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{0, 1, 9, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{1, 8, 3, 9, 8, 1, X, X, X, X, X, X, X, X, X, X}},
       {{1, 2, 10, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{0, 8, 3, 1, 2, 10, X, X, X, X, X, X, X, X, X, X}},
       {{9, 2, 10, 0, 2, 9, X, X, X, X, X, X, X, X, X, X}},
       {{2, 8, 3, 2, 10, 8, 10, 9, 8, X, X, X, X, X, X, X}},
       {{3, 11, 2, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{0, 11, 2, 8, 11, 0, X, X, X, X, X, X, X, X, X, X}},
       {{1, 9, 0, 2, 3, 11, X, X, X, X, X, X, X, X, X, X}},
       {{1, 11, 2, 1, 9, 11, 9, 8, 11, X, X, X, X, X, X, X}},
       {{3, 10, 1, 11, 10, 3, X, X, X, X, X, X, X, X, X, X}},
       {{0, 10, 1, 0, 8, 10, 8, 11, 10, X, X, X, X, X, X, X}},
       {{3, 9, 0, 3, 11, 9, 11, 10, 9, X, X, X, X, X, X, X}},
       {{9, 8, 10, 10, 8, 11, X, X, X, X, X, X, X, X, X, X}},
       {{4, 7, 8, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{4, 3, 0, 7, 3, 4, X, X, X, X, X, X, X, X, X, X}},
       {{0, 1, 9, 8, 4, 7, X, X, X, X, X, X, X, X, X, X}},
       {{4, 1, 9, 4, 7, 1, 7, 3, 1, X, X, X, X, X, X, X}},
       {{1, 2, 10, 8, 4, 7, X, X, X, X, X, X, X, X, X, X}},
       {{3, 4, 7, 3, 0, 4, 1, 2, 10, X, X, X, X, X, X, X}},
       {{9, 2, 10, 9, 0, 2, 8, 4, 7, X, X, X, X, X, X, X}},
       {{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, X, X, X, X}},
       {{8, 4, 7, 3, 11, 2, X, X, X, X, X, X, X, X, X, X}},
       {{11, 4, 7, 11, 2, 4, 2, 0, 4, X, X, X, X, X, X, X}},
       {{9, 0, 1, 8, 4, 7, 2, 3, 11, X, X, X, X, X, X, X}},
       {{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, X, X, X, X}},
       {{3, 10, 1, 3, 11, 10, 7, 8, 4, X, X, X, X, X, X, X}},
       {{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, X, X, X, X}},
       {{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, X, X, X, X}},
       {{4, 7, 11, 4, 11, 9, 9, 11, 10, X, X, X, X, X, X, X}},
       {{9, 5, 4, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{9, 5, 4, 0, 8, 3, X, X, X, X, X, X, X, X, X, X}},
       {{0, 5, 4, 1, 5, 0, X, X, X, X, X, X, X, X, X, X}},
       {{8, 5, 4, 8, 3, 5, 3, 1, 5, X, X, X, X, X, X, X}},
       {{1, 2, 10, 9, 5, 4, X, X, X, X, X, X, X, X, X, X}},
       {{3, 0, 8, 1, 2, 10, 4, 9, 5, X, X, X, X, X, X, X}},
       {{5, 2, 10, 5, 4, 2, 4, 0, 2, X, X, X, X, X, X, X}},
       {{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, X, X, X, X}},
       {{9, 5, 4, 2, 3, 11, X, X, X, X, X, X, X, X, X, X}},
       {{0, 11, 2, 0, 8, 11, 4, 9, 5, X, X, X, X, X, X, X}},
       {{0, 5, 4, 0, 1, 5, 2, 3, 11, X, X, X, X, X, X, X}},
       {{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, X, X, X, X}},
       {{10, 3, 11, 10, 1, 3, 9, 5, 4, X, X, X, X, X, X, X}},
       {{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, X, X, X, X}},
       {{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, X, X, X, X}},
       {{5, 4, 8, 5, 8, 10, 10, 8, 11, X, X, X, X, X, X, X}},
       {{9, 7, 8, 5, 7, 9, X, X, X, X, X, X, X, X, X, X}},
       {{9, 3, 0, 9, 5, 3, 5, 7, 3, X, X, X, X, X, X, X}},
       {{0, 7, 8, 0, 1, 7, 1, 5, 7, X, X, X, X, X, X, X}},
       {{1, 5, 3, 3, 5, 7, X, X, X, X, X, X, X, X, X, X}},
       {{9, 7, 8, 9, 5, 7, 10, 1, 2, X, X, X, X, X, X, X}},
       {{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, X, X, X, X}},
       {{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, X, X, X, X}},
       {{2, 10, 5, 2, 5, 3, 3, 5, 7, X, X, X, X, X, X, X}},
       {{7, 9, 5, 7, 8, 9, 3, 11, 2, X, X, X, X, X, X, X}},
       {{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, X, X, X, X}},
       {{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, X, X, X, X}},
       {{11, 2, 1, 11, 1, 7, 7, 1, 5, X, X, X, X, X, X, X}},
       {{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, X, X, X, X}},
       {{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, X}},
       {{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, X}},
       {{11, 10, 5, 7, 11, 5, X, X, X, X, X, X, X, X, X, X}},
       {{10, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{0, 8, 3, 5, 10, 6, X, X, X, X, X, X, X, X, X, X}},
       {{9, 0, 1, 5, 10, 6, X, X, X, X, X, X, X, X, X, X}},
       {{1, 8, 3, 1, 9, 8, 5, 10, 6, X, X, X, X, X, X, X}},
       {{1, 6, 5, 2, 6, 1, X, X, X, X, X, X, X, X, X, X}},
       {{1, 6, 5, 1, 2, 6, 3, 0, 8, X, X, X, X, X, X, X}},
       {{9, 6, 5, 9, 0, 6, 0, 2, 6, X, X, X, X, X, X, X}},
       {{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, X, X, X, X}},
       {{2, 3, 11, 10, 6, 5, X, X, X, X, X, X, X, X, X, X}},
       {{11, 0, 8, 11, 2, 0, 10, 6, 5, X, X, X, X, X, X, X}},
       {{0, 1, 9, 2, 3, 11, 5, 10, 6, X, X, X, X, X, X, X}},
       {{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, X, X, X, X}},
       {{6, 3, 11, 6, 5, 3, 5, 1, 3, X, X, X, X, X, X, X}},
       {{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, X, X, X, X}},
       {{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, X, X, X, X}},
       {{6, 5, 9, 6, 9, 11, 11, 9, 8, X, X, X, X, X, X, X}},
       {{5, 10, 6, 4, 7, 8, X, X, X, X, X, X, X, X, X, X}},
       {{4, 3, 0, 4, 7, 3, 6, 5, 10, X, X, X, X, X, X, X}},
       {{1, 9, 0, 5, 10, 6, 8, 4, 7, X, X, X, X, X, X, X}},
       {{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, X, X, X, X}},
       {{6, 1, 2, 6, 5, 1, 4, 7, 8, X, X, X, X, X, X, X}},
       {{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, X, X, X, X}},
       {{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, X, X, X, X}},
       {{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, X}},
       {{3, 11, 2, 7, 8, 4, 10, 6, 5, X, X, X, X, X, X, X}},
       {{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, X, X, X, X}},
       {{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, X, X, X, X}},
       {{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, X}},
       {{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, X, X, X, X}},
       {{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, X}},
       {{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, X}},
       {{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, X, X, X, X}},
       {{10, 4, 9, 6, 4, 10, X, X, X, X, X, X, X, X, X, X}},
       {{4, 10, 6, 4, 9, 10, 0, 8, 3, X, X, X, X, X, X, X}},
       {{10, 0, 1, 10, 6, 0, 6, 4, 0, X, X, X, X, X, X, X}},
       {{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, X, X, X, X}},
       {{1, 4, 9, 1, 2, 4, 2, 6, 4, X, X, X, X, X, X, X}},
       {{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, X, X, X, X}},
       {{0, 2, 4, 4, 2, 6, X, X, X, X, X, X, X, X, X, X}},
       {{8, 3, 2, 8, 2, 4, 4, 2, 6, X, X, X, X, X, X, X}},
       {{10, 4, 9, 10, 6, 4, 11, 2, 3, X, X, X, X, X, X, X}},
       {{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, X, X, X, X}},
       {{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, X, X, X, X}},
       {{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, X}},
       {{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, X, X, X, X}},
       {{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, X}},
       {{3, 11, 6, 3, 6, 0, 0, 6, 4, X, X, X, X, X, X, X}},
       {{6, 4, 8, 11, 6, 8, X, X, X, X, X, X, X, X, X, X}},
       {{7, 10, 6, 7, 8, 10, 8, 9, 10, X, X, X, X, X, X, X}},
       {{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, X, X, X, X}},
       {{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, X, X, X, X}},
       {{10, 6, 7, 10, 7, 1, 1, 7, 3, X, X, X, X, X, X, X}},
       {{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, X, X, X, X}},
       {{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, X}},
       {{7, 8, 0, 7, 0, 6, 6, 0, 2, X, X, X, X, X, X, X}},
       {{7, 3, 2, 6, 7, 2, X, X, X, X, X, X, X, X, X, X}},
       {{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, X, X, X, X}},
       {{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, X}},
       {{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, X}},
       {{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, X, X, X, X}},
       {{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, X}},
       {{0, 9, 1, 11, 6, 7, X, X, X, X, X, X, X, X, X, X}},
       {{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, X, X, X, X}},
       {{7, 11, 6, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{7, 6, 11, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{3, 0, 8, 11, 7, 6, X, X, X, X, X, X, X, X, X, X}},
       {{0, 1, 9, 11, 7, 6, X, X, X, X, X, X, X, X, X, X}},
       {{8, 1, 9, 8, 3, 1, 11, 7, 6, X, X, X, X, X, X, X}},
       {{10, 1, 2, 6, 11, 7, X, X, X, X, X, X, X, X, X, X}},
       {{1, 2, 10, 3, 0, 8, 6, 11, 7, X, X, X, X, X, X, X}},
       {{2, 9, 0, 2, 10, 9, 6, 11, 7, X, X, X, X, X, X, X}},
       {{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, X, X, X, X}},
       {{7, 2, 3, 6, 2, 7, X, X, X, X, X, X, X, X, X, X}},
       {{7, 0, 8, 7, 6, 0, 6, 2, 0, X, X, X, X, X, X, X}},
       {{2, 7, 6, 2, 3, 7, 0, 1, 9, X, X, X, X, X, X, X}},
       {{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, X, X, X, X}},
       {{10, 7, 6, 10, 1, 7, 1, 3, 7, X, X, X, X, X, X, X}},
       {{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, X, X, X, X}},
       {{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, X, X, X, X}},
       {{7, 6, 10, 7, 10, 8, 8, 10, 9, X, X, X, X, X, X, X}},
       {{6, 8, 4, 11, 8, 6, X, X, X, X, X, X, X, X, X, X}},
       {{3, 6, 11, 3, 0, 6, 0, 4, 6, X, X, X, X, X, X, X}},
       {{8, 6, 11, 8, 4, 6, 9, 0, 1, X, X, X, X, X, X, X}},
       {{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, X, X, X, X}},
       {{6, 8, 4, 6, 11, 8, 2, 10, 1, X, X, X, X, X, X, X}},
       {{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, X, X, X, X}},
       {{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, X, X, X, X}},
       {{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, X}},
       {{8, 2, 3, 8, 4, 2, 4, 6, 2, X, X, X, X, X, X, X}},
       {{0, 4, 2, 4, 6, 2, X, X, X, X, X, X, X, X, X, X}},
       {{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, X, X, X, X}},
       {{1, 9, 4, 1, 4, 2, 2, 4, 6, X, X, X, X, X, X, X}},
       {{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, X, X, X, X}},
       {{10, 1, 0, 10, 0, 6, 6, 0, 4, X, X, X, X, X, X, X}},
       {{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, X}},
       {{10, 9, 4, 6, 10, 4, X, X, X, X, X, X, X, X, X, X}},
       {{4, 9, 5, 7, 6, 11, X, X, X, X, X, X, X, X, X, X}},
       {{0, 8, 3, 4, 9, 5, 11, 7, 6, X, X, X, X, X, X, X}},
       {{5, 0, 1, 5, 4, 0, 7, 6, 11, X, X, X, X, X, X, X}},
       {{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, X, X, X, X}},
       {{9, 5, 4, 10, 1, 2, 7, 6, 11, X, X, X, X, X, X, X}},
       {{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, X, X, X, X}},
       {{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, X, X, X, X}},
       {{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, X}},
       {{7, 2, 3, 7, 6, 2, 5, 4, 9, X, X, X, X, X, X, X}},
       {{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, X, X, X, X}},
       {{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, X, X, X, X}},
       {{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, X}},
       {{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, X, X, X, X}},
       {{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, X}},
       {{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, X}},
       {{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, X, X, X, X}},
       {{6, 9, 5, 6, 11, 9, 11, 8, 9, X, X, X, X, X, X, X}},
       {{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, X, X, X, X}},
       {{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, X, X, X, X}},
       {{6, 11, 3, 6, 3, 5, 5, 3, 1, X, X, X, X, X, X, X}},
       {{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, X, X, X, X}},
       {{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, X}},
       {{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, X}},
       {{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, X, X, X, X}},
       {{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, X, X, X, X}},
       {{9, 5, 6, 9, 6, 0, 0, 6, 2, X, X, X, X, X, X, X}},
       {{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, X}},
       {{1, 5, 6, 2, 1, 6, X, X, X, X, X, X, X, X, X, X}},
       {{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, X}},
       {{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, X, X, X, X}},
       {{0, 3, 8, 5, 6, 10, X, X, X, X, X, X, X, X, X, X}},
       {{10, 5, 6, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{11, 5, 10, 7, 5, 11, X, X, X, X, X, X, X, X, X, X}},
       {{11, 5, 10, 11, 7, 5, 8, 3, 0, X, X, X, X, X, X, X}},
       {{5, 11, 7, 5, 10, 11, 1, 9, 0, X, X, X, X, X, X, X}},
       {{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, X, X, X, X}},
       {{11, 1, 2, 11, 7, 1, 7, 5, 1, X, X, X, X, X, X, X}},
       {{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, X, X, X, X}},
       {{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, X, X, X, X}},
       {{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, X}},
       {{2, 5, 10, 2, 3, 5, 3, 7, 5, X, X, X, X, X, X, X}},
       {{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, X, X, X, X}},
       {{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, X, X, X, X}},
       {{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, X}},
       {{1, 3, 5, 3, 7, 5, X, X, X, X, X, X, X, X, X, X}},
       {{0, 8, 7, 0, 7, 1, 1, 7, 5, X, X, X, X, X, X, X}},
       {{9, 0, 3, 9, 3, 5, 5, 3, 7, X, X, X, X, X, X, X}},
       {{9, 8, 7, 5, 9, 7, X, X, X, X, X, X, X, X, X, X}},
       {{5, 8, 4, 5, 10, 8, 10, 11, 8, X, X, X, X, X, X, X}},
       {{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, X, X, X, X}},
       {{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, X, X, X, X}},
       {{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, X}},
       {{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, X, X, X, X}},
       {{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, X}},
       {{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, X}},
       {{9, 4, 5, 2, 11, 3, X, X, X, X, X, X, X, X, X, X}},
       {{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, X, X, X, X}},
       {{5, 10, 2, 5, 2, 4, 4, 2, 0, X, X, X, X, X, X, X}},
       {{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, X}},
       {{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, X, X, X, X}},
       {{8, 4, 5, 8, 5, 3, 3, 5, 1, X, X, X, X, X, X, X}},
       {{0, 4, 5, 1, 0, 5, X, X, X, X, X, X, X, X, X, X}},
       {{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, X, X, X, X}},
       {{9, 4, 5, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{4, 11, 7, 4, 9, 11, 9, 10, 11, X, X, X, X, X, X, X}},
       {{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, X, X, X, X}},
       {{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, X, X, X, X}},
       {{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, X}},
       {{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, X, X, X, X}},
       {{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, X}},
       {{11, 7, 4, 11, 4, 2, 2, 4, 0, X, X, X, X, X, X, X}},
       {{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, X, X, X, X}},
       {{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, X, X, X, X}},
       {{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, X}},
       {{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, X}},
       {{1, 10, 2, 8, 7, 4, X, X, X, X, X, X, X, X, X, X}},
       {{4, 9, 1, 4, 1, 7, 7, 1, 3, X, X, X, X, X, X, X}},
       {{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, X, X, X, X}},
       {{4, 0, 3, 7, 4, 3, X, X, X, X, X, X, X, X, X, X}},
       {{4, 8, 7, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{9, 10, 8, 10, 11, 8, X, X, X, X, X, X, X, X, X, X}},
       {{3, 0, 9, 3, 9, 11, 11, 9, 10, X, X, X, X, X, X, X}},
       {{0, 1, 10, 0, 10, 8, 8, 10, 11, X, X, X, X, X, X, X}},
       {{3, 1, 10, 11, 3, 10, X, X, X, X, X, X, X, X, X, X}},
       {{1, 2, 11, 1, 11, 9, 9, 11, 8, X, X, X, X, X, X, X}},
       {{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, X, X, X, X}},
       {{0, 2, 11, 8, 0, 11, X, X, X, X, X, X, X, X, X, X}},
       {{3, 2, 11, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{2, 3, 8, 2, 8, 10, 10, 8, 9, X, X, X, X, X, X, X}},
       {{9, 10, 2, 0, 9, 2, X, X, X, X, X, X, X, X, X, X}},
       {{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, X, X, X, X}},
       {{1, 10, 2, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{1, 3, 8, 9, 1, 8, X, X, X, X, X, X, X, X, X, X}},
       {{0, 9, 1, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{0, 3, 8, X, X, X, X, X, X, X, X, X, X, X, X, X}},
       {{X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X}}}};

    // vertices of each line
    static constexpr ndarray<unsigned int, n_lines, 2> line_to_vertex_table = {
      {{{0, 1}},
       {{1, 2}},
       {{2, 3}},
       {{3, 0}},
       {{4, 5}},
       {{5, 6}},
       {{6, 7}},
       {{7, 4}},
       {{0, 4}},
       {{1, 5}},
       {{2, 6}},
       {{3, 7}}}};

    // run dimension-independent code
    internal::process_sub_cell<3,
                               n_vertices,
                               n_sub_vertices,
                               n_configurations,
                               n_lines,
                               16>(cut_line_table,
                                   new_line_table,
                                   line_to_vertex_table,
                                   ls_values,
                                   points,
                                   mask,
                                   iso_level,
                                   tolerance,
                                   vertices,
                                   cells,
                                   write_back_cell_data);
  }

} /* namespace GridTools */


// explicit instantiations
#include "grid/grid_tools.inst"

DEAL_II_NAMESPACE_CLOSE

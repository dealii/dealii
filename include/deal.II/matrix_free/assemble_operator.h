// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
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

#pragma once

#include <deal.II/lac/vector.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/graph/smallest_last_ordering.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/shared_array_property_map.hpp>

#include <set>
#include <vector>

// TODO distributed graph

namespace dealii
{
  namespace MatrixFreeUtilities
  {
    // Typedefs for boost graph types.
    namespace
    {
      using Graph =
        boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
      using vertex_descriptor  = boost::graph_traits<Graph>::vertex_descriptor;
      using vertices_size_type = boost::graph_traits<Graph>::vertices_size_type;
      using vertex_index_map =
        boost::property_map<Graph, boost::vertex_index_t>::const_type;
      using Edge = std::pair<vertices_size_type, vertices_size_type>;
    } // namespace

    /**
     * Cache the graph generated from a sparsity pattern for later computations.
     * This is useful if the matrix does not change its pattern.
     * // TODO cache entries
     */
    struct GraphCache
    {
      Graph                           graph;
      std::vector<vertices_size_type> color_vec;
      boost::iterator_property_map<vertices_size_type *, vertex_index_map>
                                                   coloring;
      vertices_size_type                           num_colors;
      std::vector<std::vector<vertices_size_type>> entries_per_color;
      unsigned int                                 max_entries_per_row = 0;

      /**
       * Create a boost graph (V,E) from the provided sparsity pattern.
       * The vertices V are given by the column indices.
       * An edge (v1,v2) is an element of E, if there is a row r,
       * such that (r,v1) and (r,v2) are both entries of the sparsity pattern.
       *
       * @param <tt>pattern</tt> the sparsity pattern for which the graph should be built.
       * Pattern needs to provide: n_rows(), n_cols(), begin(r), end(r),
       * iterator->row(), iterator->colum().
       */
      template <typename Pattern>
      void
      make_graph(const Pattern &pattern)
      {
        std::set<Edge> edges;
        for (unsigned int r = 0; r < pattern.n_rows(); r++)
          {
            for (auto it = pattern.begin(r); it != pattern.end(r); it++)
              {
                for (auto it2 = it; it2 != pattern.end(r); it2++)
                  {
                    AssertDimension(r, it->row());
                    AssertDimension(r, it2->row());
                    auto c1 = it->column();
                    auto c2 = it2->column();
                    if (c1 > c2)
                      std::swap(c1, c2);
                    edges.insert(std::make_pair(c1, c2));
                  }
              }
          }

        // Construct a graph from the edge-set.
        for (const Edge &e : edges)
          add_edge(e.first, e.second, graph);

        AssertDimension(pattern.n_cols(), boost::num_vertices(graph));

        max_entries_per_row = pattern.max_entries_per_row();
      }

      /**
       * Color the existing graph using sequential vertex coloring.
       * Different orderings of the vertices can give different results.
       * By default, we do not reorder them, as the benefit does not seem to be
       * huge.
       *
       * @param <tt>order</tt> defines a reording of the vertices.
       * The given type needs to satisfy the boost::property_map concept.
       * Example:
       *   boost::vector_property_map<vertices_size_type> order;
       *   boost::smallest_last_vertex_ordering(graph, order);
       */
      template <
        typename Order = boost::typed_identity_property_map<vertices_size_type>>
      void
      make_coloring(const Order &&order =
                      boost::typed_identity_property_map<vertices_size_type>{})
      {
        color_vec.resize(boost::num_vertices(graph));
        coloring =
          boost::make_iterator_property_map(&color_vec.front(),
                                            get(boost::vertex_index, graph));
        num_colors = boost::sequential_vertex_coloring(graph, order, coloring);

        entries_per_color.resize(num_colors);
        for (vertices_size_type color = 0; color < num_colors; color++)
          {
            // Columns of color c.
            entries_per_color[color].reserve(max_entries_per_row);
            for (vertices_size_type i = 0; i < color_vec.size(); i++)
              if (color_vec[i] == color)
                entries_per_color[color].push_back(i);
          }
      }
    };

    /**
     * Assemble a given operator into a matrix, using the provided sparsity
     pattern.
     *
     * @param <tt>output</tt> resulting matrix. Needs to provide a function set(row, column, value) that allows to write to each entry of the
     sparsity pattern.

     * @param <tt>op</tt> the operator to assemble. Needs to provide a function vmult(RangeVector&, const DomainVector&).

     * @param <tt>pattern</tt> the underlying sparsity pattern of <tt>op</tt>.
     * Note that <tt>pattern</tt> can be a subset of the sparsity pattern from
     <tt>output</tt>

     * @param <tt>dummy_domain</tt> a vector used to initialize vectors in the domain of <tt>op</tt> by copying the structure of
     <tt>dummy_domain</tt>.
     * The type <tt>DomainVector</tt> needs to provide reinit(const
     DomainVector&, bool=false) to create a zero vector.

     * @param <tt>dummy_range</tt> see <tt>dummy_domain</tt>.

     * @param <tt>cache</tt> an already created CacheGraph object storing intermediate results.
     * This allows to reuse the graph in case the sparsity pattern does not
     change.
     * If cache is empty (nullptr), it will be created.
     */
    template <typename Matrix,
              typename Operator,
              typename Pattern,
              typename DomainVector,
              typename RangeVector>
    void
    assemble_operator(
      Matrix &                     output,
      const Operator &             op,
      const Pattern &              pattern,
      const DomainVector &         dummy_domain,
      const RangeVector &          dummy_range,
      std::shared_ptr<GraphCache> &cache = std::shared_ptr<GraphCache>())
    {
      // Early exit for empty pattern.
      if (pattern.n_nonzero_elements() == 0)
        return;

      // If not cache is given, initialize a new object.
      if (!cache)
        {
          cache = std::make_shared<GraphCache>();

          // Create a set with all the edges defined by the sparsity pattern.
          // An edge E = (c1,c2) is contained iff there exists a row r,
          // such that (r,c1) and (r,c2) are both contained in the pattern.
          cache->make_graph(pattern);
          cache->make_coloring();
        }

      // Coloring: for each edge (v1,v2) the colors of the vertices is
      // different, i.e. c(v1) != c(v2). To obtain the matrix from the operator,
      // we need to test (multiply) with one vector e per color. For each color
      // c, e[i] = 1 if c(v[i]) == c else 0 (entries = {i : v[i] == c}). This
      // gives us a vector col = sum_i A e_i = sum_i column_i(A). We need to
      // decompose col into its contributions from the different columns. This
      // is possible, as all the columns listed in entries have no non-zero row
      // in common.
      for (vertices_size_type color = 0; color < cache->num_colors; color++)
        {
          // Columns of color c.
          std::vector<vertices_size_type> entries;
          entries.reserve(pattern.max_entries_per_row());
          for (vertices_size_type i = 0; i < pattern.n_cols(); i++)
            if (cache->color_vec[i] == color)
              entries.push_back(i);

          // Construct the test-vector e = sum e_i for i in entries.
          DomainVector e;
          e.reinit(dummy_domain, false);
          for (auto i : entries)
            e[i] = 1.0;

          // Apply to operator. image = sum A*e_i for i in entries.
          RangeVector image;
          image.reinit(dummy_range, false);
          op.vmult(image, e);

          // Split the result into its contributions from the single unit
          // vectors, i.e. A*e_i. This is possible, since for each row, only one
          // column in entry can store a value.
          for (vertices_size_type r = 0; r < pattern.n_rows(); r++)
            {
              if (image[r] != 0.0)
                {
#ifdef DEBUG
                  bool found = false;
#endif
                  for (auto c : entries)
                    {
                      if (pattern.exists(r, c))
                        {
                          output.set(r, c, image[r]);
#ifdef DEBUG
                          // Verify that we can uniquely identify the column.
                          for (auto j : entries)
                            Assert(c == j || !pattern.exists(r, j),
                                   ExcInternalError());

                          found = true;
#endif
                          break;
                        }
                    }

                  // Make sure that each non-zero entry found its way into the
                  // matrix.
                  Assert(
                    found,
                    ExcMessage(
                      "The provided pattern does not include all required entries."));
                }
            }
        }
    }

    /**
     * Variant of assemble_operator using Vector as domain and range vectors.
     */
    template <typename Matrix, typename Operator, typename Pattern>
    void
    assemble_operator(
      Matrix &                     output,
      const Operator &             op,
      const Pattern &              pattern,
      std::shared_ptr<GraphCache> &cache = std::shared_ptr<GraphCache>())
    {
      Vector<typename Matrix::value_type> domain(pattern.n_cols());
      Vector<typename Matrix::value_type> range(pattern.n_rows());
      assemble_operator(output, op, pattern, domain, range, cache);
    }
  } // namespace MatrixFreeUtilities
} // namespace dealii

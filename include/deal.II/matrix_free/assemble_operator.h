// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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

#ifndef dealii_assemble_operator_h
#define dealii_assemble_operator_h

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

DEAL_II_NAMESPACE_OPEN

// TODO distributed graph

namespace MatrixFreeTools
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
   */
  struct GraphCache
  {
    // The graph object.
    Graph graph;

    // Vector mapping each column index to its color.
    std::vector<vertices_size_type> color_vec;

    // Number of required colors.
    // Equals max(color_vec)+1.
    vertices_size_type num_colors;

    // For each color = 0,.., num_colors-1, entries_per_color[color] stores a
    // list of column indices.
    std::vector<std::vector<vertices_size_type>> entries_per_color;

    /**
     * Create a boost graph (V,E) from the provided sparsity pattern. The
     * vertices V are given by the column indices. An edge (v1,v2) is an element
     * of E, if there is a row r, such that (r,v1) and (r,v2) are both entries
     * of the sparsity pattern.
     *
     * @param pattern The sparsity pattern for which the graph should be built. Pattern needs to provide: n_rows(), n_cols(), begin(r), end(r),
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
                  edges.insert({std::min(it->column(), it2->column()),
                                std::max(it->column(), it2->column())});
                }
            }
        }

      // Construct a graph from the edge-set.
      for (const Edge &e : edges)
        add_edge(e.first, e.second, graph);

      AssertDimension(pattern.n_cols(), boost::num_vertices(graph));
    }

    /**
     * Color the existing graph using sequential vertex coloring. Different
     * orderings of the vertices can give different results. By default, we do
     * not reorder them, as the benefit does not seem to be huge.
     *
     * @param order Defines a reording of the vertices.
     * The given type needs to satisfy the boost::property_map concept.
     * Example:
     *   boost::vector_property_map<vertices_size_type> order;
     *   boost::smallest_last_vertex_ordering(graph, order);
     */
    template <
      typename Order = boost::typed_identity_property_map<vertices_size_type>>
    void
    make_coloring(const Order &order =
                    boost::typed_identity_property_map<vertices_size_type>{})
    {
      using Coloring =
        boost::iterator_property_map<vertices_size_type *, vertex_index_map>;

      // Color the graph.
      color_vec.resize(boost::num_vertices(graph));
      Coloring coloring =
        boost::make_iterator_property_map(&color_vec.front(),
                                          get(boost::vertex_index, graph));
      num_colors = boost::sequential_vertex_coloring(graph, order, coloring);
      AssertDimension(*std::max_element(color_vec.begin(), color_vec.end()) + 1,
                      num_colors);

      // Create list of column indices for each color.
      entries_per_color.resize(num_colors);
      for (vertices_size_type color = 0; color < num_colors; color++)
        {
          // Columns of color c.
          for (vertices_size_type i = 0; i < color_vec.size(); i++)
            if (color_vec[i] == color)
              entries_per_color[color].push_back(i);
        }
    }
  };

  /**
   * Assemble a given linear operator into a matrix, using the provided sparsity
   pattern.
   *
   * @param output_matrix Resulting matrix. Needs to provide a function set(row, column, value) that allows to write to each entry of the sparsity
   pattern.

   * @param linear_operator The operator to assemble. Needs to provide a function vmult(RangeVector&, const DomainVector&).

   * @param pattern The underlying sparsity pattern of @p linear_operator. Note that @p pattern can be a subset of the sparsity pattern from @p output_matrix.

   * @param domain_vector, @param range_vector Vectors used for the operator application. Need to be properly setup to be used as arguments to @p linear_operator.vmult(range_vector, domain_vector).
   *
   * @param cache An already created CacheGraph object storing intermediate results. This allows to reuse the graph in case the sparsity pattern
   does not change. If cache is empty (nullptr), it will be created.
   */
  template <typename Matrix,
            typename Operator,
            typename Pattern,
            typename DomainVector,
            typename RangeVector>
  void
  assemble_operator(Matrix &                     output_matrix,
                    const Operator &             linear_operator,
                    const Pattern &              pattern,
                    DomainVector &               domain_vector,
                    RangeVector &                range_vector,
                    std::shared_ptr<GraphCache> &cache)
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

    // Coloring: for each edge (v1,v2) the colors of the vertices is different,
    // i.e. c(v1) != c(v2). To obtain the matrix from the operator, we need to
    // test (multiply) with one vector e per color. For each color c, e[i] = 1
    // if c(v[i]) == c else 0 (entries = {i : v[i] == c}). This gives us a
    // vector col = sum_i A e_i = sum_i column_i(A). We need to decompose col
    // into its contributions from the different columns. This is possible, as
    // all the columns listed in entries have no non-zero row in common.
    for (vertices_size_type color = 0; color < cache->num_colors; color++)
      {
        // Columns of color c.
        std::vector<vertices_size_type> &entries =
          cache->entries_per_color[color];

        // Construct the test-vector e = sum e_i for i in entries.
        domain_vector = 0;
        for (auto i : entries)
          domain_vector[i] = 1.0;

        // Apply to operator. range_vector = sum A*e_i for i in entries.
        linear_operator.vmult(range_vector, domain_vector);

        // Split the result into its contributions from the single unit vectors,
        // i.e. A*e_i. This is possible, since for each row, only one column in
        // entry can store a value.
        for (vertices_size_type r = 0; r < pattern.n_rows(); r++)
          {
            if (range_vector[r] != 0.0)
              {
#ifdef DEBUG
                bool found = false;
#endif
                for (auto c : entries)
                  {
                    if (pattern.exists(r, c))
                      {
                        output_matrix.set(r, c, range_vector[r]);
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
  assemble_operator(Matrix &                     output_matrix,
                    const Operator &             linear_operator,
                    const Pattern &              pattern,
                    std::shared_ptr<GraphCache> &cache)
  {
    Vector<typename Matrix::value_type> domain(pattern.n_cols());
    Vector<typename Matrix::value_type> range(pattern.n_rows());
    assemble_operator(
      output_matrix, linear_operator, pattern, domain, range, cache);
  }

  /**
   * Variant of assemble_operator without cache.
   */
  template <typename Matrix, typename Operator, typename Pattern>
  std::shared_ptr<GraphCache>
  assemble_operator(Matrix &        output_matrix,
                    const Operator &linear_operator,
                    const Pattern & pattern)
  {
    Vector<typename Matrix::value_type> domain_vector(pattern.n_cols());
    Vector<typename Matrix::value_type> range_vector(pattern.n_rows());
    std::shared_ptr<GraphCache>         cache;
    assemble_operator(output_matrix,
                      linear_operator,
                      pattern,
                      domain_vector,
                      range_vector,
                      cache);
    return cache;
  }

  /**
   * Variant of assemble_operator without cache.
   */
  template <typename Matrix,
            typename Operator,
            typename Pattern,
            typename DomainVector,
            typename RangeVector>
  std::shared_ptr<GraphCache>
  assemble_operator(Matrix &        output_matrix,
                    const Operator &linear_operator,
                    const Pattern & pattern,
                    DomainVector &  domain_vector,
                    RangeVector &   range_vector)
  {
    std::shared_ptr<GraphCache> cache;
    assemble_operator(output_matrix,
                      linear_operator,
                      pattern,
                      domain_vector,
                      range_vector,
                      cache);
    return cache;
  }
} // namespace MatrixFreeTools

DEAL_II_NAMESPACE_CLOSE

#endif
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

struct GraphCache
{
  Graph                           graph;
  std::vector<vertices_size_type> color_vec;
  boost::iterator_property_map<vertices_size_type *, vertex_index_map> coloring;
  vertices_size_type                           num_colors;
  std::vector<std::vector<vertices_size_type>> entries_per_color;
  unsigned int                                 max_entries_per_row = 0;

  void
  make_coloring()
  {
    boost::typed_identity_property_map<vertices_size_type> order;
    // boost::vector_property_map<vertices_size_type> order;
    // boost::smallest_last_vertex_ordering(graph, order);

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
                if (c1 == c2)
                  continue;
                if (c1 > c2)
                  std::swap(c1, c2);
                edges.insert(std::make_pair(c1, c2));
              }
          }
      }

    // Construct a graph from the edge-set.
    for (const Edge &e : edges)
      add_edge(e.first, e.second, graph);

    max_entries_per_row = pattern.max_entries_per_row();
  }
};

template <typename Matrix,
          typename Operator,
          typename Pattern,
          typename DomainVector,
          typename RangeVector>
std::shared_ptr<GraphCache>
assemble_operator(Matrix &                          output,
                  const Operator &                  op,
                  const Pattern &                   pattern,
                  const DomainVector &              dummy_domain,
                  const RangeVector &               dummy_range,
                  const std::shared_ptr<GraphCache> cache)
{
  // Coloring: for each edge (v1, v2) the colors of the vertices is different,
  // i.e. c(v1) != c(v2). To obtain the matrix from the operator, we need to
  // test (multiply) with one vector e per color. For each color c, e[i] = 1 if
  // c(v[i]) == c else 0 (entries = {i : v[i] == c}). This gives us a vector col
  // = sum_i A e_i = sum_i column_i(A). We need to decompose col into its
  // contributions from the different columns. This is possible, as all the
  // columns listed in entries have no non-zero row in common.
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

      // Split the result into its contributions from the single unit vectors,
      // i.e. A*e_i. This is possible, since for each row, only one column in
      // entry can store a value.
      for (vertices_size_type r = 0; r < op.n(); r++)
        {
          if (image[r] != 0.0)
            {
              for (auto c : entries)
                {
                  if (pattern.exists(r, c))
                    {
                      output.set(r, c, image[r]);
#ifdef DEBUG
                      // Verify that we can uniquely identify the column.
                      for (auto j : entries)
                        Assert(c == j || !pattern.exists(r, j),
                               dealii::ExcInternalError());
#endif
                      break;
                    }
                }
            }
        }
    }

  return cache;
}

template <typename Matrix,
          typename Operator,
          typename Pattern,
          typename DomainVector,
          typename RangeVector>
std::shared_ptr<GraphCache>
assemble_operator(Matrix &            output,
                  const Operator &    op,
                  const Pattern &     pattern,
                  const DomainVector &dummy_domain,
                  const RangeVector & dummy_range)
{
  auto cache = std::make_shared<GraphCache>();

  // Early exit for empty pattern.
  if (pattern.n_nonzero_elements() == 0)
    return cache;

  // Create a set with all the edges defined by the sparsity pattern.
  // An edge E = (c1,c2) is contained iff there exists a row r,
  // such that (r,c1) and (r,c2) are both contained in the pattern.
  cache->make_graph(pattern);
  cache->make_coloring();

  AssertDimension(cache->color_vec.size(), pattern.n_rows());

  return assemble_operator(
    output, op, pattern, dummy_domain, dummy_range, cache);
}

template <typename Matrix, typename Operator, typename Pattern>
std::shared_ptr<GraphCache>
assemble_operator(Matrix &                          output,
                  const Operator &                  op,
                  const Pattern &                   pattern,
                  const std::shared_ptr<GraphCache> cache)
{
  dealii::Vector<typename Matrix::value_type> domain(pattern.n_cols());
  dealii::Vector<typename Matrix::value_type> range(pattern.n_rows());
  return assemble_operator(output, op, pattern, domain, range, cache);
}

template <typename Matrix, typename Operator, typename Pattern>
std::shared_ptr<GraphCache>
assemble_operator(Matrix &output, const Operator &op, const Pattern &pattern)
{
  dealii::Vector<typename Matrix::value_type> domain(pattern.n_cols());
  dealii::Vector<typename Matrix::value_type> range(pattern.n_rows());
  return assemble_operator(output, op, pattern, domain, range);
}

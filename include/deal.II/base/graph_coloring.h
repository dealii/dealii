
// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

#ifndef __deal2__graph_coloring_h
#define __deal2__graph_coloring_h


#include <deal.II/base/config.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/std_cxx11/function.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <set>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace containing functions that can color graphs.
 */
namespace GraphColoring
{
  namespace internal
  {
    /**
     * Given two sets of indices that are assumed to be sorted, determine
     * whether they will have a nonempty intersection. The actual intersection is
     * not computed.
     * @param indices1 A set of indices, assumed sorted.
     * @param indices2 A set of indices, assumed sorted.
     * @return Whether the two sets of indices do have a nonempty intersection.
     */
    inline
    bool
    have_nonempty_intersection (const std::vector<types::global_dof_index> &indices1,
                                const std::vector<types::global_dof_index> &indices2)
    {
      // we assume that both arrays are sorted, so we can walk
      // them in lockstep and see if we encounter an index that's
      // in both arrays. once we reach the end of either array,
      // we know that there is no intersection
      std::vector<types::global_dof_index>::const_iterator
      p = indices1.begin(),
      q = indices2.begin();
      while ((p != indices1.end()) && (q != indices2.end()))
        {
          if (*p < *q)
            ++p;
          else if (*p > *q)
            ++q;
          else
            // conflict found!
            return true;
        }

      // no conflict found!
      return false;
    }


    /**
     * Create a partitioning of the given range of iterators using a simplified
     * version of the Cuthill-McKee algorithm (Breadth First Search algorithm).
     * The function creates partitions that contain "zones" of iterators
     * where the first partition contains the first iterator, the second
     * zone contains all those iterators that have conflicts with the single
     * element in the first zone, the third zone contains those iterators that
     * have conflicts with the iterators of the second zone and have not previously
     * been assigned to a zone, etc. If the iterators represent cells, then this
     * generates partitions that are like onion shells around the very first
     * cell. Note that elements in each zone may conflict with other elements in
     * the same zone.
     *
     * The question whether two iterators conflict is determined by a user-provided
     * function. The meaning of this function is discussed in the documentation of
     * the GraphColoring::make_graph_coloring() function.
     *
     * @param[in] begin The first element of a range of iterators for which a
     *      partitioning is sought.
     * @param[in] end The element past the end of the range of iterators.
     * @param[in] get_conflict_indices A user defined function object returning
     *      a set of indicators that are descriptive of what represents a
     *      conflict. See above for a more thorough discussion.
     * @return A set of sets of iterators (where sets are represented by
     *      std::vector for efficiency). Each element of the outermost set
     *      corresponds to the iterators pointing to objects that are in the
     *      same partition (i.e., the same zone).
     *
     * @author Martin Kronbichler, Bruno Turcksin
     */
    template <typename Iterator>
    std::vector<std::vector<Iterator> >
    create_partitioning(const Iterator &begin,
                        const typename identity<Iterator>::type &end,
                        const std_cxx11::function<std::vector<types::global_dof_index> (const Iterator &)> &get_conflict_indices)
    {
      // Number of iterators.
      unsigned int n_iterators = 0;

      // Create a map from conflict indices to iterators
      boost::unordered_map<types::global_dof_index,std::vector<Iterator> > indices_to_iterators;
      for (Iterator it=begin; it!=end; ++it)
        {
          const std::vector<types::global_dof_index> conflict_indices = get_conflict_indices(it);
          const unsigned int n_conflict_indices = conflict_indices.size();
          for (unsigned int i=0; i<n_conflict_indices; ++i)
            indices_to_iterators[conflict_indices[i]].push_back(it);
          ++n_iterators;
        }

      // create the very first zone which contains only the first
      // iterator. then create the other zones. keep track of all the
      // iterators that have already been assigned to a zone
      std::vector<std::vector<Iterator> > zones(1,std::vector<Iterator> (1,begin));
      std::set<Iterator> used_it;
      used_it.insert(begin);
      while (used_it.size()!=n_iterators)
        {
          // loop over the elements of the previous zone. for each element of
          // the previous zone, get the conflict indices and from there get
          // those iterators that are conflicting with the current element
          typename std::vector<Iterator>::iterator previous_zone_it(zones.back().begin());
          typename std::vector<Iterator>::iterator previous_zone_end(zones.back().end());
          std::vector<Iterator> new_zone;
          for (; previous_zone_it!=previous_zone_end; ++previous_zone_it)
            {
              const std::vector<types::global_dof_index>
              conflict_indices = get_conflict_indices(*previous_zone_it);

              const unsigned int n_conflict_indices(conflict_indices.size());
              for (unsigned int i=0; i<n_conflict_indices; ++i)
                {
                  const std::vector<Iterator> &conflicting_elements
                    = indices_to_iterators[conflict_indices[i]];
                  for (unsigned int j=0; j<conflicting_elements.size(); ++j)
                    {
                      // check that the iterator conflicting with the current one is not
                      // associated to a zone yet and if so, assign it to the current
                      // zone. mark it as used
                      //
                      // we can shortcut this test if the conflicting iterator is the
                      // current iterator
                      if ((conflicting_elements[j] != *previous_zone_it)
                          &&
                          (used_it.count(conflicting_elements[j])==0))
                        {
                          new_zone.push_back(conflicting_elements[j]);
                          used_it.insert(conflicting_elements[j]);
                        }
                    }
                }
            }

          // If there are iterators in the new zone, then the zone is added to the
          // partition. Otherwise, the graph is disconnected and we need to find
          // an iterator on the other part of the graph. start the whole process again
          // with the first iterator that hasn't been assigned to a zone yet
          if (new_zone.size()!=0)
            zones.push_back(new_zone);
          else
            for (Iterator it=begin; it!=end; ++it)
              if (used_it.count(it)==0)
                {
                  zones.push_back(std::vector<Iterator> (1,it));
                  used_it.insert(it);
                  break;
                }
        }

      return zones;
    }



    /**
     * This function uses DSATUR (Degree SATURation) to color the elements of
     * a set. DSATUR works as follows:
     *   -# Arrange the vertices by decreasing order of degrees.
     *   -# Color a vertex of maximal degree with color 1.
     *   -# Choose a vertex with a maximal saturation degree. If there is equality,
     *      choose any vertex of maximal degree in the uncolored subgraph.
     *   -# Color the chosen vertex with the least possible (lowest numbered) color.
     *   -# If all the vertices are colored, stop. Otherwise, return to 3.
     *
     * @param[in] partition The set of iterators that should be colored.
     * @param[in] get_conflict_indices A user defined function object returning
     *      a set of indicators that are descriptive of what represents a
     *      conflict. See above for a more thorough discussion.
     * @param[out] partition_coloring A set of sets of iterators (where sets are represented by
     *      std::vector for efficiency). Each element of the outermost set
     *      corresponds to the iterators pointing to objects that are in the
     *      same partition (have the same color) and consequently do not
     *      conflict. The elements of different sets may conflict.
     */
    template <typename Iterator>
    void
    make_dsatur_coloring(std::vector<Iterator> &partition,
                         const std_cxx11::function<std::vector<types::global_dof_index> (const Iterator &)> &get_conflict_indices,
                         std::vector<std::vector<Iterator> > &partition_coloring)
    {
      partition_coloring.clear ();

      // Number of zones composing the partitioning.
      const unsigned int partition_size(partition.size());
      std::vector<unsigned int> sorted_vertices(partition_size);
      std::vector<int> degrees(partition_size);
      std::vector<std::vector<types::global_dof_index> > conflict_indices(partition_size);
      std::vector<std::vector<unsigned int> > graph(partition_size);

      // Get the conflict indices associated to each iterator. The conflict_indices have to
      // be sorted so we can more easily find conflicts later on
      for (unsigned int i=0; i<partition_size; ++i)
        {
          conflict_indices[i] = get_conflict_indices(partition[i]);
          std::sort(conflict_indices[i].begin(), conflict_indices[i].end());
        }

      // Compute the degree of each vertex of the graph using the
      // intersection of the conflict indices.
      for (unsigned int i=0; i<partition_size; ++i)
        for (unsigned int j=i+1; j<partition_size; ++j)
          // If the two iterators share indices then we increase the degree of the
          // vertices and create an ''edge'' in the graph.
          if (have_nonempty_intersection (conflict_indices[i], conflict_indices[j]))
            {
              ++degrees[i];
              ++degrees[j];
              graph[i].push_back(j);
              graph[j].push_back(i);
            }

      // Sort the vertices by decreasing degree.
      std::vector<int>::iterator degrees_it;
      for (unsigned int i=0; i<partition_size; ++i)
        {
          // Find the largest element.
          degrees_it = std::max_element(degrees.begin(),degrees.end());
          sorted_vertices[i] = degrees_it-degrees.begin();
          // Put the largest element to -1 so it cannot be chosen again.
          *degrees_it = -1;
        }

      // Color the graph.
      std::vector<boost::unordered_set<unsigned int> > colors_used;
      for (unsigned int i=0; i<partition_size; ++i)
        {
          const unsigned int current_vertex(sorted_vertices[i]);
          bool new_color(true);
          // Try to use an existing color, i.e., try to find a color which is not
          // associated to one of the vertices linked to current_vertex.
          // Loop over the color.
          for (unsigned int j=0; j<partition_coloring.size(); ++j)
            {
              // Loop on the vertices linked to current_vertex. If one vertex linked
              // to current_vertex is already using the color j, this color cannot
              // be used anymore.
              bool unused_color(true);
              for (unsigned int k=0; k<graph[current_vertex].size(); ++k)
                if (colors_used[j].count(graph[current_vertex][k])==1)
                  {
                    unused_color = false;
                    break;
                  }
              if (unused_color)
                {
                  partition_coloring[j].push_back(partition[current_vertex]);
                  colors_used[j].insert(current_vertex);
                  new_color = false;
                  break;
                }
            }
          // Add a new color.
          if (new_color)
            {
              partition_coloring.push_back(std::vector<Iterator> (1,
                                                                  partition[current_vertex]));
              boost::unordered_set<unsigned int> tmp;
              tmp.insert(current_vertex);
              colors_used.push_back(tmp);
            }
        }
    }



    /**
     * Given a partition-coloring graph, i.e., a set of zones (partitions) each
     * of which is colored, produce a combined coloring for the entire set of
     * iterators. This is possible because any color on an
     * even (resp. odd) zone does not conflict with any color of any other even
     * (resp. odd) zone. Consequently, we can combine colors from all even and all
     * odd zones. This function tries to create colors of similar number of elements.
     */
    template <typename Iterator>
    std::vector<std::vector<Iterator> >
    gather_colors(const std::vector<std::vector<std::vector<Iterator> > > &partition_coloring)
    {
      std::vector<std::vector<Iterator> > coloring;

      // Count the number of iterators in each color.
      const unsigned int partition_size(partition_coloring.size());
      std::vector<std::vector<unsigned int> > colors_counter(partition_size);
      for (unsigned int i=0; i<partition_size; ++i)
        {
          const unsigned int n_colors(partition_coloring[i].size());
          colors_counter[i].resize(n_colors);
          for (unsigned int j=0; j<n_colors; ++j)
            colors_counter[i][j] = partition_coloring[i][j].size();
        }

      // Find the partition with the largest number of colors for the even partition.
      unsigned int i_color(0);
      unsigned int max_even_n_colors(0);
      const unsigned int colors_size(colors_counter.size());
      for (unsigned int i=0; i<colors_size; i+=2)
        {
          if (max_even_n_colors<colors_counter[i].size())
            {
              max_even_n_colors = colors_counter[i].size();
              i_color = i;
            }
        }
      coloring.resize(max_even_n_colors);
      for (unsigned int j=0; j<colors_counter[i_color].size(); ++j)
        coloring[j] = partition_coloring[i_color][j];

      for (unsigned int i=0; i<partition_size; i+=2)
        {
          if (i!=i_color)
            {
              boost::unordered_set<unsigned int> used_k;
              for (unsigned int j=0; j<colors_counter[i].size(); ++j)
                {
                  // Find the color in the current partition with the largest number of
                  // iterators.
                  std::vector<unsigned int>::iterator it;
                  it = std::max_element(colors_counter[i].begin(),colors_counter[i].end());
                  unsigned int min_iterators(-1);
                  unsigned int pos(0);
                  // Find the color of coloring with the least number of colors among
                  // the colors that have not been used yet.
                  for (unsigned int k=0; k<max_even_n_colors; ++k)
                    if (used_k.count(k)==0)
                      if (colors_counter[i_color][k]<min_iterators)
                        {
                          min_iterators = colors_counter[i_color][k];
                          pos = k;
                        }
                  colors_counter[i_color][pos] += *it;
                  // Concatenate the current color with the existing coloring.
                  coloring[pos].insert(coloring[pos].end(),
                                       partition_coloring[i][it-colors_counter[i].begin()].begin(),
                                       partition_coloring[i][it-colors_counter[i].begin()].end());
                  used_k.insert(pos);
                  // Put the number of iterators to the current color to zero.
                  *it = 0;
                }
            }
        }

      // If there is more than one partition, do the same thing that we did for the even partitions
      // to the odd partitions
      if (partition_size>1)
        {
          unsigned int max_odd_n_colors(0);
          for (unsigned int i=1; i<partition_size; i+=2)
            {
              if (max_odd_n_colors<colors_counter[i].size())
                {
                  max_odd_n_colors = colors_counter[i].size();
                  i_color = i;
                }
            }
          coloring.resize(max_even_n_colors+max_odd_n_colors);
          for (unsigned int j=0; j<colors_counter[i_color].size(); ++j)
            coloring[max_even_n_colors+j] = partition_coloring[i_color][j];

          for (unsigned int i=1; i<partition_size; i+=2)
            {
              if (i!=i_color)
                {
                  boost::unordered_set<unsigned int> used_k;
                  for (unsigned int j=0; j<colors_counter[i].size(); ++j)
                    {
                      // Find the color in the current partition with the largest number of
                      // iterators.
                      std::vector<unsigned int>::iterator it;
                      it = std::max_element(colors_counter[i].begin(),colors_counter[i].end());
                      unsigned int min_iterators(-1);
                      unsigned int pos(0);
                      // Find the color of coloring with the least number of colors among
                      // the colors that have not been used yet.
                      for (unsigned int k=0; k<max_odd_n_colors; ++k)
                        if (used_k.count(k)==0)
                          if (colors_counter[i_color][k]<min_iterators)
                            {
                              min_iterators = colors_counter[i_color][k];
                              pos = k;
                            }
                      colors_counter[i_color][pos] += *it;
                      // Concatenate the current color with the existing coloring.
                      coloring[max_even_n_colors+pos].insert(coloring[max_even_n_colors+pos].end(),
                                                             partition_coloring[i][it-colors_counter[i].begin()].begin(),
                                                             partition_coloring[i][it-colors_counter[i].begin()].end());
                      used_k.insert(pos);
                      // Put the number of iterators to the current color to zero.
                      *it = 0;
                    }
                }
            }
        }

      return coloring;
    }
  }


  /**
   * Create a partitioning of the given range of iterators so that
   * iterators that point to conflicting objects will be placed
   * into different partitions, where the question whether two objects conflict
   * is determined by a user-provided function.
   *
   * This function can also be considered as a graph coloring: each object
   * pointed to by an iterator is considered to be a node and there is an
   * edge between each two nodes that conflict. The graph coloring algorithm
   * then assigns a color to each node in such a way that two nodes connected
   * by an edge do not have the same color.
   *
   * A typical use case for this function is in assembling a matrix in parallel.
   * There, one would like to assemble local contributions on different cells
   * at the same time (an operation that is purely local and so requires
   * no synchronization) but then we need to add these local contributions
   * to the global matrix. In general, the contributions from different cells
   * may be to the same matrix entries if the cells share degrees of freedom
   * and, consequently, can not happen at the same time unless we want to
   * risk a race condition (see http://en.wikipedia.org/wiki/Race_condition ).
   * Thus, we call these two cells in conflict, and we can only allow operations
   * in parallel from cells that do not conflict. In other words, two cells
   * are in conflict if the set of matrix entries (for example characterized
   * by the rows) have a nonempty intersection.
   *
   * In this generality, computing the graph of conflicts would require calling
   * a function that determines whether two iterators (or the two objects they
   * represent) conflict, and calling it for every pair of iterators, i.e.,
   * $\frac 12 N (N-1)$ times. This is too expensive in general. A better
   * approach is to require a user-defined function that returns for every
   * iterator it is called for a set of indicators of some kind that characterize
   * a conflict; two iterators are in conflict if their conflict indicator sets
   * have a nonempty intersection. In the example of assembling a matrix,
   * the conflict indicator set would contain the indices of all degrees of
   * freedom on the cell pointed to (in the case of continuous Galerkin methods)
   * or the union of indices of degree of freedom on the current cell and all
   * cells adjacent to the faces of the current cell (in the case of
   * discontinuous Galerkin methods, because there one computes face integrals
   * coupling the degrees of freedom connected by a common face -- see step-12).
   *
   * @note The conflict set returned by the user defined function passed as
   * third argument needs to accurately describe <i>all</i> degrees of
   * freedom for which anything is written into the matrix or right hand side.
   * In other words, if the writing happens through a function like
   * ConstraintMatrix::copy_local_to_global(), then the set of conflict indices
   * must actually contain not only the degrees of freedom on the current
   * cell, but also those they are linked to by constraints such as hanging
   * nodes.
   *
   * In other situations, the conflict indicator sets may represent
   * something different altogether -- it is up to the caller of this function
   * to describe what it means for two iterators to conflict. Given this,
   * computing conflict graph edges can be done significantly more cheaply
   * than with ${\cal O}(N^2)$ operations.
   *
   * In any case, the result of the function will be so that iterators whose
   * conflict indicator sets have overlap will not be assigned to the same
   * color.
   *
   * @note The algorithm used in this function is described in a paper by
   * Turcksin, Kronbichler and Bangerth, see @ref workstream_paper .
   *
   * @param[in] begin The first element of a range of iterators for which a
   *      coloring is sought.
   * @param[in] end The element past the end of the range of iterators.
   * @param[in] get_conflict_indices A user defined function object returning
   *      a set of indicators that are descriptive of what represents a
   *      conflict. See above for a more thorough discussion.
   * @return A set of sets of iterators (where sets are represented by
   *      std::vector for efficiency). Each element of the outermost set
   *      corresponds to the iterators pointing to objects that are in the
   *      same partition (have the same color) and consequently do not
   *      conflict. The elements of different sets may conflict.
   *
   * @author Martin Kronbichler, Bruno Turcksin
   */
  template <typename Iterator>
  std::vector<std::vector<Iterator> >
  make_graph_coloring(const Iterator &begin,
                      const typename identity<Iterator>::type &end,
                      const std_cxx11::function<std::vector<types::global_dof_index> (const typename identity<Iterator>::type &)> &get_conflict_indices)
  {
    Assert (begin != end, ExcMessage ("GraphColoring is not prepared to deal with empty ranges!"));

    // Create the partitioning.
    std::vector<std::vector<Iterator> >
    partitioning = internal::create_partitioning (begin,
                                                  end,
                                                  get_conflict_indices);

    // Color the iterators within each partition.
    // Run the coloring algorithm on each zone in parallel
    const unsigned int partitioning_size(partitioning.size());
    std::vector<std::vector<std::vector<Iterator> > >
    partition_coloring(partitioning_size);

    Threads::TaskGroup<> tasks;
    for (unsigned int i=0; i<partitioning_size; ++i)
      tasks += Threads::new_task (&internal::make_dsatur_coloring<Iterator>,
                                  partitioning[i],
                                  get_conflict_indices,
                                  partition_coloring[i]);
    tasks.join_all();

    // Gather the colors together.
    return internal::gather_colors(partition_coloring);
  }

} // End graph_coloring namespace

DEAL_II_NAMESPACE_CLOSE


//----------------------------   graph_coloring.h     ---------------------------
// end of #ifndef __deal2__graph_coloring_h
#endif
//----------------------------   graph_coloring.h     ---------------------------

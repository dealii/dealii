
// ---------------------------------------------------------------------
// $Id: graph_coloring.h 30494 2013-08-26 10:04:44Z kronbichler $
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
#include <deal.II/dofs/dof_handler.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Currently graph_coloring only works only for continuous FE.
namespace graph_coloring {

  /**
   * Create the partitioning using a simplified version of the Cuthill-McKee
   * algorithm (Breadth First Search algorithm).
   */
  template <typename Iterator>
  std::vector<std::vector<Iterator> > create_partitioning(Iterator const &begin,
      typename identity<Iterator>::type const &end)
  {
    std::vector<std::vector<Iterator> > partitioning(1,std::vector<Iterator> (1,begin));
    // Create a map from dofs to cells
    boost::unordered_map<types::global_dof_index,std::vector<Iterator> > dof_to_cell;
    for (Iterator it=begin; it!=end; ++it)
    {
      // Need to check that it works with hp
      const unsigned int dofs_per_cell(it->get_fe().dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      it->get_dof_indices(local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        dof_to_cell[local_dof_indices[i]].push_back(it);
    }

    bool done(false);
    // TO CHANGE: std::set is slow
    std::set<Iterator> used_it;
    used_it.insert(begin);
    while (!done)
    {
      typename std::vector<Iterator>::iterator vector_it(partitioning.back().begin());
      typename std::vector<Iterator>::iterator vector_end(partitioning.back().end());
      std::vector<Iterator> new_zone;
      for (; vector_it!=vector_end; ++vector_it)
      {
        const unsigned int dofs_per_cell((*vector_it)->get_fe().dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        (*vector_it)->get_dof_indices(local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          std::vector<Iterator> iterator_vector(dof_to_cell[local_dof_indices[i]]);
          for (unsigned int j=0; j<iterator_vector.size(); ++j)
          {
            // Check that the iterator is not associated to a zone yet.
            if (used_it.count(iterator_vector[j])==0)
            {
              new_zone.push_back(iterator_vector[j]);
              used_it.insert(iterator_vector[j]);
            }
          }
        }
      }
      // If there are iterators in the new zone, then the zone is added to the
      // partition. Otherwise, we are done and exit the while.
      if (new_zone.size()!=0)
        partitioning.push_back(new_zone);
      else
        done = true;
    }

    return partitioning;
  }

  /**
   * This function uses DSATUR (Degree SATURation) to color the graph:
   *   -# Arrange the vertices by decreasing order of degrees.
   *   -# Color a vertex of maximal degree with color 1.
   *   -# Choose a vertex with a maximal saturation degree. If there is equality,
   *      choose any vertex of maximal degree in the uncolored subgraph.
   *   -# Color the chosen vertex with the least possible (lowest numbered) color.
   *   -# If all the vertices are colored, stop. Otherwise, return to 3.
   */
  template <typename Iterator>
  std::vector<std::vector<Iterator> > make_dsatur_coloring(std::vector<Iterator> const &partition)
  {
    std::vector<std::vector<Iterator> > partition_coloring;
    // Number of zones composing the partitioning.
    const unsigned int partition_size(partition.size());
    std::vector<unsigned int> sorted_vertices(partition_size);
    std::vector<unsigned int> degrees(partition_size);
    std::vector<std::vector<types::global_dof_index> > local_dof_indices(partition_size);
    std::vector<std::vector<unsigned int> > graph(partition_size);

    // Get the dofs associated to each cell. The dof_indices have to be sorted so
    // set_intersection can be used later.
    for (unsigned int i=0; i<partition_size; ++i)
    {
      local_dof_indices[i].resize(partition[i]->get_fe().dofs_per_cell);
      partition[i]->get_dof_indices(local_dof_indices[i]);
      std::sort(local_dof_indices[i].begin(),local_dof_indices[i].end());
    }

    // Compute the degree of each vertex of the graph (cell) using the
    // intersection of dofs
    std::vector<types::global_dof_index> dofs_intersection;
    std::vector<types::global_dof_index>::iterator intersection_it;
    for (unsigned int i=0; i<partition_size; ++i)
      for (unsigned int j=i+1; j<partition_size; ++j)
      {
        dofs_intersection.resize(std::max(local_dof_indices[i].size(),local_dof_indices[j].size()));
        intersection_it = std::set_intersection(local_dof_indices[i].begin(),
            local_dof_indices[i].end(),local_dof_indices[j].begin(),
            local_dof_indices[j].end(),dofs_intersection.begin());
        // If the two cells share dofs then we increase the degree of the
        // vertices and create an ''edge'' in the graph.
        if (intersection_it!=dofs_intersection.begin())
        {
          ++degrees[i];
          ++degrees[j];
          graph[i].push_back(j);
          graph[j].push_back(i);
        }
      }

    // Sort the vertices by decreasing degree.
    std::vector<unsigned int>::iterator degrees_it;
    for (unsigned int i=0; i<partition_size; ++i)
    {
      // Find the largest element.
      degrees_it = std::max_element(degrees.begin(),degrees.end());
      sorted_vertices[i] = degrees_it-degrees.begin();
      // Zero the largest element.
      *degrees_it = 0;
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

    return partition_coloring;
  }

  /**
   * Given a partition-coloring graph, gather the colors together. All the
   * colors on even (resp. odd) partition can be executed simultaneously. This
   * function tries to create colors of similar number of elements.
   */
  template <typename Iterator>
  std::vector<std::vector<Iterator> >
  gather_colors(std::vector<std::vector<std::vector<Iterator> > > const &partition_coloring)
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
          // cells.
          std::vector<unsigned int>::iterator it;
          it = std::max_element(colors_counter[i].begin(),colors_counter[i].end());
          unsigned int min_cells(-1);
          unsigned int pos(0);
          // Find the color of coloring with the least number of colors among
          // the colors that have not been used yet.
          for (unsigned int k=0; k<max_even_n_colors; ++k)
            if (used_k.count(k)==0)
              if (colors_counter[i_color][k]<min_cells)
              {
                min_cells = colors_counter[i_color][k];
                pos = k;
              }
          colors_counter[i_color][pos] += *it;
          // Concatenate the current color with the existing coloring.
          coloring[pos].insert(coloring[pos].end(),
              partition_coloring[i][it-colors_counter[i].begin()].begin(),
              partition_coloring[i][it-colors_counter[i].begin()].end());
          used_k.insert(pos);
          // Put the number of cells to the current color to zero.
          *it = 0;
        }
      }
    }

    // Do the same thing that we did for the even partitions to the odd
    // partitions
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
          // cells.
          std::vector<unsigned int>::iterator it;
          it = std::max_element(colors_counter[i].begin(),colors_counter[i].end());
          unsigned int min_cells(-1);
          unsigned int pos(0);
          // Find the color of coloring with the least number of colors among
          // the colors that have not been used yet.
          for (unsigned int k=0; k<max_odd_n_colors; ++k)
            if (used_k.count(k)==0)
              if (colors_counter[i_color][k]<min_cells)
              {
                min_cells = colors_counter[i_color][k];
                pos = k;
              }
          colors_counter[i_color][pos] += *it;
          // Concatenate the current color with the existing coloring.
          coloring[max_even_n_colors+pos].insert(coloring[max_even_n_colors+pos].end(),
              partition_coloring[i][it-colors_counter[i].begin()].begin(),
              partition_coloring[i][it-colors_counter[i].begin()].end());
          used_k.insert(pos);
          // Put the number of cells to the current color to zero.
          *it = 0;
        }
      }
    }

    return coloring;
  }

  /**
   * This function creates a coloring given two iterators on the DoFHandler.
   */
  template <typename Iterator>
  std::vector<std::vector<Iterator> > 
  make_graph_coloring(Iterator const &begin,typename identity<Iterator>::type const &end) 
  {
    // Create the partitioning.
    std::vector<std::vector<Iterator> > partitioning = create_partitioning(begin,end);

    // Color the cells within each partition.
    const unsigned int partitioning_size(partitioning.size());
    std::vector<std::vector<std::vector<Iterator> > > partition_coloring(
        partitioning_size);
    for (unsigned int i=0; i<partitioning_size; ++i)
    {
      // Compute the coloring of the graph using the DSATUR algorithm
      partition_coloring[i] = make_dsatur_coloring(partitioning[i]);
    }

    // Gather the colors together. 
    std::vector<std::vector<Iterator> > coloring = gather_colors(partition_coloring);

    return coloring;
  }

} // End graph_coloring namespace

DEAL_II_NAMESPACE_CLOSE


//----------------------------   graph_coloring.h     ---------------------------
// end of #ifndef __deal2__graph_coloring_h
#endif
//----------------------------   graph_coloring.h     ---------------------------

// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



#include <deal.II/base/bounding_box.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_tools.h>

#ifdef DEAL_II_WITH_P4EST

#include <vector>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace GridTools
  {
    template < int dim, int spacedim >
    BoundingBox < spacedim >
    compute_cell_locally_owned_bounding_box
    (const typename parallel::distributed::Triangulation< dim, spacedim >::cell_iterator &parent_cell,
     bool &has_locally_owned)
    {
      //Finding all active cells descendants of the current one (or the current one if it is active)
      std::vector< typename parallel::distributed::Triangulation< dim, spacedim >::active_cell_iterator >
      local_cells = dealii::GridTools::get_active_child_cells < parallel::distributed::Triangulation< dim, spacedim > >
                    (parent_cell);

      //If the current cell is active but not locally owned
      //or none of its children is active:
      if (local_cells.size()==0)
        {
          has_locally_owned = false;
          BoundingBox < spacedim > bbox;
          return bbox;
        }

      Point<spacedim> minp;
      Point<spacedim> maxp;

      bool flag_first = true;
      for (unsigned int i=0; i<local_cells.size(); ++i)
        if ( local_cells[i]->is_locally_owned())
          {
            if (flag_first)
              {
                minp = local_cells[i]->vertex(0);
                maxp = local_cells[i]->vertex(0);
                flag_first = false;
              }

            for (unsigned int v=0; v<GeometryInfo<spacedim>::vertices_per_cell; ++v)
              for ( unsigned int d=0; d<spacedim; ++d)
                {
                  minp[d] = std::min( minp[d], local_cells[i]->vertex(v)[d]);
                  maxp[d] = std::max( maxp[d], local_cells[i]->vertex(v)[d]);
                }
          }

      if (minp != maxp )
        {
          has_locally_owned = true;
          BoundingBox < spacedim > bbox(std::make_pair(minp,maxp));
          return bbox;
        }
      else
        {
          has_locally_owned = false;
          //If none of the local cells is locally owned then
          //minp and maxp are never updated:
          BoundingBox < spacedim > bbox;
          return bbox;
        }
    }


    template < int dim, int spacedim>
    std::vector< BoundingBox<spacedim> >
    compute_locally_owned_bounding_box
    (const parallel::distributed::Triangulation< dim, spacedim > &distributed_tria,
     const unsigned int &refinement_level, const unsigned int &max_boxes, const double &err)
    {
      // Algorithm brief description: we begin by using refinement_lefel (and coarser levels)
      // to create Bounding Boxes of locally owned cells. Then these cells are merged.

      Assert( refinement_level < distributed_tria.n_levels(),
              ExcMessage ( "Error: refinement level is higher then total levels in the triangulation!") );

      std::vector< BoundingBox < spacedim > > bounding_boxes;

      // Creating a bounding box for all active cell on coarser level: the refinement_level should be
      // low enough to avoid the creation  of Bounding Boxes here
      for (unsigned int i=0; i < refinement_level; ++i)
        for (typename Triangulation< dim, spacedim >::cell_iterator
             cell: distributed_tria.active_cell_iterators_on_level(i))
          {
            bool has_locally_owned = false;
            BoundingBox < spacedim > bbox = compute_cell_locally_owned_bounding_box <dim,spacedim> (cell,has_locally_owned);
            if (has_locally_owned)
              bounding_boxes.push_back(bbox);
          }

      // Creating a Bounding Box for all cells on the chosen refinement_level
      for (typename Triangulation< dim, spacedim >::cell_iterator
           cell: distributed_tria.cell_iterators_on_level(refinement_level))
        {
          bool has_locally_owned = false;
          BoundingBox < spacedim > bbox = compute_cell_locally_owned_bounding_box <dim,spacedim> (cell,has_locally_owned);
          if (has_locally_owned)
            bounding_boxes.push_back(bbox);
        }


      // Part 1: merging neighbours
      bool found_neighbours = true;
      // This array stores the indices of arrays we have already merged
      std::vector<unsigned int> merged_boxes;

      while (found_neighbours)
        {
          found_neighbours = false;
          for (unsigned int i=0; i<bounding_boxes.size()-1; ++i)
            {
              if ( std::find(merged_boxes.begin(),merged_boxes.end(),i) == merged_boxes.end())
                for (unsigned int j=i+1; j<bounding_boxes.size(); ++j)
                  if ( std::find(merged_boxes.begin(),merged_boxes.end(),j) == merged_boxes.end()
                       && bounding_boxes[i].is_neighbour(bounding_boxes[j], err) == 2  )
                    {
                      bounding_boxes[i].merge_with(bounding_boxes[j]);
                      merged_boxes.push_back(j);
                      found_neighbours = true;
                    }
            }
        }

      // Part 2: if there are too many bounding boxes, merging smaller boxes
      std::vector< BoundingBox < spacedim > > merged_b_boxes;
      for (unsigned int i=0; i<bounding_boxes.size(); ++i)
        if (std::find(merged_boxes.begin(),merged_boxes.end(),i) == merged_boxes.end())
          merged_b_boxes.push_back(bounding_boxes[i]);



      if (merged_b_boxes.size() > max_boxes)
        {
          std::vector<double> volumes;
          for ( BoundingBox<spacedim> box: merged_b_boxes)
            volumes.push_back(box.volume());

          while ( merged_b_boxes.size() > max_boxes)
            {
              unsigned int min_idx = std::min_element(volumes.begin(),volumes.end()) -
                                     volumes.begin();
              volumes.erase(volumes.begin() + min_idx);
              //Finding a neighbour
              bool not_removed = true;
              for (unsigned int i=0; i<merged_b_boxes.size() && not_removed; ++i)
                if ( i != min_idx && merged_b_boxes[i].is_neighbour(merged_b_boxes[min_idx]) > 0 )
                  {
                    merged_b_boxes[i].merge_with(merged_b_boxes[min_idx]);
                    merged_b_boxes.erase(merged_b_boxes.begin() + min_idx);
                    not_removed = false;
                  }
              Assert( !not_removed,
                      ExcMessage ( "Error: couldn't reach target number of bounding boxes!") );
            }
        }

      return merged_b_boxes;
    }
  }
}

#include "grid_tools.inst"

DEAL_II_NAMESPACE_CLOSE

#endif

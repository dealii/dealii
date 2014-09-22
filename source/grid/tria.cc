// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/table.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/std_cxx11/bind.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_levels.h>
#include <deal.II/grid/tria_faces.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/magic_numbers.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <numeric>
#include <map>
#include <list>
#include <cmath>
#include <functional>

#include <deal.II/base/std_cxx11/array.h>

DEAL_II_NAMESPACE_OPEN

bool
SubCellData::check_consistency (const unsigned int dim) const
{
  switch (dim)
    {
    case 1:
      return ((boundary_lines.size() == 0) &&
              (boundary_quads.size() == 0));
    case 2:
      return (boundary_quads.size() == 0);
    };
  return true;
}


namespace internal
{
  namespace Triangulation
  {

    NumberCache<1>::NumberCache ()
      :
      n_levels (0),
      n_lines (0),
      n_active_lines (0)
      // all other fields are
      // default constructed
    {}



    std::size_t
    NumberCache<1>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (n_levels) +
              MemoryConsumption::memory_consumption (n_lines) +
              MemoryConsumption::memory_consumption (n_lines_level) +
              MemoryConsumption::memory_consumption (n_active_lines) +
              MemoryConsumption::memory_consumption (n_active_lines_level));
    }


    NumberCache<2>::NumberCache () :
      n_quads (0),
      n_active_quads (0)
      // all other fields are
      // default constructed
    {}



    std::size_t
    NumberCache<2>::memory_consumption () const
    {
      return (NumberCache<1>::memory_consumption () +
              MemoryConsumption::memory_consumption (n_quads) +
              MemoryConsumption::memory_consumption (n_quads_level) +
              MemoryConsumption::memory_consumption (n_active_quads) +
              MemoryConsumption::memory_consumption (n_active_quads_level));
    }



    NumberCache<3>::NumberCache () :
      n_hexes (0),
      n_active_hexes (0)
      // all other fields are
      // default constructed
    {}



    std::size_t
    NumberCache<3>::memory_consumption () const
    {
      return (NumberCache<2>::memory_consumption () +
              MemoryConsumption::memory_consumption (n_hexes) +
              MemoryConsumption::memory_consumption (n_hexes_level) +
              MemoryConsumption::memory_consumption (n_active_hexes) +
              MemoryConsumption::memory_consumption (n_active_hexes_level));
    }
  }
}

// anonymous namespace for internal helper functions
namespace
{
  // return whether the given cell is
  // patch_level_1, i.e. determine
  // whether either all or none of
  // its children are further
  // refined. this function can only
  // be called for non-active cells.
  template <int dim, int spacedim>
  bool cell_is_patch_level_1 (const TriaIterator<dealii::CellAccessor<dim, spacedim> > &cell)
  {
    Assert (cell->active() == false, ExcInternalError());

    unsigned int n_active_children = 0;
    for (unsigned int i=0; i<cell->n_children(); ++i)
      if (cell->child(i)->active())
        ++n_active_children;

    return (n_active_children == 0) || (n_active_children == cell->n_children());
  }



  // return, whether a given @p cell will be
  // coarsened, which is the case if all
  // children are active and have their coarsen
  // flag set. In case only part of the coarsen
  // flags are set, remove them.
  template <int dim, int spacedim>
  bool cell_will_be_coarsened (const TriaIterator<dealii::CellAccessor<dim,spacedim> > &cell)
  {
    // only cells with children should be
    // considered for coarsening

    if (cell->has_children())
      {
        unsigned int children_to_coarsen=0;
        const unsigned int n_children=cell->n_children();

        for (unsigned int c=0; c<n_children; ++c)
          if (cell->child(c)->active() &&
              cell->child(c)->coarsen_flag_set())
            ++children_to_coarsen;
        if (children_to_coarsen==n_children)
          return true;
        else
          for (unsigned int c=0; c<n_children; ++c)
            if (cell->child(c)->active())
              cell->child(c)->clear_coarsen_flag();
      }
    // no children, so no coarsening
    // possible. however, no children also
    // means that this cell will be in the same
    // state as if it had children and was
    // coarsened. So, what should we return -
    // false or true?
    // make sure we do not have to do this at
    // all...
    Assert(cell->has_children(), ExcInternalError());
    // ... and then simply return false
    return false;
  }


  // return, whether the face @p face_no of the
  // given @p cell will be refined after the
  // current refinement step, considering
  // refine and coarsen flags and considering
  // only those refinemnts that will be caused
  // by the neighboring cell.

  // this function is used on both active cells
  // and cells with children. on cells with
  // children it also of interest to know 'how'
  // the face will be refined. thus there is an
  // additional third argument @p
  // expected_face_ref_case returning just
  // that. be aware, that this vriable will
  // only contain useful information if this
  // function is called for an active cell.
  //
  // thus, this is an internal function, users
  // should call one of the two alternatives
  // following below.
  template <int dim, int spacedim>
  bool
  face_will_be_refined_by_neighbor_internal(const TriaIterator<dealii::CellAccessor<dim,spacedim> > &cell,
                                            const unsigned int                                   face_no,
                                            RefinementCase<dim-1>                                    &expected_face_ref_case)
  {
    // first of all: set the default value for
    // expected_face_ref_case, which is no
    // refinement at all
    expected_face_ref_case=RefinementCase<dim-1>::no_refinement;

    const typename Triangulation<dim,spacedim>::cell_iterator neighbor=cell->neighbor(face_no);

    // If we are at the boundary, there is no
    // neighbor which could refine the face
    if (neighbor.state()!=IteratorState::valid)
      return false;

    if (neighbor->has_children())
      {
        // if the neighbor is refined, it may be
        // coarsened. if so, then it won't refine
        // the face, no matter what else happens
        if (cell_will_be_coarsened(neighbor))
          return false;
        else
          // if the neighor is refined, then he
          // is also refined at our current
          // face. He will stay so without
          // coarsening, so return true in that
          // case.
          {
            expected_face_ref_case=cell->face(face_no)->refinement_case();
            return true;
          }
      }

    // now, the neighbor is not refined, but
    // perhaps he will be
    const RefinementCase<dim> nb_ref_flag=neighbor->refine_flag_set();
    if (nb_ref_flag != RefinementCase<dim>::no_refinement)
      {
        // now we need to know, which of the
        // neighbors faces points towards us
        const unsigned int neighbor_neighbor=cell->neighbor_face_no(face_no);
        // check, whether the cell will be
        // refined in a way that refines our
        // face
        const RefinementCase<dim-1> face_ref_case=
          GeometryInfo<dim>::face_refinement_case(nb_ref_flag,
                                                  neighbor_neighbor,
                                                  neighbor->face_orientation(neighbor_neighbor),
                                                  neighbor->face_flip(neighbor_neighbor),
                                                  neighbor->face_rotation(neighbor_neighbor));
        if (face_ref_case != RefinementCase<dim-1>::no_refinement)
          {
            const typename Triangulation<dim,spacedim>::face_iterator neighbor_face=neighbor->face(neighbor_neighbor);
            const int this_face_index=cell->face_index(face_no);

            // there are still two basic
            // possibilities here: the neighbor
            // might be coarser or as coarse
            // as we are
            if (neighbor_face->index()==this_face_index)
              // the neighbor is as coarse as
              // we are and will be refined at
              // the face of consideration, so
              // return true
              {
                expected_face_ref_case = face_ref_case;
                return true;
              }
            else
              {

                // the neighbor is coarser.
                // this is the most complicated
                // case. It might be, that the
                // neighbor's face will be
                // refined, but that we will
                // not see this, as we are
                // refined in a similar way.

                // so, the neighbor's face must
                // have children. check, if our
                // cell's face is one of these
                // (it could also be a
                // grand_child)
                for (unsigned int c=0; c<neighbor_face->n_children(); ++c)
                  if (neighbor_face->child_index(c)==this_face_index)
                    {
                      // if the flagged refine
                      // case of the face is a
                      // subset or the same as
                      // the current refine case,
                      // then the face, as seen
                      // from our cell, won't be
                      // refined by the neighbor
                      if ((neighbor_face->refinement_case() | face_ref_case)
                          == neighbor_face->refinement_case())
                        return false;
                      else
                        {
                          // if we are active, we
                          // must be an
                          // anisotropic child
                          // and the coming
                          // face_ref_case is
                          // isotropic. Thus,
                          // from our cell we
                          // will see exactly the
                          // opposite refine case
                          // that the face has
                          // now...
                          Assert(face_ref_case==RefinementCase<dim-1>::isotropic_refinement, ExcInternalError());
                          expected_face_ref_case = ~neighbor_face->refinement_case();
                          return true;
                        }
                    }

                // so, obviously we were not
                // one of the children, but a
                // grandchild. This is only
                // possible in 3d.
                Assert(dim==3, ExcInternalError());
                // In that case, however, no
                // matter what the neighbor
                // does, he won't be finer
                // after the next refinement
                // step.
                return false;
              }
          }// if face will be refined
      }// if neighbor is flagged for refinement

    // no cases left, so the neighbor will not
    // refine the face
    return false;
  }

  // version of above function for both active
  // and non-active cells
  template <int dim, int spacedim>
  bool
  face_will_be_refined_by_neighbor(const TriaIterator<dealii::CellAccessor<dim, spacedim> > &cell,
                                   const unsigned int                                   face_no)
  {
    RefinementCase<dim-1> dummy = RefinementCase<dim-1>::no_refinement;
    return face_will_be_refined_by_neighbor_internal(cell, face_no, dummy);
  }

  // version of above function for active cells
  // only. Additionally returning the refine
  // case (to come) of the face under
  // consideration
  template <int dim, int spacedim>
  bool
  face_will_be_refined_by_neighbor(const TriaActiveIterator<dealii::CellAccessor<dim,spacedim> > &cell,
                                   const unsigned int                                         face_no,
                                   RefinementCase<dim-1>                                          &expected_face_ref_case)
  {
    return face_will_be_refined_by_neighbor_internal(cell, face_no,
                                                     expected_face_ref_case);
  }



  template <int dim, int spacedim>
  bool
  satisfies_level1_at_vertex_rule (const Triangulation<dim,spacedim> &triangulation)
  {
    std::vector<unsigned int> min_adjacent_cell_level (triangulation.n_vertices(),
                                                       triangulation.n_levels());
    std::vector<unsigned int> max_adjacent_cell_level (triangulation.n_vertices(),
                                                       0);

    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          min_adjacent_cell_level[cell->vertex_index(v)]
            = std::min<unsigned int>
              (min_adjacent_cell_level[cell->vertex_index(v)],
               cell->level());
          max_adjacent_cell_level[cell->vertex_index(v)]
            = std::max<unsigned int> (min_adjacent_cell_level[cell->vertex_index(v)],
                                      cell->level());
        }

    for (unsigned int k=0; k<triangulation.n_vertices(); ++k)
      if (triangulation.vertex_used(k))
        if (max_adjacent_cell_level[k] -
            min_adjacent_cell_level[k] > 1)
          return false;
    return true;
  }



  /**
   * Fill the vector @p line_cell_count
   * needed by @p delete_children with the
   * number of cells bounded by a given
   * line.
   */
  template <int dim, int spacedim>
  std::vector<unsigned int>
  count_cells_bounded_by_line (const Triangulation<dim,spacedim> &triangulation)
  {
    if (dim >= 2)
      {
        std::vector<unsigned int> line_cell_count(triangulation.n_raw_lines(),0);
        typename Triangulation<dim,spacedim>::cell_iterator
        cell=triangulation.begin(),
        endc=triangulation.end();
        for (; cell!=endc; ++cell)
          for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
            ++line_cell_count[cell->line_index(l)];
        return line_cell_count;
      }
    else
      return std::vector<unsigned int>();
  }




  /**
   * Fill the vector @p quad_cell_count
   * needed by @p delete_children with the
   * number of cells bounded by a given
   * quad.
   */
  template <int dim, int spacedim>
  std::vector<unsigned int>
  count_cells_bounded_by_quad (const Triangulation<dim,spacedim> &triangulation)
  {
    if (dim >= 3)
      {
        std::vector<unsigned int> quad_cell_count (triangulation.n_raw_quads(),0);
        typename Triangulation<dim,spacedim>::cell_iterator
        cell=triangulation.begin(),
        endc=triangulation.end();
        for (; cell!=endc; ++cell)
          for (unsigned int q=0; q<GeometryInfo<dim>::faces_per_cell; ++q)
            ++quad_cell_count[cell->quad_index(q)];
        return quad_cell_count;
      }
    else
      return std::vector<unsigned int>();
  }



  /**
   * A set of three functions that
   * reorder the data given to
   * create_triangulation_compatibility
   * from the "classic" to the
   * "current" format of vertex
   * numbering of cells and
   * faces. These functions do the
   * reordering of their arguments
   * in-place.
   */
  void
  reorder_compatibility (const std::vector<CellData<1> > &,
                         const SubCellData &)
  {
    // nothing to do here: the format
    // hasn't changed for 1d
  }


  void
  reorder_compatibility (std::vector<CellData<2> > &cells,
                         const SubCellData &)
  {
    for (unsigned int cell=0; cell<cells.size(); ++cell)
      std::swap(cells[cell].vertices[2],cells[cell].vertices[3]);
  }


  void
  reorder_compatibility (std::vector<CellData<3> > &cells,
                         SubCellData               &subcelldata)
  {
    unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
    for (unsigned int cell=0; cell<cells.size(); ++cell)
      {
        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
          tmp[i] = cells[cell].vertices[i];
        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
          cells[cell].vertices[GeometryInfo<3>::ucd_to_deal[i]] = tmp[i];
      }

    // now points in boundary quads
    std::vector<CellData<2> >::iterator boundary_quad
      = subcelldata.boundary_quads.begin();
    std::vector<CellData<2> >::iterator end_quad
      = subcelldata.boundary_quads.end();
    for (unsigned int quad_no=0; boundary_quad!=end_quad; ++boundary_quad, ++quad_no)
      std::swap(boundary_quad->vertices[2], boundary_quad->vertices[3]);
  }



  /**
   * Return the index of the vertex
   * in the middle of this object,
   * if it exists. In order to
   * exist, the object needs to be
   * refined - for 2D and 3D it
   * needs to be refined
   * isotropically or else the
   * anisotropic children have to
   * be refined again. If the
   * middle vertex does not exist,
   * return
   * <tt>numbers::invalid_unsigned_int</tt>.
   *
   * This function should not really be
   * used in application programs.
   */
  template <int dim, int spacedim>
  unsigned int
  middle_vertex_index(const typename Triangulation<dim,spacedim>::line_iterator &line)
  {
    if (line->has_children())
      return line->child(0)->vertex_index(1);
    return numbers::invalid_unsigned_int;
  }


  template <int dim, int spacedim>
  unsigned int
  middle_vertex_index(const typename Triangulation<dim,spacedim>::quad_iterator &quad)
  {
    switch (static_cast<unsigned char> (quad->refinement_case()))
      {
      case RefinementCase<2>::cut_x:
        return middle_vertex_index<dim,spacedim>(quad->child(0)->line(1));
        break;
      case RefinementCase<2>::cut_y:
        return middle_vertex_index<dim,spacedim>(quad->child(0)->line(3));
        break;
      case RefinementCase<2>::cut_xy:
        return quad->child(0)->vertex_index(3);
        break;
      default:
        break;
      }
    return numbers::invalid_unsigned_int;
  }


  template <int dim, int spacedim>
  unsigned int
  middle_vertex_index(const typename Triangulation<dim,spacedim>::hex_iterator &hex)
  {
    switch (static_cast<unsigned char> (hex->refinement_case()))
      {
      case RefinementCase<3>::cut_x:
        return middle_vertex_index<dim,spacedim>(hex->child(0)->quad(1));
        break;
      case RefinementCase<3>::cut_y:
        return middle_vertex_index<dim,spacedim>(hex->child(0)->quad(3));
        break;
      case RefinementCase<3>::cut_z:
        return middle_vertex_index<dim,spacedim>(hex->child(0)->quad(5));
        break;
      case RefinementCase<3>::cut_xy:
        return middle_vertex_index<dim,spacedim>(hex->child(0)->line(11));
        break;
      case RefinementCase<3>::cut_xz:
        return middle_vertex_index<dim,spacedim>(hex->child(0)->line(5));
        break;
      case RefinementCase<3>::cut_yz:
        return middle_vertex_index<dim,spacedim>(hex->child(0)->line(7));
        break;
      case RefinementCase<3>::cut_xyz:
        return hex->child(0)->vertex_index(7);
        break;
      default:
        break;
      }
    return numbers::invalid_unsigned_int;
  }


  /**
   * Collect all coarse mesh cells
   * with at least one vertex at
   * which the determinant of the
   * Jacobian is zero or
   * negative. This is the function
   * for the case dim!=spacedim,
   * where we can not determine
   * whether a cell is twisted as it
   * may, for example, discretize a
   * manifold with a twist.
   */
  template <class TRIANGULATION>
  inline
  typename TRIANGULATION::DistortedCellList
  collect_distorted_coarse_cells (const TRIANGULATION &)
  {
    return typename TRIANGULATION::DistortedCellList();
  }



  /**
   * Collect all coarse mesh cells
   * with at least one vertex at
   * which the determinant of the
   * Jacobian is zero or
   * negative. This is the function
   * for the case dim==spacedim.
   */
  template <int dim>
  inline
  typename Triangulation<dim,dim>::DistortedCellList
  collect_distorted_coarse_cells (const Triangulation<dim,dim> &triangulation)
  {
    typename Triangulation<dim,dim>::DistortedCellList distorted_cells;
    for (typename Triangulation<dim,dim>::cell_iterator
         cell = triangulation.begin(0); cell != triangulation.end(0); ++cell)
      {
        Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          vertices[i] = cell->vertex(i);

        Tensor<0,dim> determinants[GeometryInfo<dim>::vertices_per_cell];
        GeometryInfo<dim>::alternating_form_at_vertices (vertices,
                                                         determinants);

        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          if (determinants[i] <= 1e-9 * std::pow (cell->diameter(),
                                                  1.*dim))
            {
              distorted_cells.distorted_cells.push_back (cell);
              break;
            }
      }

    return distorted_cells;
  }


  /**
   * Return whether any of the
   * children of the given cell is
   * distorted or not. This is the
   * function for dim==spacedim.
   */
  template <int dim>
  bool
  has_distorted_children (const typename Triangulation<dim,dim>::cell_iterator &cell,
                          internal::int2type<dim>,
                          internal::int2type<dim>)
  {
    Assert (cell->has_children(), ExcInternalError());

    for (unsigned int c=0; c<cell->n_children(); ++c)
      {
        Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          vertices[i] = cell->child(c)->vertex(i);

        Tensor<0,dim> determinants[GeometryInfo<dim>::vertices_per_cell];
        GeometryInfo<dim>::alternating_form_at_vertices (vertices,
                                                         determinants);

        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          if (determinants[i] <= 1e-9 * std::pow (cell->child(c)->diameter(),
                                                  1.*dim))
            return true;
      }

    return false;
  }


  /**
   * Function for dim!=spacedim. As
   * for
   * collect_distorted_coarse_cells,
   * there is nothing that we can do
   * in this case.
   */
  template <int dim, int spacedim>
  bool
  has_distorted_children (const typename Triangulation<dim,spacedim>::cell_iterator &,
                          internal::int2type<dim>,
                          internal::int2type<spacedim>)
  {
    return false;
  }



  /**
   * For a given triangulation: set up the
   * neighbor information on all cells.
   */
  template <int spacedim>
  void
  update_neighbors (Triangulation<1,spacedim> &)
  {
  }


  template <int dim, int spacedim>
  void
  update_neighbors (Triangulation<dim,spacedim> &triangulation)
  {
    // each face can be neighbored on two sides
    // by cells. according to the face's
    // intrinsic normal we define the left
    // neighbor as the one for which the face
    // normal points outward, and store that
    // one first; the second one is then
    // the right neighbor for which the
    // face normal points inward. This
    // information depends on the type of cell
    // and local number of face for the
    // 'standard ordering and orientation' of
    // faces and then on the face_orientation
    // information for the real mesh. Set up a
    // table to have fast access to those
    // offsets (0 for left and 1 for
    // right). Some of the values are invalid
    // as they reference too large face
    // numbers, but we just leave them at a
    // zero value.
    //
    // Note, that in 2d for lines as faces the
    // normal direction given in the
    // GeometryInfo class is not consistent. We
    // thus define here that the normal for a
    // line points to the right if the line
    // points upwards.
    //
    // There is one more point to
    // consider, however: if we have
    // dim<spacedim, then we may have
    // cases where cells are
    // inverted. In effect, both
    // cells think they are the left
    // neighbor of an edge, for
    // example, which leads us to
    // forget neighborship
    // information (a case that shows
    // this is
    // codim_one/hanging_nodes_02). We
    // store whether a cell is
    // inverted using the
    // direction_flag, so if a cell
    // has a false direction_flag,
    // then we need to invert our
    // selection whether we are a
    // left or right neighbor in all
    // following computations.
    //
    // first index:  dimension (minus 2)
    // second index: local face index
    // third index:  face_orientation (false and true)
    static const unsigned int left_right_offset[2][6][2] =
    {
      // quadrilateral
      { {0,1}, // face 0, face_orientation = false and true
        {1,0}, // face 1, face_orientation = false and true
        {1,0}, // face 2, face_orientation = false and true
        {0,1}, // face 3, face_orientation = false and true
        {0,0}, // face 4, invalid face
        {0,0}
      },// face 5, invalid face
      // hexahedron
      { {0,1},
        {1,0},
        {0,1},
        {1,0},
        {0,1},
        {1,0}
      }
    };

    // now create a vector of the two active
    // neighbors (left and right) for each face
    // and fill it by looping over all cells. For
    // cases with anisotropic refinement and more
    // then one cell neighboring at a given side
    // of the face we will automatically get the
    // active one on the highest level as we loop
    // over cells from lower levels first.
    const typename Triangulation<dim,spacedim>::cell_iterator dummy;
    std::vector<typename Triangulation<dim,spacedim>::cell_iterator>
    adjacent_cells(2*triangulation.n_raw_faces(), dummy);

    typename Triangulation<dim,spacedim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        {
          const typename Triangulation<dim,spacedim>::face_iterator
          face=cell->face(f);

          const unsigned int
          offset = (cell->direction_flag()
                    ?
                    left_right_offset[dim-2][f][cell->face_orientation(f)]
                    :
                    1-left_right_offset[dim-2][f][cell->face_orientation(f)]);

          adjacent_cells[2*face->index() + offset] = cell;

          // if this cell is not refined, but the
          // face is, then we'll have to set our
          // cell as neighbor for the child faces
          // as well. Fortunately the normal
          // orientation of children will be just
          // the same.
          if (dim==2)
            {
              if (cell->active() && face->has_children())
                {
                  adjacent_cells[2*face->child(0)->index() + offset] = cell;
                  adjacent_cells[2*face->child(1)->index() + offset] = cell;
                }
            }
          else // -> dim == 3
            {
              // We need the same as in 2d
              // here. Furthermore, if the face is
              // refined with cut_x or cut_y then
              // those children again in the other
              // direction, and if this cell is
              // refined isotropically (along the
              // face) then the neighbor will
              // (probably) be refined as cut_x or
              // cut_y along the face. For those
              // neighboring children cells, their
              // neighbor will be the current,
              // inactive cell, as our children are
              // too fine to be neighbors. Catch that
              // case by also acting on inactive
              // cells with isotropic refinement
              // along the face. If the situation
              // described is not present, the data
              // will be overwritten later on when we
              // visit cells on finer levels, so no
              // harm will be done.
              if (face->has_children() &&
                  (cell->active() ||
                   GeometryInfo<dim>::face_refinement_case(cell->refinement_case(),f) == RefinementCase<dim-1>::isotropic_refinement))
                {

                  for (unsigned int c=0; c<face->n_children(); ++c)
                    adjacent_cells[2*face->child(c)->index() + offset] = cell;
                  if (face->child(0)->has_children())
                    {
                      adjacent_cells[2*face->child(0)->child(0)->index() + offset] = cell;
                      adjacent_cells[2*face->child(0)->child(1)->index() + offset] = cell;
                    }
                  if (face->child(1)->has_children())
                    {
                      adjacent_cells[2*face->child(1)->child(0)->index() + offset] = cell;
                      adjacent_cells[2*face->child(1)->child(1)->index() + offset] = cell;
                    }
                } // if cell active and face refined
            } // else -> dim==3
        } // for all faces of all cells

    // now loop again over all cells and set the
    // corresponding neighbor cell. Note, that we
    // have to use the opposite of the
    // left_right_offset in this case as we want
    // the offset of the neighbor, not our own.
    for (cell=triangulation.begin(); cell != endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        {
          const unsigned int
          offset = (cell->direction_flag()
                    ?
                    left_right_offset[dim-2][f][cell->face_orientation(f)]
                    :
                    1-left_right_offset[dim-2][f][cell->face_orientation(f)]);
          cell->set_neighbor(f,
                             adjacent_cells[2*cell->face(f)->index() + 1 - offset]);
        }
  }

}// end of anonymous namespace


namespace internal
{
  namespace Triangulation
  {
    // make sure that if in the following we
    // write Triangulation<dim,spacedim>
    // we mean the *class*
    // dealii::Triangulation, not the
    // enclosing namespace
    // internal::Triangulation
    using dealii::Triangulation;

    /**
     *  Exception
     * @ingroup Exceptions
     */
    DeclException0 (ExcCellShouldBeUnused);
    /**
     *  Exception
     * @ingroup Exceptions
     */
    DeclException0 (ExcTooFewVerticesAllocated);
    /**
     *  Exception
     * @ingroup Exceptions
     */
    DeclException0 (ExcUncaughtState);
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException2 (ExcGridsDoNotMatch,
                    int, int,
                    << "The present grid has " << arg1 << " active cells, "
                    << "but the one in the file had " << arg2);
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException1 (ExcGridHasInvalidCell,
                    int,
                    << "Something went wrong when making cell " << arg1
                    << ". Read the docs and the source code "
                    << "for more information.");
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException0 (ExcGridHasInvalidVertices);
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException1 (ExcInternalErrorOnCell,
                    int,
                    << "Something went wrong upon construction of cell "
                    << arg1);
    /**
     * A cell was entered which has
     * negative measure. In most
     * cases, this is due to a wrong
     * order of the vertices of the
     * cell.
     *
     * @ingroup Exceptions
     */
    DeclException1 (ExcCellHasNegativeMeasure,
                    int,
                    << "Cell " << arg1 << " has negative measure.");
    /**
     * A cell is created with a
     * vertex number exceeding the
     * vertex array.
     *
     * @ingroup Exceptions
     */
    DeclException3 (ExcInvalidVertexIndex,
                    int, int, int,
                    << "Error while creating cell " << arg1
                    << ": the vertex index " << arg2 << " must be between 0 and "
                    << arg3 << ".");
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException2 (ExcLineInexistant,
                    int, int,
                    << "When trying to give a boundary indicator to a line: "
                    << "the line with end vertices " << arg1 << " and "
                    << arg2 << " does not exist.");
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException4 (ExcQuadInexistant,
                    int, int, int, int,
                    << "When trying to give a boundary indicator to a quad: "
                    << "the quad with bounding lines " << arg1 << ", " << arg2
                    << ", " << arg3 << ", " << arg4 << " does not exist.");
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException0 (ExcInteriorLineCantBeBoundary);
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException0 (ExcInteriorQuadCantBeBoundary);
    /**
     * Exception
     * @ingroup Exceptions
     */
    DeclException2 (ExcMultiplySetLineInfoOfLine,
                    int, int,
                    << "In SubCellData the line info of the line with vertex indices "
                    << arg1 << " and " << arg2 << " is multiply set.");


    /**
     * A class into which we put many of the functions that implement
     * functionality of the Triangulation class. The main reason for this
     * class is as follows: the majority of the functions in Triangulation
     * need to be implemented differently for dim==1, dim==2, and
     * dim==3. However, their implementation is largly independent of the
     * spacedim template parameter. So we would like to write things like
     *
     * template <int spacedim>
     * void Triangulation<1,spacedim>::create_triangulation (...) {...}
     *
     * Unfortunately, C++ doesn't allow this: member functions of class
     * templates have to be either not specialized at all, or fully
     * specialized. No partial specialization is allowed. One possible
     * solution would be to just duplicate the bodies of the functions and
     * have equally implemented functions
     *
     * template <>
     * void Triangulation<1,1>::create_triangulation (...) {...}
     *
     * template <>
     * void Triangulation<1,2>::create_triangulation (...) {...}
     *
     * but that is clearly an unsatisfactory solution. Rather, what we do
     * is introduce the current Implementation class in which we can write
     * these functions as member templates over spacedim, i.e. we can have
     *
     * template <int dim_, int spacedim_>
     * template <int spacedim>
     * void Triangulation<dim_,spacedim_>::Implementation::
     *            create_triangulation (...,
     *                                  Triangulation<1,spacedim> &tria ) {...}
     *
     * The outer template parameters are here unused, only the inner
     * ones are of real interest.
     *
     * One may ask why we put these functions into an class rather
     * than an anonymous namespace, for example?
     *
     * First, these implementation functions need to be friends of the
     * Triangulation class. It is simpler to make the entire class a friend
     * rather than listing all members of an implementation namespace as
     * friends of the Triangulation class (there is no such thing as a "friend
     * namespace XXX" directive).
     *
     * Ideally, we would make this class a member class of the
     * Triangulation<dim,spacedim> class, since then our implementation functions
     * have immediate access to the typedefs and static functions of the
     * surrounding Triangulation class. I.e., we do not have to write "typename
     * Triangulation<dim,spacedim>::active_cell_iterator" but can write
     * "active_cell_iterator" right away. This is, in fact, the way it was
     * implemented first, but we ran into a bug in gcc4.0:
     * @code
     *  class Triangulation {
     *    struct Implementation;
     *    friend class TriaAccessor;
     *  };
     *
     *  class TriaAccessor {
     *    struct Implementation;
     *    friend class Triangulation;
     *  };
     * @endcode
     *
     * Here, friendship (per C++ standard) is supposed to extend to all members of
     * the befriended class, including its 'Implementation' member class. But gcc4.0
     * gets this wrong: the members of Triangulation::Implementation are not friends
     * of TriaAccessor and the other way around. Ideally, one would fix this by
     * saying
     * @code
     *  class Triangulation {
     *    struct Implementation;
     *    friend class TriaAccessor;
     *    friend class TriaAccessor::Implementation;   // **
     *  };
     *
     *  class TriaAccessor {
     *    struct Implementation;
     *    friend class Triangulation;
     *    friend class Triangulation::Implementation;
     *  };
     * @endcode
     * but that's not legal because in ** we don't know yet that TriaAccessor has
     * a member class Implementation and so we can't make it a friend. The only
     * way forward at this point was to make Implementation a class in the
     * internal namespace so that we can forward declare it and make it a friend
     * of the respective other outer class -- not quite what we wanted but the
     * only way I could see to make it work...
     */
    struct Implementation
    {
      /**
       * For a given Triangulation, update the
       * number cache for lines. For 1d, we have
       * to deal with the fact that lines have
       * levels, whereas for higher dimensions
       * they do not.
       *
       * The second argument indicates
       * for how many levels the
       * Triangulation has objects,
       * though the highest levels need
       * not contain active cells if they
       * have previously all been
       * coarsened away.
       */
      template <int dim, int spacedim>
      static
      void compute_number_cache (const Triangulation<dim,spacedim>       &triangulation,
                                 const unsigned int                       level_objects,
                                 internal::Triangulation::NumberCache<1> &number_cache)
      {
        typedef
        typename Triangulation<dim,spacedim>::line_iterator line_iterator;
        typedef
        typename Triangulation<dim,spacedim>::active_line_iterator active_line_iterator;

        number_cache.n_levels = 0;
        if (level_objects > 0)
          // find the last level
          // on which there are
          // used cells
          for (unsigned int level=0; level<level_objects; ++level)
            if (triangulation.begin(level) !=
                triangulation.end(level))
              number_cache.n_levels = level+1;

        // no cells at all?
        Assert (number_cache.n_levels > 0, ExcInternalError());

        ///////////////////////////////////
        // update the number of lines
        // on the different levels in
        // the cache
        number_cache.n_lines_level.resize (number_cache.n_levels);
        number_cache.n_lines = 0;

        number_cache.n_active_lines_level.resize (number_cache.n_levels);
        number_cache.n_active_lines = 0;

        // for 1d, lines have levels so take
        // count the objects per level and
        // globally
        if (dim == 1)
          {
            for (unsigned int level=0; level<number_cache.n_levels; ++level)
              {
                // count lines on this level
                number_cache.n_lines_level[level] = 0;

                line_iterator line = triangulation.begin_line (level),
                              endc = (level == number_cache.n_levels-1 ?
                                      line_iterator(triangulation.end_line()) :
                                      triangulation.begin_line (level+1));
                for (; line!=endc; ++line)
                  ++number_cache.n_lines_level[level];

                // update total number of lines
                number_cache.n_lines += number_cache.n_lines_level[level];
              }

            // do the update for the number of
            // active lines as well
            for (unsigned int level=0; level<number_cache.n_levels; ++level)
              {
                // count lines on this level
                number_cache.n_active_lines_level[level] = 0;

                active_line_iterator line = triangulation.begin_active_line (level),
                                     endc = triangulation.end_line ();
                for (; (line!=endc) && (line->level() == static_cast<signed int>(level)); ++line)
                  ++number_cache.n_active_lines_level[level];

                // update total number of lines
                number_cache.n_active_lines += number_cache.n_active_lines_level[level];
              }
          }
        else
          {
            // for dim>1, there are no
            // levels for lines
            {
              line_iterator line = triangulation.begin_line (),
                            endc = triangulation.end_line();
              for (; line!=endc; ++line)
                ++number_cache.n_lines;
            }

            {
              active_line_iterator line = triangulation.begin_active_line (),
                                   endc = triangulation.end_line();
              for (; line!=endc; ++line)
                ++number_cache.n_active_lines;
            }
          }
      }

      /**
       * For a given Triangulation, update the
       * number cache for quads. For 2d, we have
       * to deal with the fact that quads have
       * levels, whereas for higher dimensions
       * they do not.
       *
       * The second argument indicates
       * for how many levels the
       * Triangulation has objects,
       * though the highest levels need
       * not contain active cells if they
       * have previously all been
       * coarsened away.
       *
       * At the beginning of the function, we call the
       * respective function to update the number
       * cache for lines.
       */
      template <int dim, int spacedim>
      static
      void compute_number_cache (const Triangulation<dim,spacedim>       &triangulation,
                                 const unsigned int                       level_objects,
                                 internal::Triangulation::NumberCache<2> &number_cache)
      {
        // update lines and n_levels
        compute_number_cache (triangulation,
                              level_objects,
                              static_cast<internal::Triangulation::NumberCache<1>&>
                              (number_cache));

        typedef
        typename Triangulation<dim,spacedim>::quad_iterator quad_iterator;
        typedef
        typename Triangulation<dim,spacedim>::active_quad_iterator active_quad_iterator;

        ///////////////////////////////////
        // update the number of quads
        // on the different levels in
        // the cache
        number_cache.n_quads_level.resize (number_cache.n_levels);
        number_cache.n_quads = 0;

        number_cache.n_active_quads_level.resize (number_cache.n_levels);
        number_cache.n_active_quads = 0;

        // for 2d, quads have levels so take
        // count the objects per level and
        // globally
        if (dim == 2)
          {
            for (unsigned int level=0; level<number_cache.n_levels; ++level)
              {
                // count quads on this level
                number_cache.n_quads_level[level] = 0;

                quad_iterator quad = triangulation.begin_quad (level),
                              endc = (level == number_cache.n_levels-1 ?
                                      quad_iterator(triangulation.end_quad()) :
                                      triangulation.begin_quad (level+1));
                for (; quad!=endc; ++quad)
                  ++number_cache.n_quads_level[level];

                // update total number of quads
                number_cache.n_quads += number_cache.n_quads_level[level];
              }

            // do the update for the number of
            // active quads as well
            for (unsigned int level=0; level<number_cache.n_levels; ++level)
              {
                // count quads on this level
                number_cache.n_active_quads_level[level] = 0;

                active_quad_iterator quad = triangulation.begin_active_quad (level),
                                     endc = triangulation.end_quad ();
                for (; (quad!=endc) && (quad->level() == static_cast<signed int>(level)); ++quad)
                  ++number_cache.n_active_quads_level[level];

                // update total number of quads
                number_cache.n_active_quads += number_cache.n_active_quads_level[level];
              }
          }
        else
          {
            // for dim>2, there are no
            // levels for quads
            {
              quad_iterator quad = triangulation.begin_quad (),
                            endc = triangulation.end_quad();
              for (; quad!=endc; ++quad)
                ++number_cache.n_quads;
            }

            {
              active_quad_iterator quad = triangulation.begin_active_quad (),
                                   endc = triangulation.end_quad();
              for (; quad!=endc; ++quad)
                ++number_cache.n_active_quads;
            }
          }
      }

      /**
       * For a given Triangulation, update the
       * number cache for hexes. For 3d, we have
       * to deal with the fact that hexes have
       * levels, whereas for higher dimensions
       * they do not.
       *
       * The second argument indicates
       * for how many levels the
       * Triangulation has objects,
       * though the highest levels need
       * not contain active cells if they
       * have previously all been
       * coarsened away.
       *
       * At the end of the function, we call the
       * respective function to update the number
       * cache for quads, which will in turn call
       * the respective function for lines.
       */
      template <int dim, int spacedim>
      static
      void compute_number_cache (const Triangulation<dim,spacedim>       &triangulation,
                                 const unsigned int                       level_objects,
                                 internal::Triangulation::NumberCache<3> &number_cache)
      {
        // update quads, lines and n_levels
        compute_number_cache (triangulation,
                              level_objects,
                              static_cast<internal::Triangulation::NumberCache<2>&>
                              (number_cache));

        typedef
        typename Triangulation<dim,spacedim>::hex_iterator hex_iterator;
        typedef
        typename Triangulation<dim,spacedim>::active_hex_iterator active_hex_iterator;

        ///////////////////////////////////
        // update the number of hexes
        // on the different levels in
        // the cache
        number_cache.n_hexes_level.resize (number_cache.n_levels);
        number_cache.n_hexes = 0;

        number_cache.n_active_hexes_level.resize (number_cache.n_levels);
        number_cache.n_active_hexes = 0;

        // for 3d, hexes have levels so take
        // count the objects per level and
        // globally
        if (dim == 3)
          {
            for (unsigned int level=0; level<number_cache.n_levels; ++level)
              {
                // count hexes on this level
                number_cache.n_hexes_level[level] = 0;

                hex_iterator hex = triangulation.begin_hex (level),
                             endc = (level == number_cache.n_levels-1 ?
                                     hex_iterator(triangulation.end_hex()) :
                                     triangulation.begin_hex (level+1));
                for (; hex!=endc; ++hex)
                  ++number_cache.n_hexes_level[level];

                // update total number of hexes
                number_cache.n_hexes += number_cache.n_hexes_level[level];
              }

            // do the update for the number of
            // active hexes as well
            for (unsigned int level=0; level<number_cache.n_levels; ++level)
              {
                // count hexes on this level
                number_cache.n_active_hexes_level[level] = 0;

                active_hex_iterator hex = triangulation.begin_active_hex (level),
                                    endc = triangulation.end_hex ();
                for (; (hex!=endc) && (hex->level() == static_cast<signed int>(level)); ++hex)
                  ++number_cache.n_active_hexes_level[level];

                // update total number of hexes
                number_cache.n_active_hexes += number_cache.n_active_hexes_level[level];
              }
          }
        else
          {
            // for dim>3, there are no
            // levels for hexs
            {
              hex_iterator hex  = triangulation.begin_hex (),
                           endc = triangulation.end_hex();
              for (; hex!=endc; ++hex)
                ++number_cache.n_hexes;
            }

            {
              active_hex_iterator hex  = triangulation.begin_active_hex (),
                                  endc = triangulation.end_hex();
              for (; hex!=endc; ++hex)
                ++number_cache.n_active_hexes;
            }
          }
      }


      /**
       * Create a triangulation from
       * given data. This function does
       * this work for 1-dimensional
       * triangulations independently
       * of the actual space dimension.
       */
      template <int spacedim>
      static
      void
      create_triangulation (const std::vector<Point<spacedim> > &v,
                            const std::vector<CellData<1> >     &cells,
                            const SubCellData                   &/*subcelldata*/,
                            Triangulation<1,spacedim>           &triangulation)
      {
        AssertThrow (v.size() > 0, ExcMessage ("No vertices given"));
        AssertThrow (cells.size() > 0, ExcMessage ("No cells given"));

        // note: since no boundary
        // information can be given in one
        // dimension, the @p{subcelldata}
        // field is ignored. (only used for
        // error checking, which is a good
        // idea in any case)
        const unsigned int dim=1;

        // copy vertices
        triangulation.vertices = v;
        triangulation.vertices_used = std::vector<bool> (v.size(), true);

        // store the indices of the lines
        // which are adjacent to a given
        // vertex
        std::vector<std::vector<int> > lines_at_vertex (v.size());

        // reserve enough space
        triangulation.levels.push_back (new internal::Triangulation::TriaLevel<dim>);
        triangulation.levels[0]->reserve_space (cells.size(), dim, spacedim);
        triangulation.levels[0]->cells.reserve_space (0,cells.size());

        // make up cells
        typename Triangulation<dim,spacedim>::raw_line_iterator
        next_free_line = triangulation.begin_raw_line ();
        for (unsigned int cell=0; cell<cells.size(); ++cell)
          {
            while (next_free_line->used())
              ++next_free_line;

            next_free_line->set (internal::Triangulation
                                 ::TriaObject<1> (cells[cell].vertices[0],
                                                  cells[cell].vertices[1]));
            next_free_line->set_used_flag ();
            next_free_line->set_material_id (cells[cell].material_id);
            next_free_line->set_manifold_id (cells[cell].manifold_id);
            next_free_line->clear_user_data ();
            next_free_line->set_subdomain_id (0);

            // note that this cell is
            // adjacent to these vertices
            lines_at_vertex[cells[cell].vertices[0]].push_back (cell);
            lines_at_vertex[cells[cell].vertices[1]].push_back (cell);
          }


        // some security tests
        {
          unsigned int boundary_nodes = 0;
          for (unsigned int i=0; i<lines_at_vertex.size(); ++i)
            switch (lines_at_vertex[i].size())
              {
              case 1:
                // this vertex has only
                // one adjacent line
                ++boundary_nodes;
                break;
              case 2:
                break;
              default:
                // a node must have one
                // or two adjacent
                // lines
                AssertThrow (false, ExcInternalError());
              }

          // assert there are no more
          // than two boundary
          // nodes. note that if the
          // space dimension is
          // bigger than 1, then we
          // can have fewer than 2
          // nodes (for example a
          // ring of cells -- no end
          // points at all)
          AssertThrow (((spacedim == 1) && (boundary_nodes == 2))
                       ||
                       (spacedim > 1),
                       ExcMessage("The Triangulation has too many end points"));
        }



        // update neighborship info
        typename Triangulation<dim,spacedim>::active_line_iterator
        line = triangulation.begin_active_line ();
        // for all lines
        for (; line!=triangulation.end(); ++line)
          // for each of the two vertices
          for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
            // if first cell adjacent to
            // this vertex is the present
            // one, then the neighbor is
            // the second adjacent cell and
            // vice versa
            if (lines_at_vertex[line->vertex_index(vertex)][0] == line->index())
              if (lines_at_vertex[line->vertex_index(vertex)].size() == 2)
                {
                  const typename Triangulation<dim,spacedim>::cell_iterator
                  neighbor (&triangulation,
                            0,              // level
                            lines_at_vertex[line->vertex_index(vertex)][1]);
                  line->set_neighbor (vertex, neighbor);
                }
              else
                // no second adjacent cell
                // entered -> cell at
                // boundary
                line->set_neighbor (vertex, triangulation.end());
            else
              // present line is not first
              // adjacent one -> first
              // adjacent one is neighbor
              {
                const typename Triangulation<dim,spacedim>::cell_iterator
                neighbor (&triangulation,
                          0,              // level
                          lines_at_vertex[line->vertex_index(vertex)][0]);
                line->set_neighbor (vertex, neighbor);
              }

        // finally set the
        // vertex_to_boundary_id_map_1d
        // and vertex_to_manifold_id_map_1d
        // maps
        triangulation.vertex_to_boundary_id_map_1d->clear();
        triangulation.vertex_to_manifold_id_map_1d->clear();
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
             cell = triangulation.begin_active();
             cell != triangulation.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            {
              (*triangulation
               .vertex_to_manifold_id_map_1d)[cell->face(f)->vertex_index()]
                = numbers::invalid_manifold_id;

              if (cell->at_boundary(f))
                (*triangulation
                 .vertex_to_boundary_id_map_1d)[cell->face(f)->vertex_index()]
                  = f;
            }
      }


      /**
       * Create a triangulation from
       * given data. This function does
       * this work for 2-dimensional
       * triangulations independently
       * of the actual space dimension.
       */
      template <int spacedim>
      static
      void
      create_triangulation (const std::vector<Point<spacedim> > &v,
                            const std::vector<CellData<2> >     &cells,
                            const SubCellData                   &subcelldata,
                            Triangulation<2,spacedim>           &triangulation)
      {
        AssertThrow (v.size() > 0, ExcMessage ("No vertices given"));
        AssertThrow (cells.size() > 0, ExcMessage ("No cells given"));

        const unsigned int dim=2;

        // copy vertices
        triangulation.vertices = v;
        triangulation.vertices_used = std::vector<bool> (v.size(), true);

        // make up a list of the needed
        // lines each line is a pair of
        // vertices. The list is kept
        // sorted and it is guaranteed that
        // each line is inserted only once.
        // While the key of such an entry
        // is the pair of vertices, the
        // thing it points to is an
        // iterator pointing to the line
        // object itself. In the first run,
        // these iterators are all invalid
        // ones, but they are filled
        // afterwards
        std::map<std::pair<int,int>,
            typename Triangulation<dim,spacedim>::line_iterator> needed_lines;
        for (unsigned int cell=0; cell<cells.size(); ++cell)
          {
            for (unsigned int vertex=0; vertex<4; ++vertex)
              AssertThrow (cells[cell].vertices[vertex] < triangulation.vertices.size(),
                           ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex],
                                                  triangulation.vertices.size()));

            for (unsigned int line=0; line<GeometryInfo<dim>::faces_per_cell; ++line)
              {
                // given a line vertex number
                // (0,1) on a specific line we
                // get the cell vertex number
                // (0-4) through the
                // line_to_cell_vertices
                // function
                std::pair<int,int> line_vertices(
                  cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
                  cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);

                // assert that the line was
                // not already inserted in
                // reverse order. This
                // happens in spite of the
                // vertex rotation above,
                // if the sense of the cell
                // was incorrect.
                //
                // Here is what usually
                // happened when this
                // exception is thrown:
                // consider these two cells
                // and the vertices
                //  3---4---5
                //  |   |   |
                //  0---1---2
                // If in the input vector
                // the two cells are given
                // with vertices <0 1 4 3>
                // and <4 1 2 5>, in the
                // first cell the middle
                // line would have
                // direction 1->4, while in
                // the second it would be
                // 4->1.  This will cause
                // the exception.
                AssertThrow (needed_lines.find(std::make_pair(line_vertices.second,
                                                              line_vertices.first))
                             ==
                             needed_lines.end(),
                             ExcGridHasInvalidCell(cell));

                // insert line, with
                // invalid iterator if line
                // already exists, then
                // nothing bad happens here
                needed_lines[line_vertices] = triangulation.end_line();
              }
          }


        // check that every vertex has at
        // least two adjacent lines
        {
          std::vector<unsigned short int> vertex_touch_count (v.size(), 0);
          typename std::map<std::pair<int,int>,
                   typename Triangulation<dim,spacedim>::line_iterator>::iterator i;
          for (i=needed_lines.begin(); i!=needed_lines.end(); i++)
            {
              // touch the vertices of
              // this line
              ++vertex_touch_count[i->first.first];
              ++vertex_touch_count[i->first.second];
            }

          // assert minimum touch count
          // is at least two. if not so,
          // then clean triangulation and
          // exit with an exception
          AssertThrow (* (std::min_element(vertex_touch_count.begin(),
                                           vertex_touch_count.end())) >= 2,
                       ExcGridHasInvalidVertices());
        }

        // reserve enough space
        triangulation.levels.push_back (new internal::Triangulation::TriaLevel<dim>);
        triangulation.faces = new internal::Triangulation::TriaFaces<dim>;
        triangulation.levels[0]->reserve_space (cells.size(), dim, spacedim);
        triangulation.faces->lines.reserve_space (0,needed_lines.size());
        triangulation.levels[0]->cells.reserve_space (0,cells.size());

        // make up lines
        {
          typename Triangulation<dim,spacedim>::raw_line_iterator
          line = triangulation.begin_raw_line();
          typename std::map<std::pair<int,int>,
                   typename Triangulation<dim,spacedim>::line_iterator>::iterator i;
          for (i = needed_lines.begin();
               line!=triangulation.end_line(); ++line, ++i)
            {
              line->set (internal::Triangulation::TriaObject<1>(i->first.first,
                                                                i->first.second));
              line->set_used_flag ();
              line->clear_user_flag ();
              line->clear_user_data ();
              i->second = line;
            }
        }


        // store for each line index
        // the adjacent cells
        std::map<int,std::vector<typename Triangulation<dim,spacedim>::cell_iterator> >
        adjacent_cells;

        // finally make up cells
        {
          typename Triangulation<dim,spacedim>::raw_cell_iterator
          cell = triangulation.begin_raw_quad();
          for (unsigned int c=0; c<cells.size(); ++c, ++cell)
            {
              typename Triangulation<dim,spacedim>::line_iterator
              lines[GeometryInfo<dim>::lines_per_cell];
              for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                lines[line]=needed_lines[std::make_pair(
                                           cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
                                           cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)])];

              cell->set (internal::Triangulation::TriaObject<2> (lines[0]->index(),
                                                                 lines[1]->index(),
                                                                 lines[2]->index(),
                                                                 lines[3]->index()));

              cell->set_used_flag ();
              cell->set_material_id (cells[c].material_id);
              cell->set_manifold_id (cells[c].manifold_id);
              cell->clear_user_data ();
              cell->set_subdomain_id (0);

              // note that this cell is
              // adjacent to the four
              // lines
              for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                adjacent_cells[lines[line]->index()].push_back (cell);
            }
        }


        for (typename Triangulation<dim,spacedim>::line_iterator
             line=triangulation.begin_line();
             line!=triangulation.end_line(); ++line)
          {
            const unsigned int n_adj_cells = adjacent_cells[line->index()].size();
            // assert that every line has
            // one or two adjacent cells
            AssertThrow ((n_adj_cells >= 1) &&
                         (n_adj_cells <= 2),
                         ExcInternalError());

            // if only one cell: line is at
            // boundary -> give it the
            // boundary indicator zero by
            // default
            if (n_adj_cells == 1)
              line->set_boundary_indicator (0);
            else
              // interior line -> numbers::internal_face_boundary_id
              line->set_boundary_indicator (numbers::internal_face_boundary_id);
            line->set_manifold_id(numbers::flat_manifold_id);
          }

        // set boundary indicators where
        // given
        std::vector<CellData<1> >::const_iterator boundary_line
          = subcelldata.boundary_lines.begin();
        std::vector<CellData<1> >::const_iterator end_boundary_line
          = subcelldata.boundary_lines.end();
        for (; boundary_line!=end_boundary_line; ++boundary_line)
          {
            typename Triangulation<dim,spacedim>::line_iterator line;
            std::pair<int,int> line_vertices(std::make_pair(boundary_line->vertices[0],
                                                            boundary_line->vertices[1]));
            if (needed_lines.find(line_vertices) != needed_lines.end())
              // line found in this
              // direction
              line = needed_lines[line_vertices];
            else
              {
                // look whether it exists
                // in reverse direction
                std::swap (line_vertices.first, line_vertices.second);
                if (needed_lines.find(line_vertices) != needed_lines.end())
                  line = needed_lines[line_vertices];
                else
                  // line does not exist
                  AssertThrow (false, ExcLineInexistant(line_vertices.first,
                                                        line_vertices.second));
              }

            // assert that we only set
            // boundary info once
            AssertThrow (! (line->boundary_indicator() != 0 &&
                            line->boundary_indicator() != numbers::internal_face_boundary_id),
                         ExcMultiplySetLineInfoOfLine(line_vertices.first,
                                                      line_vertices.second));

            // Assert that only exterior lines
            // are given a boundary indicator
            AssertThrow (! (line->boundary_indicator() == numbers::internal_face_boundary_id),
                         ExcInteriorLineCantBeBoundary());

            line->set_boundary_indicator (boundary_line->boundary_id);
            line->set_manifold_id (boundary_line->manifold_id);
          }


        // finally update neighborship info
        for (typename Triangulation<dim,spacedim>::cell_iterator
             cell=triangulation.begin(); cell!=triangulation.end(); ++cell)
          for (unsigned int side=0; side<4; ++side)
            if (adjacent_cells[cell->line(side)->index()][0] == cell)
              // first adjacent cell is
              // this one
              {
                if (adjacent_cells[cell->line(side)->index()].size() == 2)
                  // there is another
                  // adjacent cell
                  cell->set_neighbor (side,
                                      adjacent_cells[cell->line(side)->index()][1]);
              }
        // first adjacent cell is not this
        // one, -> it must be the neighbor
        // we are looking for
            else
              cell->set_neighbor (side,
                                  adjacent_cells[cell->line(side)->index()][0]);
      }


      /**
       * Invent an object which compares two internal::Triangulation::TriaObject<2>
       * against each other. This comparison is needed in order to establish a map
       * of TriaObject<2> to iterators in the Triangulation<3,3>::create_triangulation
       * function.
       *
       * Since this comparison is not canonical, we do not include it into the
       * general internal::Triangulation::TriaObject<2> class.
       */
      struct QuadComparator
      {
        inline bool operator () (const internal::Triangulation::TriaObject<2> &q1,
                                 const internal::Triangulation::TriaObject<2> &q2) const
        {
          // here is room to
          // optimize the repeated
          // equality test of the
          // previous lines; the
          // compiler will probably
          // take care of most of
          // it anyway
          if ((q1.face(0) < q2.face(0))          ||
              ((q1.face(0) == q2.face(0)) &&
               (q1.face(1) <  q2.face(1)))       ||
              ((q1.face(0) == q2.face(0)) &&
               (q1.face(1) == q2.face(1)) &&
               (q1.face(2) <  q2.face(2)))       ||
              ((q1.face(0) == q2.face(0)) &&
               (q1.face(1) == q2.face(1)) &&
               (q1.face(2) == q2.face(2)) &&
               (q1.face(3) <  q2.face(3))))
            return true;
          else
            return false;
        }
      };



      template <int spacedim>
      static
      void
      create_triangulation (const std::vector<Point<spacedim> > &v,
                            const std::vector<CellData<3> >     &cells,
                            const SubCellData                   &subcelldata,
                            Triangulation<3,spacedim>           &triangulation)
      {
        AssertThrow (v.size() > 0, ExcMessage ("No vertices given"));
        AssertThrow (cells.size() > 0, ExcMessage ("No cells given"));

        const unsigned int dim=3;

        // copy vertices
        triangulation.vertices = v;
        triangulation.vertices_used = std::vector<bool> (v.size(), true);

        // check that all cells have
        // positive volume. if not call the
        // invert_all_cells_of_negative_grid
        // and reorder_cells function of
        // GridReordering before creating
        // the triangulation
        for (unsigned int cell_no=0; cell_no<cells.size(); ++cell_no)
          AssertThrow (dealii::GridTools::cell_measure(triangulation.vertices,
                                                       cells[cell_no].vertices) >= 0,
                       ExcGridHasInvalidCell(cell_no));

        ///////////////////////////////////////
        // first set up some collections of data
        //
        // make up a list of the needed
        // lines
        //
        // each line is a pair of
        // vertices. The list is kept
        // sorted and it is guaranteed that
        // each line is inserted only once.
        // While the key of such an entry
        // is the pair of vertices, the
        // thing it points to is an
        // iterator pointing to the line
        // object itself. In the first run,
        // these iterators are all invalid
        // ones, but they are filled
        // afterwards same applies for the
        // quads
        typename std::map<std::pair<int,int>,
                 typename Triangulation<dim,spacedim>::line_iterator> needed_lines;
        for (unsigned int cell=0; cell<cells.size(); ++cell)
          {
            // check whether vertex indices
            // are valid ones
            for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
              AssertThrow (cells[cell].vertices[vertex] < triangulation.vertices.size(),
                           ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex],
                                                  triangulation.vertices.size()));

            for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
              {
                // given a line vertex number
                // (0,1) on a specific line we
                // get the cell vertex number
                // (0-7) through the
                // line_to_cell_vertices
                // function
                std::pair<int,int> line_vertices(
                  cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
                  cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);

                // if that line was already inserted
                // in reverse order do nothing, else
                // insert the line
                if ( (needed_lines.find(std::make_pair(line_vertices.second,
                                                       line_vertices.first))
                      ==
                      needed_lines.end()))
                  {
                    // insert line, with
                    // invalid iterator. if line
                    // already exists, then
                    // nothing bad happens here
                    needed_lines[line_vertices] = triangulation.end_line();
                  }
              }
          }


        /////////////////////////////////
        // now for some sanity-checks:
        //
        // check that every vertex has at
        // least tree adjacent lines
        {
          std::vector<unsigned short int> vertex_touch_count (v.size(), 0);
          typename std::map<std::pair<int,int>,
                   typename Triangulation<dim,spacedim>::line_iterator>::iterator i;
          for (i=needed_lines.begin(); i!=needed_lines.end(); i++)
            {
              // touch the vertices of
              // this line
              ++vertex_touch_count[i->first.first];
              ++vertex_touch_count[i->first.second];
            }

          // assert minimum touch count
          // is at least three. if not so,
          // then clean triangulation and
          // exit with an exception
          AssertThrow (* (std::min_element(vertex_touch_count.begin(),
                                           vertex_touch_count.end())) >= 3,
                       ExcGridHasInvalidVertices());
        }


        ///////////////////////////////////
        // actually set up data structures
        // for the lines
        // reserve enough space
        triangulation.levels.push_back (new internal::Triangulation::TriaLevel<dim>);
        triangulation.faces = new internal::Triangulation::TriaFaces<dim>;
        triangulation.levels[0]->reserve_space (cells.size(), dim, spacedim);
        triangulation.faces->lines.reserve_space (0,needed_lines.size());

        // make up lines
        {
          typename Triangulation<dim,spacedim>::raw_line_iterator
          line = triangulation.begin_raw_line();
          typename std::map<std::pair<int,int>,
                   typename Triangulation<dim,spacedim>::line_iterator>::iterator i;
          for (i = needed_lines.begin(); line!=triangulation.end_line(); ++line, ++i)
            {
              line->set (internal::Triangulation::TriaObject<1>(i->first.first,
                                                                i->first.second));
              line->set_used_flag ();
              line->clear_user_flag ();
              line->clear_user_data ();

              // now set the iterator for
              // this line
              i->second = line;
            }
        }


        ///////////////////////////////////////////
        // make up the quads of this triangulation
        //
        // same thing: the iterators are
        // set to the invalid value at
        // first, we only collect the data
        // now

        // the bool array stores, whether the lines
        // are in the standard orientation or not

        // note that QuadComparator is a
        // class declared and defined in
        // this file
        std::map<internal::Triangulation::TriaObject<2>,
            std::pair<typename Triangulation<dim,spacedim>::quad_iterator,
            std_cxx11::array<bool,GeometryInfo<dim>::lines_per_face> >,
            QuadComparator>
            needed_quads;
        for (unsigned int cell=0; cell<cells.size(); ++cell)
          {
            // the faces are quads which
            // consist of four numbers
            // denoting the index of the
            // four lines bounding the
            // quad. we can get this index
            // by asking @p{needed_lines}
            // for an iterator to this
            // line, dereferencing it and
            // thus return an iterator into
            // the @p{lines} array of the
            // triangulation, which is
            // already set up. we can then
            // ask this iterator for its
            // index within the present
            // level (the level is zero, of
            // course)
            //
            // to make things easier, we
            // don't create the lines
            // (pairs of their vertex
            // indices) in place, but
            // before they are really
            // needed.
            std::pair<int,int> line_list[GeometryInfo<dim>::lines_per_cell],
                inverse_line_list[GeometryInfo<dim>::lines_per_cell];
            unsigned int face_line_list[GeometryInfo<dim>::lines_per_face];
            std_cxx11::array<bool,GeometryInfo<dim>::lines_per_face> orientation;

            for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
              {
                line_list[line]=std::pair<int,int> (
                                  cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
                                  cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);
                inverse_line_list[line]=std::pair<int,int> (
                                          cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)],
                                          cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)]);
              }

            for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
              {
                // set up a list of the lines to be
                // used for this face. check the
                // direction for each line
                //
                // given a face line number (0-3) on
                // a specific face we get the cell
                // line number (0-11) through the
                // face_to_cell_lines function
                for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
                  if (needed_lines.find (inverse_line_list[GeometryInfo<dim>::
                                                           face_to_cell_lines(face,l)]) == needed_lines.end())
                    {
                      face_line_list[l]=needed_lines[line_list[GeometryInfo<dim>::
                                                               face_to_cell_lines(face,l)]]->index();
                      orientation[l]=true;
                    }
                  else
                    {
                      face_line_list[l]=needed_lines[inverse_line_list[GeometryInfo<dim>::
                                                                       face_to_cell_lines(face,l)]]->index();
                      orientation[l]=false;
                    }


                internal::Triangulation::TriaObject<2>
                quad(face_line_list[0],
                     face_line_list[1],
                     face_line_list[2],
                     face_line_list[3]);

                // insert quad, with
                // invalid iterator
                //
                // if quad already exists,
                // then nothing bad happens
                // here, as this will then
                // simply become an
                // interior face of the
                // triangulation. however,
                // we will run into major
                // trouble if the face was
                // already inserted in the
                // opposite
                // direction. there are
                // really only two
                // orientations for a face
                // to be in, since the edge
                // directions are already
                // set. thus, vertex 0 is
                // the one from which two
                // edges originate, and
                // vertex 3 is the one to
                // which they converge. we
                // are then left with
                // orientations 0-1-2-3 and
                // 2-3-0-1 for the order of
                // lines. the
                // corresponding quad can
                // be easily constructed by
                // exchanging lines. we do
                // so here, just to check
                // that that flipped quad
                // isn't already in the
                // triangulation. if it is,
                // then don't insert the
                // new one and instead
                // later set the
                // face_orientation flag
                const internal::Triangulation::TriaObject<2>
                test_quad_1(quad.face(2), quad.face(3),
                            quad.face(0), quad.face(1)),//face_orientation=false, face_flip=false, face_rotation=false
                                      test_quad_2(quad.face(0), quad.face(1),
                                                  quad.face(3), quad.face(2)),//face_orientation=false, face_flip=false, face_rotation=true
                                      test_quad_3(quad.face(3), quad.face(2),
                                                  quad.face(1), quad.face(0)),//face_orientation=false, face_flip=true,  face_rotation=false
                                      test_quad_4(quad.face(1), quad.face(0),
                                                  quad.face(2), quad.face(3)),//face_orientation=false, face_flip=true,  face_rotation=true
                                      test_quad_5(quad.face(2), quad.face(3),
                                                  quad.face(1), quad.face(0)),//face_orientation=true,  face_flip=false, face_rotation=true
                                      test_quad_6(quad.face(1), quad.face(0),
                                                  quad.face(3), quad.face(2)),//face_orientation=true,  face_flip=true,  face_rotation=false
                                      test_quad_7(quad.face(3), quad.face(2),
                                                  quad.face(0), quad.face(1));//face_orientation=true,  face_flip=true,  face_rotation=true
                if (needed_quads.find (test_quad_1) == needed_quads.end() &&
                    needed_quads.find (test_quad_2) == needed_quads.end() &&
                    needed_quads.find (test_quad_3) == needed_quads.end() &&
                    needed_quads.find (test_quad_4) == needed_quads.end() &&
                    needed_quads.find (test_quad_5) == needed_quads.end() &&
                    needed_quads.find (test_quad_6) == needed_quads.end() &&
                    needed_quads.find (test_quad_7) == needed_quads.end())
                  needed_quads[quad] = std::make_pair(triangulation.end_quad(),orientation);
              }
          }


        /////////////////////////////////
        // enter the resulting quads into
        // the arrays of the Triangulation
        //
        // first reserve enough space
        triangulation.faces->quads.reserve_space (0,needed_quads.size());

        {
          typename Triangulation<dim,spacedim>::raw_quad_iterator
          quad = triangulation.begin_raw_quad();
          typename std::map<internal::Triangulation::TriaObject<2>,
                   std::pair<typename Triangulation<dim,spacedim>::quad_iterator,
                   std_cxx11::array<bool,GeometryInfo<dim>::lines_per_face> >,
                   QuadComparator>
                   ::iterator q;
          for (q = needed_quads.begin(); quad!=triangulation.end_quad(); ++quad, ++q)
            {
              quad->set (q->first);
              quad->set_used_flag ();
              quad->clear_user_flag ();
              quad->clear_user_data ();
              // set the line orientation
              quad->set_line_orientation(0,q->second.second[0]);
              quad->set_line_orientation(1,q->second.second[1]);
              quad->set_line_orientation(2,q->second.second[2]);
              quad->set_line_orientation(3,q->second.second[3]);


              // now set the iterator for
              // this quad
              q->second.first = quad;
            }
        }

        /////////////////////////////////
        // finally create the cells
        triangulation.levels[0]->cells.reserve_space (cells.size());

        // store for each quad index the
        // adjacent cells
        std::map<int,std::vector<typename Triangulation<dim,spacedim>::cell_iterator> >
        adjacent_cells;

        // finally make up cells
        {
          typename Triangulation<dim,spacedim>::raw_cell_iterator
          cell = triangulation.begin_raw_hex();
          for (unsigned int c=0; c<cells.size(); ++c, ++cell)
            {
              // first find for each of
              // the cells the quad
              // iterator of the
              // respective faces.
              //
              // to this end, set up the
              // lines of this cell and
              // find the quads that are
              // bounded by these lines;
              // these are then the faces
              // of the present cell
              std::pair<int,int> line_list[GeometryInfo<dim>::lines_per_cell],
                  inverse_line_list[GeometryInfo<dim>::lines_per_cell];
              unsigned int face_line_list[4];
              for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                {
                  line_list[line]=std::make_pair(
                                    cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
                                    cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);
                  inverse_line_list[line]=std::pair<int,int> (
                                            cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)],
                                            cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)]);
                }

              // get the iterators
              // corresponding to the
              // faces. also store
              // whether they are
              // reversed or not
              typename Triangulation<dim,spacedim>::quad_iterator
              face_iterator[GeometryInfo<dim>::faces_per_cell];
              bool face_orientation[GeometryInfo<dim>::faces_per_cell];
              bool face_flip[GeometryInfo<dim>::faces_per_cell];
              bool face_rotation[GeometryInfo<dim>::faces_per_cell];
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                {
                  for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
                    if (needed_lines.find (inverse_line_list[GeometryInfo<dim>::
                                                             face_to_cell_lines(face,l)]) == needed_lines.end())
                      face_line_list[l]=needed_lines[line_list[GeometryInfo<dim>::
                                                               face_to_cell_lines(face,l)]]->index();
                    else
                      face_line_list[l]=needed_lines[inverse_line_list[GeometryInfo<dim>::
                                                                       face_to_cell_lines(face,l)]]->index();

                  internal::Triangulation::TriaObject<2>
                  quad(face_line_list[0],
                       face_line_list[1],
                       face_line_list[2],
                       face_line_list[3]);

                  if (needed_quads.find (quad) != needed_quads.end())
                    {
                      // face is in standard
                      // orientation (and not
                      // flipped or rotated). this
                      // must be true for at least
                      // one of the two cells
                      // containing this face
                      // (i.e. for the cell which
                      // originally inserted the
                      // face)
                      face_iterator[face] = needed_quads[quad].first;
                      face_orientation[face] = true;
                      face_flip[face]=false;
                      face_rotation[face]=false;
                    }
                  else
                    {
                      // face must be available in
                      // reverse order
                      // then. construct all
                      // possibilities and check
                      // them one after the other
                      const internal::Triangulation::TriaObject<2>
                      test_quad_1(quad.face(2), quad.face(3),
                                  quad.face(0), quad.face(1)),//face_orientation=false, face_flip=false, face_rotation=false
                                            test_quad_2(quad.face(0), quad.face(1),
                                                        quad.face(3), quad.face(2)),//face_orientation=false, face_flip=false, face_rotation=true
                                            test_quad_3(quad.face(3), quad.face(2),
                                                        quad.face(1), quad.face(0)),//face_orientation=false, face_flip=true,  face_rotation=false
                                            test_quad_4(quad.face(1), quad.face(0),
                                                        quad.face(2), quad.face(3)),//face_orientation=false, face_flip=true,  face_rotation=true
                                            test_quad_5(quad.face(2), quad.face(3),
                                                        quad.face(1), quad.face(0)),//face_orientation=true,  face_flip=false, face_rotation=true
                                            test_quad_6(quad.face(1), quad.face(0),
                                                        quad.face(3), quad.face(2)),//face_orientation=true,  face_flip=true,  face_rotation=false
                                            test_quad_7(quad.face(3), quad.face(2),
                                                        quad.face(0), quad.face(1));//face_orientation=true,  face_flip=true,  face_rotation=true
                      if (needed_quads.find (test_quad_1) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_1].first;
                          face_orientation[face] = false;
                          face_flip[face]=false;
                          face_rotation[face]=false;
                        }
                      else if (needed_quads.find (test_quad_2) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_2].first;
                          face_orientation[face] = false;
                          face_flip[face]=false;
                          face_rotation[face]=true;
                        }
                      else if (needed_quads.find (test_quad_3) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_3].first;
                          face_orientation[face] = false;
                          face_flip[face]=true;
                          face_rotation[face]=false;
                        }
                      else if (needed_quads.find (test_quad_4) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_4].first;
                          face_orientation[face] = false;
                          face_flip[face]=true;
                          face_rotation[face]=true;
                        }
                      else if (needed_quads.find (test_quad_5) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_5].first;
                          face_orientation[face] = true;
                          face_flip[face]=false;
                          face_rotation[face]=true;
                        }
                      else if (needed_quads.find (test_quad_6) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_6].first;
                          face_orientation[face] = true;
                          face_flip[face]=true;
                          face_rotation[face]=false;
                        }
                      else if (needed_quads.find (test_quad_7) != needed_quads.end())
                        {
                          face_iterator[face] = needed_quads[test_quad_7].first;
                          face_orientation[face] = true;
                          face_flip[face]=true;
                          face_rotation[face]=true;
                        }

                      else
                        // we didn't find the
                        // face in any direction,
                        // so something went
                        // wrong above
                        Assert(false,ExcInternalError());

                    }
                }// for all faces

              // make the cell out of
              // these iterators
              cell->set (internal::Triangulation
                         ::TriaObject<3> (face_iterator[0]->index(),
                                          face_iterator[1]->index(),
                                          face_iterator[2]->index(),
                                          face_iterator[3]->index(),
                                          face_iterator[4]->index(),
                                          face_iterator[5]->index()));

              cell->set_used_flag ();
              cell->set_material_id (cells[c].material_id);
              cell->set_manifold_id (cells[c].manifold_id);
              cell->clear_user_flag ();
              cell->clear_user_data ();
              cell->set_subdomain_id (0);

              // set orientation flag for
              // each of the faces
              for (unsigned int quad=0; quad<GeometryInfo<dim>::faces_per_cell; ++quad)
                {
                  cell->set_face_orientation (quad, face_orientation[quad]);
                  cell->set_face_flip (quad, face_flip[quad]);
                  cell->set_face_rotation (quad, face_rotation[quad]);
                }


              // note that this cell is
              // adjacent to the six
              // quads
              for (unsigned int quad=0; quad<GeometryInfo<dim>::faces_per_cell; ++quad)
                adjacent_cells[face_iterator[quad]->index()].push_back (cell);

#ifdef DEBUG
              // make some checks on the
              // lines and their
              // ordering

              // first map all cell lines
              // to the two face lines
              // which should
              // coincide. all face lines
              // are included with a cell
              // line number (0-11)
              // key. At the end all keys
              // will be included twice
              // (for each of the two
              // coinciding lines once)
              std::multimap<unsigned int, std::pair<unsigned int, unsigned int> >
              cell_to_face_lines;
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_face; ++line)
                  cell_to_face_lines.insert(
                    std::pair<unsigned int, std::pair<unsigned int, unsigned int> > (
                      GeometryInfo<dim>::face_to_cell_lines(face,line),
                      std::pair<unsigned int, unsigned int> (face,line)));
              std::multimap<unsigned int, std::pair<unsigned int, unsigned int> >::const_iterator
              map_iter=cell_to_face_lines.begin();

              for (; map_iter!=cell_to_face_lines.end(); ++map_iter)
                {
                  const unsigned int cell_line=map_iter->first;
                  const unsigned int face1=map_iter->second.first;
                  const unsigned int line1=map_iter->second.second;
                  ++map_iter;
                  Assert(map_iter!=cell_to_face_lines.end(), ExcInternalErrorOnCell(c));
                  Assert(map_iter->first==cell_line, ExcInternalErrorOnCell(c));
                  const unsigned int face2=map_iter->second.first;
                  const unsigned int line2=map_iter->second.second;

                  // check that the pair
                  // of lines really
                  // coincide. Take care
                  // about the face
                  // orientation;
                  Assert (face_iterator[face1]->line(GeometryInfo<dim>::standard_to_real_face_line(
                                                       line1,
                                                       face_orientation[face1],
                                                       face_flip[face1],
                                                       face_rotation[face1])) ==
                          face_iterator[face2]->line(GeometryInfo<dim>::standard_to_real_face_line(
                                                       line2,
                                                       face_orientation[face2],
                                                       face_flip[face2],
                                                       face_rotation[face2])),
                          ExcInternalErrorOnCell(c));
                }
#endif
            }
        }


        /////////////////////////////////////////
        // find those quads which are at the
        // boundary and mark them appropriately
        for (typename Triangulation<dim,spacedim>::quad_iterator
             quad=triangulation.begin_quad(); quad!=triangulation.end_quad(); ++quad)
          {
            const unsigned int n_adj_cells = adjacent_cells[quad->index()].size();
            // assert that every quad has
            // one or two adjacent cells
            AssertThrow ((n_adj_cells >= 1) &&
                         (n_adj_cells <= 2),
                         ExcInternalError());

            // if only one cell: quad is at
            // boundary -> give it the
            // boundary indicator zero by
            // default
            if (n_adj_cells == 1)
              quad->set_boundary_indicator (0);
            else
              // interior quad -> numbers::internal_face_boundary_id
              quad->set_boundary_indicator (numbers::internal_face_boundary_id);
            // Manifold ids are set
            // independently of where
            // they are
            quad->set_manifold_id(numbers::flat_manifold_id);
          }

        /////////////////////////////////////////
        // next find those lines which are at
        // the boundary and mark all others as
        // interior ones
        //
        // for this: first mark all lines as interior. use this loop
        // to also set all manifold ids of all lines
        for (typename Triangulation<dim,spacedim>::line_iterator
             line=triangulation.begin_line(); line!=triangulation.end_line(); ++line)
          {
            line->set_boundary_indicator (numbers::internal_face_boundary_id);
            line->set_manifold_id(numbers::flat_manifold_id);
          }

        // next reset all lines bounding
        // boundary quads as on the
        // boundary also. note that since
        // we are in 3d, there are cases
        // where one or more lines of a
        // quad that is not on the
        // boundary, are actually boundary
        // lines. they will not be marked
        // when visiting this
        // face. however, since we do not
        // support dim-2 dimensional
        // boundaries (i.e. internal lines
        // constituting boundaries), every
        // such line is also part of a face
        // that is actually on the
        // boundary, so sooner or later we
        // get to mark that line for being
        // on the boundary
        for (typename Triangulation<dim,spacedim>::quad_iterator
             quad=triangulation.begin_quad(); quad!=triangulation.end_quad(); ++quad)
          if (quad->at_boundary())
            for (unsigned int l=0; l<4; ++l)
              quad->line(l)->set_boundary_indicator (0);

        ///////////////////////////////////////
        // now set boundary indicators
        // where given
        //
        // first do so for lines
        std::vector<CellData<1> >::const_iterator boundary_line
          = subcelldata.boundary_lines.begin();
        std::vector<CellData<1> >::const_iterator end_boundary_line
          = subcelldata.boundary_lines.end();
        for (; boundary_line!=end_boundary_line; ++boundary_line)
          {
            typename Triangulation<dim,spacedim>::line_iterator line;
            std::pair <int, int> line_vertices(std::make_pair(boundary_line->vertices[0],
                                                              boundary_line->vertices[1]));
            if (needed_lines.find(line_vertices) != needed_lines.end())
              // line found in this
              // direction
              line = needed_lines[line_vertices];

            else
              {
                // look whether it exists in
                // reverse direction
                std::swap (line_vertices.first, line_vertices.second);
                if (needed_lines.find(line_vertices) != needed_lines.end())
                  line = needed_lines[line_vertices];
                else
                  // line does not exist
                  AssertThrow (false, ExcLineInexistant(line_vertices.first,
                                                        line_vertices.second));
              }
            // Assert that only exterior
            // lines are given a boundary
            // indicator
            AssertThrow (line->at_boundary(),
                         ExcInteriorLineCantBeBoundary());

            // and make sure that we don't
            // attempt to reset the
            // boundary indicator to a
            // different than the
            // previously set value
            if (line->boundary_indicator() != 0)
              AssertThrow (line->boundary_indicator() == boundary_line->boundary_id,
                           ExcMessage ("Duplicate boundary lines are only allowed "
                                       "if they carry the same boundary indicator."));

            line->set_boundary_indicator (boundary_line->boundary_id);
            // Set manifold id if given
            line->set_manifold_id(boundary_line->manifold_id);
          }


        // now go on with boundary faces
        std::vector<CellData<2> >::const_iterator boundary_quad
          = subcelldata.boundary_quads.begin();
        std::vector<CellData<2> >::const_iterator end_boundary_quad
          = subcelldata.boundary_quads.end();
        for (; boundary_quad!=end_boundary_quad; ++boundary_quad)
          {
            typename Triangulation<dim,spacedim>::quad_iterator quad;
            typename Triangulation<dim,spacedim>::line_iterator line[4];

            // first find the lines that
            // are made up of the given
            // vertices, then build up a
            // quad from these lines
            // finally use the find
            // function of the map template
            // to find the quad
            for (unsigned int i=0; i<4; ++i)
              {
                std::pair<int, int> line_vertices(
                  boundary_quad->vertices[GeometryInfo<dim-1>::line_to_cell_vertices(i,0)],
                  boundary_quad->vertices[GeometryInfo<dim-1>::line_to_cell_vertices(i,1)]);

                // check whether line
                // already exists
                if (needed_lines.find(line_vertices) != needed_lines.end())
                  line[i] = needed_lines[line_vertices];
                else
                  // look whether it exists
                  // in reverse direction
                  {
                    std::swap (line_vertices.first, line_vertices.second);
                    if (needed_lines.find(line_vertices) != needed_lines.end())
                      line[i] = needed_lines[line_vertices];
                    else
                      // line does
                      // not exist
                      AssertThrow (false, ExcLineInexistant(line_vertices.first,
                                                            line_vertices.second));
                  }
              }


            // Set up 2 quads that are
            // built up from the lines for
            // reasons of comparison to
            // needed_quads.  The second
            // quad is the reversed version
            // of the first quad in order
            // find the quad regardless of
            // its orientation.  This is
            // introduced for convenience
            // and because boundary quad
            // orientation does not carry
            // any information.
            internal::Triangulation::TriaObject<2>
            quad_compare_1(line[0]->index(), line[1]->index(),
                           line[2]->index(), line[3]->index());
            internal::Triangulation::TriaObject<2>
            quad_compare_2(line[2]->index(), line[3]->index(),
                           line[0]->index(), line[1]->index());

            // try to find the quad with
            // lines situated as
            // constructed above.  if it
            // could not be found, rotate
            // the boundary lines 3 times
            // until it is found or it does
            // not exist.

            // mapping from counterclock to
            // lexicographic ordering of
            // quad lines
            static const unsigned int lex2cclock[4]= {3,1,0,2};
            // copy lines from
            // lexicographic to
            // counterclock ordering, as
            // rotation is much simpler in
            // counterclock ordering
            typename Triangulation<dim,spacedim>::line_iterator
            line_counterclock[4];
            for (unsigned int i=0; i<4; ++i)
              line_counterclock[lex2cclock[i]]=line[i];
            unsigned int n_rotations=0;
            bool not_found_quad_1;
            while ( (not_found_quad_1=(needed_quads.find(quad_compare_1) == needed_quads.end())) &&
                    (                  needed_quads.find(quad_compare_2) == needed_quads.end()) &&
                    (n_rotations<4))
              {
                // use the rotate defined
                // in <algorithms>
                rotate(line_counterclock, line_counterclock+1, line_counterclock+4);
                // update the quads with
                // rotated lines (i runs in
                // lexicographic ordering)
                for (unsigned int i=0; i<4; ++i)
                  {
                    quad_compare_1.set_face(i,       line_counterclock[lex2cclock[i]]->index());
                    quad_compare_2.set_face((i+2)%4, line_counterclock[lex2cclock[i]]->index());
                  }

                ++n_rotations;
              }

            AssertThrow (n_rotations!=4,
                         ExcQuadInexistant(line[0]->index(), line[1]->index(),
                                           line[2]->index(), line[3]->index()));

            if (not_found_quad_1)
              quad = needed_quads[quad_compare_2].first;
            else
              quad = needed_quads[quad_compare_1].first;

            // check whether this face is
            // really an exterior one
            AssertThrow (quad->at_boundary(),
                         ExcInteriorQuadCantBeBoundary());

            // and make sure that we don't
            // attempt to reset the
            // boundary indicator to a
            // different than the
            // previously set value
            if (quad->boundary_indicator() != 0)
              AssertThrow (quad->boundary_indicator() == boundary_quad->boundary_id,
                           ExcMessage ("Duplicate boundary quads are only allowed "
                                       "if they carry the same boundary indicator."));

            quad->set_boundary_indicator (boundary_quad->boundary_id);
            quad->set_manifold_id (boundary_quad->manifold_id);
          }


        /////////////////////////////////////////
        // finally update neighborship info
        for (typename Triangulation<dim,spacedim>::cell_iterator
             cell=triangulation.begin(); cell!=triangulation.end(); ++cell)
          for (unsigned int face=0; face<6; ++face)
            if (adjacent_cells[cell->quad(face)->index()][0] == cell)
              // first adjacent cell is
              // this one
              {
                if (adjacent_cells[cell->quad(face)->index()].size() == 2)
                  // there is another
                  // adjacent cell
                  cell->set_neighbor (face,
                                      adjacent_cells[cell->quad(face)->index()][1]);
              }
        // first adjacent cell is not this
        // one, -> it must be the neighbor
        // we are looking for
            else
              cell->set_neighbor (face,
                                  adjacent_cells[cell->quad(face)->index()][0]);
      }


      /**
       * Distort a 1d triangulation in
       * some random way.
       */
      template <int spacedim>
      static
      void
      distort_random (const double               factor,
                      const bool                 keep_boundary,
                      Triangulation<1,spacedim> &triangulation)
      {
        const unsigned int dim = 1;

        // if spacedim>1 we need to
        // make sure that we perturb
        // points but keep them on
        // the manifold
        Assert (spacedim == 1,
                ExcNotImplemented());

        // this function is mostly
        // equivalent to that for the
        // general dimensional case the
        // only difference being the
        // correction for split faces which
        // is not necessary in 1D

        // find the smallest length of the
        // lines adjacent to the
        // vertex. take the initial value
        // to be larger than anything that
        // might be found: the diameter of
        // the triangulation, here computed
        // by adding up the diameters of
        // the coarse grid cells.
        double almost_infinite_length = 0;
        for (typename Triangulation<dim,spacedim>::cell_iterator
             cell=triangulation.begin(0); cell!=triangulation.end(0); ++cell)
          almost_infinite_length += cell->diameter();

        std::vector<double> minimal_length (triangulation.vertices.size(),
                                            almost_infinite_length);
        std::vector<bool>   at_boundary (triangulation.vertices.size(), false);

        for (typename Triangulation<dim,spacedim>::active_line_iterator
             line=triangulation.begin_active_line();
             line != triangulation.end_line(); ++line)
          {
            minimal_length[line->vertex_index(0)]
              = std::min(line->diameter(),
                         minimal_length[line->vertex_index(0)]);
            minimal_length[line->vertex_index(1)]
              = std::min(line->diameter(),
                         minimal_length[line->vertex_index(1)]);
          }

        // also note if a vertex is at the boundary if we are asked to
        // keep boundary vertices untouched
        if (keep_boundary)
          for (typename Triangulation<dim,spacedim>::active_line_iterator
               line=triangulation.begin_active_line();
               line != triangulation.end_line(); ++line)
            for (unsigned int vertex=0; vertex<2; ++vertex)
              if (line->at_boundary(vertex) == true)
                at_boundary[line->vertex_index(vertex)] = true;


        const unsigned int n_vertices = triangulation.vertices.size();
        Point<spacedim> shift_vector;

        for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
          {
            // ignore this vertex if we
            // whall keep the boundary and
            // this vertex *is* at the
            // boundary
            if (keep_boundary && at_boundary[vertex])
              continue;

            // first compute a random shift
            // vector
            for (unsigned int d=0; d<spacedim; ++d)
              shift_vector(d) = std::rand()*2.0/RAND_MAX-1;

            shift_vector *= factor * minimal_length[vertex] /
                            std::sqrt(shift_vector.square());

            // finally move the vertex
            triangulation.vertices[vertex] += shift_vector;
          }
      }


      /**
       * Distort a triangulation in
       * some random way. This is the
       * function taken for the case
       * dim>1.
       */
      template <int dim, int spacedim>
      static
      void
      distort_random (const double                 factor,
                      const bool                   keep_boundary,
                      Triangulation<dim,spacedim> &triangulation)
      {
//TODO:[?]Implement the random distortion in Triangulation for hanging nodes as well
//    Hanging nodes need to be reset to the correct mean value
//    at the end, which is simple for 2D but difficult for 3D. Maybe take
//    a look at how we get to the original location of the point in the
//    execute_refinement function and copy the relevant lines.

        // this function is mostly
        // equivalent to that for the
        // general dimensional case the
        // only difference being the
        // correction for split faces which
        // is not necessary in 1D
        //
        // if you change something here,
        // don't forget to do so there as
        // well

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
        for (typename Triangulation<dim,spacedim>::cell_iterator
             cell=triangulation.begin(0); cell!=triangulation.end(0); ++cell)
          almost_infinite_length += cell->diameter();

        std::vector<double> minimal_length (triangulation.vertices.size(),
                                            almost_infinite_length);

        // also note if a vertex is at the
        // boundary
        std::vector<bool>   at_boundary (triangulation.vertices.size(), false);

        for (typename Triangulation<dim,spacedim>::active_line_iterator
             line=triangulation.begin_active_line();
             line != triangulation.end_line(); ++line)
          {
            if (keep_boundary && line->at_boundary())
              {
                at_boundary[line->vertex_index(0)] = true;
                at_boundary[line->vertex_index(1)] = true;
              }

            minimal_length[line->vertex_index(0)]
              = std::min(line->diameter(),
                         minimal_length[line->vertex_index(0)]);
            minimal_length[line->vertex_index(1)]
              = std::min(line->diameter(),
                         minimal_length[line->vertex_index(1)]);
          }


        const unsigned int n_vertices = triangulation.vertices.size();
        Point<spacedim> shift_vector;

        for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
          {
            // ignore this vertex if we
            // whall keep the boundary and
            // this vertex *is* at the
            // boundary
            if (keep_boundary && at_boundary[vertex])
              continue;

            // first compute a random shift
            // vector
            for (unsigned int d=0; d<spacedim; ++d)
              shift_vector(d) = std::rand()*2.0/RAND_MAX-1;

            shift_vector *= factor * minimal_length[vertex] /
                            std::sqrt(shift_vector.square());

            // finally move the vertex
            triangulation.vertices[vertex] += shift_vector;
          }


        // finally correct hanging nodes
        // again. The following is not
        // necessary for 1D
        typename Triangulation<dim,spacedim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              // this lines has children,
              // thus there are restricted
              // nodes
              {
                // not implemented at
                // present for dim=3 or
                // higher
                // TODO: you can steal the code to do this from GridTools::transform()
                Assert (dim<=2, ExcNotImplemented());

                // compute where the common
                // point of the two child
                // lines will lie and reset
                // it to the correct value
                triangulation.vertices[cell->face(face)->child(0)->vertex_index(1)]
                  = (cell->face(face)->vertex(0) +
                     cell->face(face)->vertex(1)) / 2;
              }
      }


      /**
       * Actually delete a cell, or rather all
       * its children, which is the main step for
       * the coarsening process.  This is the
       * dimension dependent part of @p
       * execute_coarsening. The second argument
       * is a vector which gives for each line
       * index the number of cells containing
       * this line. This information is needed to
       * decide whether a refined line may be
       * coarsened or not in 3D. In 1D and 2D
       * this argument is not needed and thus
       * ignored. The same applies for the last
       * argument and quads instead of lines.
       */
      template <int spacedim>
      static
      void
      delete_children (Triangulation<1,spacedim> &triangulation,
                       typename Triangulation<1,spacedim>::cell_iterator &cell,
                       std::vector<unsigned int> &,
                       std::vector<unsigned int> &)
      {
        const unsigned int dim = 1;

        // first we need to reset the
        // neighbor pointers of the
        // neighbors of this cell's
        // children to this cell. This is
        // different for one dimension,
        // since there neighbors can have a
        // refinement level differing from
        // that of this cell's children by
        // more than one level.

        Assert (!cell->child(0)->has_children() && !cell->child(1)->has_children(),
                ExcInternalError());

        // first do it for the cells to the
        // left
        if (cell->neighbor(0).state() == IteratorState::valid)
          if (cell->neighbor(0)->has_children())
            {
              typename Triangulation<dim,spacedim>::cell_iterator
              neighbor = cell->neighbor(0);
              Assert (neighbor->level() == cell->level(), ExcInternalError());

              // right child
              neighbor = neighbor->child(1);
              while (1)
                {
                  Assert (neighbor->neighbor(1) == cell->child(0),
                          ExcInternalError());
                  neighbor->set_neighbor (1, cell);

                  // move on to further
                  // children on the
                  // boundary between this
                  // cell and its neighbor
                  if (neighbor->has_children())
                    neighbor = neighbor->child(1);
                  else
                    break;
                }
            }

        // now do it for the cells to the
        // left
        if (cell->neighbor(1).state() == IteratorState::valid)
          if (cell->neighbor(1)->has_children())
            {
              typename Triangulation<dim,spacedim>::cell_iterator
              neighbor = cell->neighbor(1);
              Assert (neighbor->level() == cell->level(), ExcInternalError());

              // left child
              neighbor = neighbor->child(0);
              while (1)
                {
                  Assert (neighbor->neighbor(0) == cell->child(1),
                          ExcInternalError());
                  neighbor->set_neighbor (0, cell);

                  // move on to further
                  // children on the
                  // boundary between this
                  // cell and its neighbor
                  if (neighbor->has_children())
                    neighbor = neighbor->child(0);
                  else
                    break;
                }
            }


        // delete the vertex which will not
        // be needed anymore. This vertex
        // is the second of the first child
        triangulation.vertices_used[cell->child(0)->vertex_index(1)] = false;

        // invalidate children.  clear user
        // pointers, to avoid that they may
        // appear at unwanted places later
        // on...
        for (unsigned int child=0; child<cell->n_children(); ++child)
          {
            cell->child(child)->clear_user_data();
            cell->child(child)->clear_user_flag();
            cell->child(child)->clear_used_flag();
          }


        // delete pointer to children
        cell->clear_children ();
        cell->clear_user_flag();
      }



      template <int spacedim>
      static
      void
      delete_children (Triangulation<2,spacedim> &triangulation,
                       typename Triangulation<2,spacedim>::cell_iterator &cell,
                       std::vector<unsigned int> &line_cell_count,
                       std::vector<unsigned int> &)
      {
        const unsigned int dim=2;
        const RefinementCase<dim> ref_case=cell->refinement_case();

        Assert(line_cell_count.size()==triangulation.n_raw_lines(), ExcInternalError());

        // vectors to hold all lines which
        // may be deleted
        std::vector<typename Triangulation<dim,spacedim>::line_iterator>
        lines_to_delete(0);

        lines_to_delete.reserve(4*2+4);

        // now we decrease the counters for
        // lines contained in the child
        // cells
        for (unsigned int c=0; c<cell->n_children(); ++c)
          {
            typename Triangulation<dim,spacedim>::cell_iterator
            child=cell->child(c);
            for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
              --line_cell_count[child->line_index(l)];
          }


        // delete the vertex which will not
        // be needed anymore. This vertex
        // is the second of the second line
        // of the first child, if the cell
        // is refined with cut_xy, else there
        // is no inner vertex.
        // additionally delete unneeded inner
        // lines
        if (ref_case==RefinementCase<dim>::cut_xy)
          {
            triangulation.vertices_used[cell->child(0)->line(1)->vertex_index(1)] = false;

            lines_to_delete.push_back(cell->child(0)->line(1));
            lines_to_delete.push_back(cell->child(0)->line(3));
            lines_to_delete.push_back(cell->child(3)->line(0));
            lines_to_delete.push_back(cell->child(3)->line(2));
          }
        else
          {
            unsigned int inner_face_no=ref_case==RefinementCase<dim>::cut_x ? 1 : 3;

            // the inner line will not be
            // used any more
            lines_to_delete.push_back(cell->child(0)->line(inner_face_no));
          }

        // invalidate children
        for (unsigned int child=0; child<cell->n_children(); ++child)
          {
            cell->child(child)->clear_user_data();
            cell->child(child)->clear_user_flag();
            cell->child(child)->clear_used_flag();
          }


        // delete pointer to children
        cell->clear_children ();
        cell->clear_refinement_case();
        cell->clear_user_flag();

        // look at the refinement of outer
        // lines. if nobody needs those
        // anymore we can add them to the
        // list of lines to be deleted.
        for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
          {
            typename Triangulation<dim,spacedim>::line_iterator
            line=cell->line(line_no);

            if (line->has_children())
              {
                // if one of the cell counters is
                // zero, the other has to be as well

                Assert((line_cell_count[line->child_index(0)] == 0 &&
                        line_cell_count[line->child_index(1)] == 0) ||
                       (line_cell_count[line->child_index(0)] > 0 &&
                        line_cell_count[line->child_index(1)] > 0),
                       ExcInternalError());

                if (line_cell_count[line->child_index(0)]==0)
                  {
                    for (unsigned int c=0; c<2; ++c)
                      Assert (!line->child(c)->has_children(),
                              ExcInternalError());

                    // we may delete the line's
                    // children and the middle vertex
                    // as no cell references them
                    // anymore
                    triangulation.vertices_used[line->child(0)->vertex_index(1)] = false;

                    lines_to_delete.push_back(line->child(0));
                    lines_to_delete.push_back(line->child(1));

                    line->clear_children();
                  }
              }
          }

        // finally, delete unneeded lines

        // clear user pointers, to avoid that
        // they may appear at unwanted places
        // later on...
        // same for user flags, then finally
        // delete the lines
        typename std::vector<typename Triangulation<dim,spacedim>::line_iterator>::iterator
        line=lines_to_delete.begin(),
        endline=lines_to_delete.end();
        for (; line!=endline; ++line)
          {
            (*line)->clear_user_data();
            (*line)->clear_user_flag();
            (*line)->clear_used_flag();
          }
      }



      template <int spacedim>
      static
      void
      delete_children (Triangulation<3,spacedim> &triangulation,
                       typename Triangulation<3,spacedim>::cell_iterator &cell,
                       std::vector<unsigned int> &line_cell_count,
                       std::vector<unsigned int> &quad_cell_count)
      {
        const unsigned int dim=3;

        Assert(line_cell_count.size()==triangulation.n_raw_lines(), ExcInternalError());
        Assert(quad_cell_count.size()==triangulation.n_raw_quads(), ExcInternalError());

        // first of all, we store the RefineCase of
        // this cell
        const RefinementCase<dim> ref_case=cell->refinement_case();
        // vectors to hold all lines and quads which
        // may be deleted
        std::vector<typename Triangulation<dim,spacedim>::line_iterator>
        lines_to_delete(0);
        std::vector<typename Triangulation<dim,spacedim>::quad_iterator>
        quads_to_delete(0);

        lines_to_delete.reserve(12*2+6*4+6);
        quads_to_delete.reserve(6*4+12);

        // now we decrease the counters for lines and
        // quads contained in the child cells
        for (unsigned int c=0; c<cell->n_children(); ++c)
          {
            typename Triangulation<dim,spacedim>::cell_iterator
            child=cell->child(c);
            for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
              --line_cell_count[child->line_index(l)];
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              --quad_cell_count[child->quad_index(f)];
          }

        ///////////////////////////////////////
        // delete interior quads and lines and the
        // interior vertex, depending on the
        // refinement case of the cell
        //
        // for append quads and lines: only append
        // them to the list of objects to be deleted

        switch (ref_case)
          {
          case RefinementCase<dim>::cut_x:
            quads_to_delete.push_back(cell->child(0)->face(1));
            break;
          case RefinementCase<dim>::cut_y:
            quads_to_delete.push_back(cell->child(0)->face(3));
            break;
          case RefinementCase<dim>::cut_z:
            quads_to_delete.push_back(cell->child(0)->face(5));
            break;
          case RefinementCase<dim>::cut_xy:
            quads_to_delete.push_back(cell->child(0)->face(1));
            quads_to_delete.push_back(cell->child(0)->face(3));
            quads_to_delete.push_back(cell->child(3)->face(0));
            quads_to_delete.push_back(cell->child(3)->face(2));

            lines_to_delete.push_back(cell->child(0)->line(11));
            break;
          case RefinementCase<dim>::cut_xz:
            quads_to_delete.push_back(cell->child(0)->face(1));
            quads_to_delete.push_back(cell->child(0)->face(5));
            quads_to_delete.push_back(cell->child(3)->face(0));
            quads_to_delete.push_back(cell->child(3)->face(4));

            lines_to_delete.push_back(cell->child(0)->line(5));
            break;
          case RefinementCase<dim>::cut_yz:
            quads_to_delete.push_back(cell->child(0)->face(3));
            quads_to_delete.push_back(cell->child(0)->face(5));
            quads_to_delete.push_back(cell->child(3)->face(2));
            quads_to_delete.push_back(cell->child(3)->face(4));

            lines_to_delete.push_back(cell->child(0)->line(7));
            break;
          case RefinementCase<dim>::cut_xyz:
            quads_to_delete.push_back(cell->child(0)->face(1));
            quads_to_delete.push_back(cell->child(2)->face(1));
            quads_to_delete.push_back(cell->child(4)->face(1));
            quads_to_delete.push_back(cell->child(6)->face(1));

            quads_to_delete.push_back(cell->child(0)->face(3));
            quads_to_delete.push_back(cell->child(1)->face(3));
            quads_to_delete.push_back(cell->child(4)->face(3));
            quads_to_delete.push_back(cell->child(5)->face(3));

            quads_to_delete.push_back(cell->child(0)->face(5));
            quads_to_delete.push_back(cell->child(1)->face(5));
            quads_to_delete.push_back(cell->child(2)->face(5));
            quads_to_delete.push_back(cell->child(3)->face(5));

            lines_to_delete.push_back(cell->child(0)->line(5));
            lines_to_delete.push_back(cell->child(0)->line(7));
            lines_to_delete.push_back(cell->child(0)->line(11));
            lines_to_delete.push_back(cell->child(7)->line(0));
            lines_to_delete.push_back(cell->child(7)->line(2));
            lines_to_delete.push_back(cell->child(7)->line(8));
            // delete the vertex which will not
            // be needed anymore. This vertex
            // is the vertex at the heart of
            // this cell, which is the sixth of
            // the first child
            triangulation.vertices_used[cell->child(0)->vertex_index(7)] = false;
            break;
          default:
            // only remaining case is
            // no_refinement, thus an error
            Assert(false, ExcInternalError());
            break;
          }


        // invalidate children
        for (unsigned int child=0; child<cell->n_children(); ++child)
          {
            cell->child(child)->clear_user_data();
            cell->child(child)->clear_user_flag();

            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                // set flags denoting deviations from
                // standard orientation of faces back
                // to initialization values
                cell->child(child)->set_face_orientation (f, true);
                cell->child(child)->set_face_flip(f,false);
                cell->child(child)->set_face_rotation(f,false);
              }

            cell->child(child)->clear_used_flag();
          }


        // delete pointer to children
        cell->clear_children ();
        cell->clear_refinement_case ();
        cell->clear_user_flag();

        // so far we only looked at inner quads,
        // lines and vertices. Now we have to
        // consider outer ones as well. here, we have
        // to check, whether there are other cells
        // still needing these objects. oherwise we
        // can delete them. first for quads (and
        // their inner lines).

        for (unsigned int quad_no=0; quad_no<GeometryInfo<dim>::faces_per_cell; ++quad_no)
          {
            typename Triangulation<dim,spacedim>::quad_iterator
            quad=cell->face(quad_no);

            Assert((GeometryInfo<dim>::face_refinement_case(ref_case,quad_no) && quad->has_children()) ||
                   GeometryInfo<dim>::face_refinement_case(ref_case,quad_no)==RefinementCase<dim-1>::no_refinement,
                   ExcInternalError());

            switch (quad->refinement_case())
              {
              case RefinementCase<dim-1>::no_refinement:
                // nothing to do as the quad
                // is not refined
                break;
              case RefinementCase<dim-1>::cut_x:
              case RefinementCase<dim-1>::cut_y:
              {
                // if one of the cell counters is
                // zero, the other has to be as
                // well
                Assert((quad_cell_count[quad->child_index(0)] == 0 &&
                        quad_cell_count[quad->child_index(1)] == 0) ||
                       (quad_cell_count[quad->child_index(0)] > 0 &&
                        quad_cell_count[quad->child_index(1)] > 0),
                       ExcInternalError());
                // it might be, that the quad is
                // refined twice anisotropically,
                // first check, whether we may
                // delete possible grand_children
                unsigned int deleted_grandchildren=0;
                unsigned int number_of_child_refinements=0;

                for (unsigned int c=0; c<2; ++c)
                  if (quad->child(c)->has_children())
                    {
                      ++number_of_child_refinements;
                      // if one of the cell counters is
                      // zero, the other has to be as
                      // well
                      Assert((quad_cell_count[quad->child(c)->child_index(0)] == 0 &&
                              quad_cell_count[quad->child(c)->child_index(1)] == 0) ||
                             (quad_cell_count[quad->child(c)->child_index(0)] > 0 &&
                              quad_cell_count[quad->child(c)->child_index(1)] > 0),
                             ExcInternalError());
                      if (quad_cell_count[quad->child(c)->child_index(0)]==0)
                        {
                          // Assert, that the two
                          // anisotropic
                          // refinements add up to
                          // isotropic refinement
                          Assert(quad->refinement_case()+quad->child(c)->refinement_case()==RefinementCase<dim>::cut_xy,
                                 ExcInternalError());
                          // we may delete the
                          // quad's children and
                          // the inner line as no
                          // cell references them
                          // anymore
                          quads_to_delete.push_back(quad->child(c)->child(0));
                          quads_to_delete.push_back(quad->child(c)->child(1));
                          if (quad->child(c)->refinement_case()==RefinementCase<2>::cut_x)
                            lines_to_delete.push_back(quad->child(c)->child(0)->line(1));
                          else
                            lines_to_delete.push_back(quad->child(c)->child(0)->line(3));
                          quad->child(c)->clear_children();
                          quad->child(c)->clear_refinement_case();
                          ++deleted_grandchildren;
                        }
                    }
                // if no grandchildren are left, we
                // may as well delete the
                // refinement of the inner line
                // between our children and the
                // corresponding vertex
                if (number_of_child_refinements>0 &&
                    deleted_grandchildren==number_of_child_refinements)
                  {
                    typename Triangulation<dim,spacedim>::line_iterator
                    middle_line;
                    if (quad->refinement_case()==RefinementCase<2>::cut_x)
                      middle_line=quad->child(0)->line(1);
                    else
                      middle_line=quad->child(0)->line(3);

                    lines_to_delete.push_back(middle_line->child(0));
                    lines_to_delete.push_back(middle_line->child(1));
                    triangulation.vertices_used[middle_vertex_index<dim,spacedim>(middle_line)]
                      = false;
                    middle_line->clear_children();
                  }

                // now consider the direct children
                // of the given quad
                if (quad_cell_count[quad->child_index(0)]==0)
                  {
                    // we may delete the quad's
                    // children and the inner line
                    // as no cell references them
                    // anymore
                    quads_to_delete.push_back(quad->child(0));
                    quads_to_delete.push_back(quad->child(1));
                    if (quad->refinement_case()==RefinementCase<2>::cut_x)
                      lines_to_delete.push_back(quad->child(0)->line(1));
                    else
                      lines_to_delete.push_back(quad->child(0)->line(3));

                    // if the counters just dropped
                    // to zero, otherwise the
                    // children would have been
                    // deleted earlier, then this
                    // cell's children must have
                    // contained the anisotropic
                    // quad children. thus, if
                    // those have again anisotropic
                    // children, which are in
                    // effect isotropic children of
                    // the original quad, those are
                    // still needed by a
                    // neighboring cell and we
                    // cannot delete them. instead,
                    // we have to reset this quad's
                    // refine case to isotropic and
                    // set the children
                    // accordingly.
                    if (quad->child(0)->has_children())
                      if (quad->refinement_case()==RefinementCase<2>::cut_x)
                        {
                          // now evereything is
                          // quite complicated. we
                          // have the children
                          // numbered according to
                          //
                          // *---*---*
                          // |n+1|m+1|
                          // *---*---*
                          // | n | m |
                          // *---*---*
                          //
                          // from the original
                          // anisotropic
                          // refinement. we have to
                          // reorder them as
                          //
                          // *---*---*
                          // | m |m+1|
                          // *---*---*
                          // | n |n+1|
                          // *---*---*
                          //
                          // for isotropic refinement.
                          //
                          // this is a bit ugly, of
                          // course: loop over all
                          // cells on all levels
                          // and look for faces n+1
                          // (switch_1) and m
                          // (switch_2).
                          const typename Triangulation<dim,spacedim>::quad_iterator
                          switch_1=quad->child(0)->child(1),
                          switch_2=quad->child(1)->child(0);

                          Assert(!switch_1->has_children(), ExcInternalError());
                          Assert(!switch_2->has_children(), ExcInternalError());

                          const int switch_1_index=switch_1->index();
                          const int switch_2_index=switch_2->index();
                          for (unsigned int l=0; l<triangulation.levels.size(); ++l)
                            for (unsigned int h=0; h<triangulation.levels[l]->cells.cells.size(); ++h)
                              for (unsigned int q=0; q<GeometryInfo<dim>::faces_per_cell; ++q)
                                {
                                  const int index=triangulation.levels[l]->cells.cells[h].face(q);
                                  if (index==switch_1_index)
                                    triangulation.levels[l]->cells.cells[h].set_face(q,switch_2_index);
                                  else if (index==switch_2_index)
                                    triangulation.levels[l]->cells.cells[h].set_face(q,switch_1_index);
                                }
                          // now we have to copy
                          // all information of the
                          // two quads
                          const int switch_1_lines[4]=
                          {
                            static_cast<signed int>(switch_1->line_index(0)),
                            static_cast<signed int>(switch_1->line_index(1)),
                            static_cast<signed int>(switch_1->line_index(2)),
                            static_cast<signed int>(switch_1->line_index(3))
                          };
                          const bool switch_1_line_orientations[4]=
                          {
                            switch_1->line_orientation(0),
                            switch_1->line_orientation(1),
                            switch_1->line_orientation(2),
                            switch_1->line_orientation(3)
                          };
                          const types::boundary_id switch_1_boundary_indicator=switch_1->boundary_indicator();
                          const unsigned int switch_1_user_index=switch_1->user_index();
                          const bool switch_1_user_flag=switch_1->user_flag_set();

                          switch_1->set(internal::Triangulation::TriaObject<2>(switch_2->line_index(0),
                                                                               switch_2->line_index(1),
                                                                               switch_2->line_index(2),
                                                                               switch_2->line_index(3)));
                          switch_1->set_line_orientation(0, switch_2->line_orientation(0));
                          switch_1->set_line_orientation(1, switch_2->line_orientation(1));
                          switch_1->set_line_orientation(2, switch_2->line_orientation(2));
                          switch_1->set_line_orientation(3, switch_2->line_orientation(3));
                          switch_1->set_boundary_indicator(switch_2->boundary_indicator());
                          switch_1->set_manifold_id(switch_2->manifold_id());
                          switch_1->set_user_index(switch_2->user_index());
                          if (switch_2->user_flag_set())
                            switch_1->set_user_flag();
                          else
                            switch_1->clear_user_flag();

                          switch_2->set(internal::Triangulation::TriaObject<2>(switch_1_lines[0],
                                                                               switch_1_lines[1],
                                                                               switch_1_lines[2],
                                                                               switch_1_lines[3]));
                          switch_2->set_line_orientation(0, switch_1_line_orientations[0]);
                          switch_2->set_line_orientation(1, switch_1_line_orientations[1]);
                          switch_2->set_line_orientation(2, switch_1_line_orientations[2]);
                          switch_2->set_line_orientation(3, switch_1_line_orientations[3]);
                          switch_2->set_boundary_indicator(switch_1_boundary_indicator);
                          switch_2->set_manifold_id(switch_1->manifold_id());
                          switch_2->set_user_index(switch_1_user_index);
                          if (switch_1_user_flag)
                            switch_2->set_user_flag();
                          else
                            switch_2->clear_user_flag();

                          const unsigned int child_0=quad->child(0)->child_index(0);
                          const unsigned int child_2=quad->child(1)->child_index(0);
                          quad->clear_children();
                          quad->clear_refinement_case();
                          quad->set_refinement_case(RefinementCase<2>::cut_xy);
                          quad->set_children(0,child_0);
                          quad->set_children(2,child_2);
                          std::swap(quad_cell_count[child_0+1],quad_cell_count[child_2]);
                        }
                      else
                        {
                          // the face was refined
                          // with cut_y, thus the
                          // children are already
                          // in correct order. we
                          // only have to set them
                          // correctly, deleting
                          // the indirection of two
                          // anisotropic refinement
                          // and going directly
                          // from the quad to
                          // isotropic children
                          const unsigned int child_0=quad->child(0)->child_index(0);
                          const unsigned int child_2=quad->child(1)->child_index(0);
                          quad->clear_children();
                          quad->clear_refinement_case();
                          quad->set_refinement_case(RefinementCase<2>::cut_xy);
                          quad->set_children(0,child_0);
                          quad->set_children(2,child_2);
                        }
                    else
                      {
                        quad->clear_children();
                        quad->clear_refinement_case();
                      }


                  }
                break;
              }
              case RefinementCase<dim-1>::cut_xy:
              {
                // if one of the cell counters is
                // zero, the others have to be as
                // well

                Assert((quad_cell_count[quad->child_index(0)] == 0 &&
                        quad_cell_count[quad->child_index(1)] == 0 &&
                        quad_cell_count[quad->child_index(2)] == 0 &&
                        quad_cell_count[quad->child_index(3)] == 0) ||
                       (quad_cell_count[quad->child_index(0)] > 0 &&
                        quad_cell_count[quad->child_index(1)] > 0 &&
                        quad_cell_count[quad->child_index(2)] > 0 &&
                        quad_cell_count[quad->child_index(3)] > 0),
                       ExcInternalError());

                if (quad_cell_count[quad->child_index(0)]==0)
                  {
                    // we may delete the quad's
                    // children, the inner lines
                    // and the middle vertex as no
                    // cell references them anymore
                    lines_to_delete.push_back(quad->child(0)->line(1));
                    lines_to_delete.push_back(quad->child(3)->line(0));
                    lines_to_delete.push_back(quad->child(0)->line(3));
                    lines_to_delete.push_back(quad->child(3)->line(2));

                    for (unsigned int child=0; child<quad->n_children(); ++child)
                      quads_to_delete.push_back(quad->child(child));

                    triangulation.vertices_used[quad->child(0)->vertex_index(3)] = false;

                    quad->clear_children();
                    quad->clear_refinement_case();
                  }
              }
              break;

              default:
                Assert(false, ExcInternalError());
                break;
              }

          }

        // now we repeat a similar procedure
        // for the outer lines of this cell.

        // if in debug mode: check that each
        // of the lines for which we consider
        // deleting the children in fact has
        // children (the bits/coarsening_3d
        // test tripped over this initially)
        for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
          {
            typename Triangulation<dim,spacedim>::line_iterator
            line=cell->line(line_no);

            Assert((GeometryInfo<dim>::line_refinement_case(ref_case,line_no) && line->has_children()) ||
                   GeometryInfo<dim>::line_refinement_case(ref_case,line_no)==RefinementCase<1>::no_refinement,
                   ExcInternalError());

            if (line->has_children())
              {
                // if one of the cell counters is
                // zero, the other has to be as well

                Assert((line_cell_count[line->child_index(0)] == 0 &&
                        line_cell_count[line->child_index(1)] == 0) ||
                       (line_cell_count[line->child_index(0)] > 0 &&
                        line_cell_count[line->child_index(1)] > 0),
                       ExcInternalError());

                if (line_cell_count[line->child_index(0)]==0)
                  {
                    for (unsigned int c=0; c<2; ++c)
                      Assert (!line->child(c)->has_children(),
                              ExcInternalError());

                    // we may delete the line's
                    // children and the middle vertex
                    // as no cell references them
                    // anymore
                    triangulation.vertices_used[line->child(0)->vertex_index(1)] = false;

                    lines_to_delete.push_back(line->child(0));
                    lines_to_delete.push_back(line->child(1));

                    line->clear_children();
                  }
              }
          }

        // finally, delete unneeded quads and lines

        // clear user pointers, to avoid that
        // they may appear at unwanted places
        // later on...
        // same for user flags, then finally
        // delete the quads and lines
        typename std::vector<typename Triangulation<dim,spacedim>::line_iterator>::iterator
        line=lines_to_delete.begin(),
        endline=lines_to_delete.end();
        for (; line!=endline; ++line)
          {
            (*line)->clear_user_data();
            (*line)->clear_user_flag();
            (*line)->clear_used_flag();
          }

        typename std::vector<typename Triangulation<dim,spacedim>::quad_iterator>::iterator
        quad=quads_to_delete.begin(),
        endquad=quads_to_delete.end();
        for (; quad!=endquad; ++quad)
          {
            (*quad)->clear_user_data();
            (*quad)->clear_children();
            (*quad)->clear_refinement_case();
            (*quad)->clear_user_flag();
            (*quad)->clear_used_flag();
          }
      }


      /**
       * Create the children of a 2d
       * cell. The arguments indicate
       * the next free spots in the
       * vertices, lines, and cells
       * arrays.
       *
       * The faces of the cell have to
       * be refined already, whereas
       * the inner lines in 2D will be
       * created in this
       * function. Therefore iterator
       * pointers into the vectors of
       * lines, quads and cells have to
       * be passed, which point at (or
       * "before") the reserved space.
       */
      template <int spacedim>
      static
      void
      create_children (Triangulation<2,spacedim> &triangulation,
                       unsigned int &next_unused_vertex,
                       typename Triangulation<2,spacedim>::raw_line_iterator &next_unused_line,
                       typename Triangulation<2,spacedim>::raw_cell_iterator &next_unused_cell,
                       typename Triangulation<2,spacedim>::cell_iterator &cell)
      {
        const unsigned int dim=2;
        // clear refinement flag
        const RefinementCase<dim> ref_case=cell->refine_flag_set();
        cell->clear_refine_flag ();

        /* For the refinement process: since we go the levels up from the lowest, there
           are (unlike above) only two possibilities: a neighbor cell is on the same
           level or one level up (in both cases, it may or may not be refined later on,
           but we don't care here).

           First:
           Set up an array of the 3x3 vertices, which are distributed on the cell
           (the array consists of indices into the @p{vertices} std::vector

           2--7--3
           |  |  |
           4--8--5
           |  |  |
           0--6--1

           note: in case of cut_x or cut_y not all these vertices are needed for the new
           cells

           Second:
           Set up an array of the new lines (the array consists of iterator pointers
           into the lines arrays)

           .-6-.-7-.         The directions are:  .->-.->-.
           1   9   3                              ^   ^   ^
           .-10.11-.                             .->-.->-.
           0   8   2                              ^   ^   ^
           .-4-.-5-.                              .->-.->-.

           cut_x:
           .-4-.-5-.
           |   |   |
           0   6   1
           |   |   |
           .-2-.-3-.

           cut_y:
           .---5---.
           1       3
           .---6---.
           0       2
           .---4---.


           Third:
           Set up an array of neighbors:

           6  7
           .--.--.
           1|  |  |3
           .--.--.
           0|  |  |2
           .--.--.
           4   5

           We need this array for two reasons: first to get the lines which will
           bound the four subcells (if the neighboring cell is refined, these
           lines already exist), and second to update neighborship information.
           Since if a neighbor is not refined, its neighborship record only
           points to the present, unrefined, cell rather than the children we
           are presently creating, we only need the neighborship information
           if the neighbor cells are refined. In all other cases, we store
           the unrefined neighbor address

           We also need for every neighbor (if refined) which number among its
           neighbors the present (unrefined) cell has, since that number is to
           be replaced and because that also is the number of the subline which
           will be the interface between that neighbor and the to be created cell.
           We will store this number (between 0 and 3) in the field
           @p{neighbors_neighbor}.

           It would be sufficient to use the children of the common line to the
           neighbor, if we only wanted to get the new sublines and the new vertex,
           but because we need to update the neighborship information of the
           two refined subcells of the neighbor, we need to search these anyway.

           Convention:
           The created children are numbered like this:

           .--.--.
           |2 . 3|
           .--.--.
           |0 | 1|
           .--.--.
        */
        // collect the
        // indices of the
        // eight
        // surrounding
        // vertices
        //   2--7--3
        //   |  |  |
        //   4--9--5
        //   |  |  |
        //   0--6--1
        int new_vertices[9];
        for (unsigned int vertex_no=0; vertex_no<4; ++vertex_no)
          new_vertices[vertex_no]=cell->vertex_index(vertex_no);
        for (unsigned int line_no=0; line_no<4; ++line_no)
          if (cell->line(line_no)->has_children())
            new_vertices[4+line_no]=cell->line(line_no)->child(0)->vertex_index(1);

        if (ref_case==RefinementCase<dim>::cut_xy)
          {

            // find the next
            // unused vertex and
            // allocate it for
            // the new vertex we
            // need here
            while (triangulation.vertices_used[next_unused_vertex] == true)
              ++next_unused_vertex;
            Assert (next_unused_vertex < triangulation.vertices.size(),
                    ExcTooFewVerticesAllocated());
            triangulation.vertices_used[next_unused_vertex] = true;

            new_vertices[8] = next_unused_vertex;

            // if this quad lives
            // in 2d, then we can
            // compute the new
            // central vertex
            // location just from
            // the surrounding
            // ones. If this is
            // not the case, then
            // we need to ask a
            // boundary object
            if (dim == spacedim)
              {
                // triangulation.vertices[next_unused_vertex] = new_point;
                triangulation.vertices[next_unused_vertex] = cell->center(true);

                // if the user_flag is set, i.e. if the
                // cell is at the boundary, use a
                // different calculation of the middle
                // vertex here. this is of advantage, if
                // the boundary is strongly curved and
                // the cell has a high aspect ratio. this
                // can happen for example, if it was
                // refined anisotropically before.
                if (cell->user_flag_set())
                  {
                    // first reset the user_flag
                    cell->clear_user_flag();
                    // the user flag indicates: at least
                    // one face is at the boundary. if it
                    // is only one, set the new middle
                    // vertex in a different way to avoid
                    // some mis-shaped elements if the
                    // new point on the boundary is not
                    // where we expect it, especially if
                    // it is to far inside the current
                    // cell
                    unsigned int boundary_face=GeometryInfo<dim>::faces_per_cell;
                    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                      if (cell->face(face)->at_boundary())
                        {
                          if (boundary_face == GeometryInfo<dim>::faces_per_cell)
                            // no boundary face found so
                            // far, so set it now
                            boundary_face=face;
                          else
                            // there is another boundary
                            // face, so reset boundary_face to
                            // invalid value as a flag to
                            // do nothing in the following
                            boundary_face=GeometryInfo<dim>::faces_per_cell+1;
                        }

                    if (boundary_face<GeometryInfo<dim>::faces_per_cell)
                      // reset the cell's middle vertex to the middle
                      // of the straight connection between the new
                      // points on this face and on the opposite face,
                      // as returned by the underlying manifold
                      // object.
                      {
                        std::vector<Point<spacedim> > ps(2);
                        std::vector<double> ws(2, 0.5);
                        ps[0] = cell->face(boundary_face)
                                ->child(0)->vertex(1);
                        ps[1] = cell->face(GeometryInfo<dim>
                                           ::opposite_face[boundary_face])
                                ->child(0)->vertex(1);
                        Quadrature<spacedim> qs(ps,ws);
                        triangulation.vertices[next_unused_vertex]
                          = cell->get_manifold().get_new_point(qs);
                      }
                  }
              }
            else
              {
                // if this quad lives in a higher dimensional space
                // then we don't need to worry if it is at the
                // boundary of the manifold -- we always have to use
                // the boundary object anyway; so ignore whether the
                // user flag is set or not
                cell->clear_user_flag();

                // An assert to make sure that the static_cast in the
                // next line has the chance to give reasonable
                // results.
                Assert(cell->material_id()<= std::numeric_limits<types::material_id>::max(),
                       ExcIndexRange(cell->material_id(),0,std::numeric_limits<types::material_id>::max()));

                // new vertex is placed on the surface according to
                // the information stored in the boundary class
                triangulation.vertices[next_unused_vertex] =
                  cell->center(true);
              }
          }


        // Now the lines:
        typename Triangulation<dim,spacedim>::raw_line_iterator new_lines[12];
        unsigned int lmin=8;
        unsigned int lmax=12;
        if (ref_case!=RefinementCase<dim>::cut_xy)
          {
            lmin=6;
            lmax=7;
          }

        for (unsigned int l=lmin; l<lmax; ++l)
          {
            while (next_unused_line->used() == true)
              ++next_unused_line;
            new_lines[l] = next_unused_line;
            ++next_unused_line;

            Assert (new_lines[l]->used() == false,
                    ExcCellShouldBeUnused());
          }

        if (ref_case==RefinementCase<dim>::cut_xy)
          {
            //   .-6-.-7-.
            //   1   9   3
            //   .-10.11-.
            //   0   8   2
            //   .-4-.-5-.

            // lines 0-7 already exist, create only the four interior
            // lines 8-11
            unsigned int l=0;
            for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
              for (unsigned int c=0; c<2; ++c, ++l)
                new_lines[l]=cell->line(face_no)->child(c);
            Assert(l==8, ExcInternalError());

            new_lines[8] ->set (internal::Triangulation::
                                TriaObject<1>(new_vertices[6], new_vertices[8]));
            new_lines[9] ->set (internal::Triangulation::
                                TriaObject<1>(new_vertices[8], new_vertices[7]));
            new_lines[10]->set (internal::Triangulation::
                                TriaObject<1>(new_vertices[4], new_vertices[8]));
            new_lines[11]->set (internal::Triangulation::
                                TriaObject<1>(new_vertices[8], new_vertices[5]));
          }
        else if (ref_case==RefinementCase<dim>::cut_x)
          {
            //   .-4-.-5-.
            //   |   |   |
            //   0   6   1
            //   |   |   |
            //   .-2-.-3-.
            new_lines[0]=cell->line(0);
            new_lines[1]=cell->line(1);
            new_lines[2]=cell->line(2)->child(0);
            new_lines[3]=cell->line(2)->child(1);
            new_lines[4]=cell->line(3)->child(0);
            new_lines[5]=cell->line(3)->child(1);
            new_lines[6]->set (internal::Triangulation::
                               TriaObject<1>(new_vertices[6], new_vertices[7]));
          }
        else
          {
            Assert(ref_case==RefinementCase<dim>::cut_y, ExcInternalError());
            //   .---5---.
            //   1       3
            //   .---6---.
            //   0       2
            //   .---4---.
            new_lines[0]=cell->line(0)->child(0);
            new_lines[1]=cell->line(0)->child(1);
            new_lines[2]=cell->line(1)->child(0);
            new_lines[3]=cell->line(1)->child(1);
            new_lines[4]=cell->line(2);
            new_lines[5]=cell->line(3);
            new_lines[6]->set (internal::Triangulation::
                               TriaObject<1>(new_vertices[4], new_vertices[5]));
          }

        for (unsigned int l=lmin; l<lmax; ++l)
          {
            new_lines[l]->set_used_flag();
            new_lines[l]->clear_user_flag();
            new_lines[l]->clear_user_data();
            new_lines[l]->clear_children();
            // interior line
            new_lines[l]->set_boundary_indicator(numbers::internal_face_boundary_id);
            new_lines[l]->set_manifold_id(cell->manifold_id());
          }

        // Now add the four (two)
        // new cells!
        typename Triangulation<dim,spacedim>::raw_cell_iterator
        subcells[GeometryInfo<dim>::max_children_per_cell];
        while (next_unused_cell->used() == true)
          ++next_unused_cell;

        const unsigned int n_children=
          GeometryInfo<dim>::n_children(ref_case);
        for (unsigned int i=0; i<n_children; ++i)
          {
            Assert (next_unused_cell->used() == false,
                    ExcCellShouldBeUnused());
            subcells[i] = next_unused_cell;
            ++next_unused_cell;
            if (i%2==1 && i<n_children-1)
              while (next_unused_cell->used() == true)
                ++next_unused_cell;
          }

        if (ref_case==RefinementCase<dim>::cut_xy)
          {
            // children:
            //   .--.--.
            //   |2 . 3|
            //   .--.--.
            //   |0 | 1|
            //   .--.--.
            // lines:
            //   .-6-.-7-.
            //   1   9   3
            //   .-10.11-.
            //   0   8   2
            //   .-4-.-5-.
            subcells[0]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[0]->index(),
                                            new_lines[8]->index(),
                                            new_lines[4]->index(),
                                            new_lines[10]->index()));
            subcells[1]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[8]->index(),
                                            new_lines[2]->index(),
                                            new_lines[5]->index(),
                                            new_lines[11]->index()));
            subcells[2]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[1]->index(),
                                            new_lines[9]->index(),
                                            new_lines[10]->index(),
                                            new_lines[6]->index()));
            subcells[3]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[9]->index(),
                                            new_lines[3]->index(),
                                            new_lines[11]->index(),
                                            new_lines[7]->index()));
          }
        else if (ref_case==RefinementCase<dim>::cut_x)
          {
            // children:
            //   .--.--.
            //   |  .  |
            //   .0 . 1.
            //   |  |  |
            //   .--.--.
            // lines:
            //   .-4-.-5-.
            //   |   |   |
            //   0   6   1
            //   |   |   |
            //   .-2-.-3-.
            subcells[0]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[0]->index(),
                                            new_lines[6]->index(),
                                            new_lines[2]->index(),
                                            new_lines[4]->index()));
            subcells[1]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[6]->index(),
                                            new_lines[1]->index(),
                                            new_lines[3]->index(),
                                            new_lines[5]->index()));
          }
        else
          {
            Assert(ref_case==RefinementCase<dim>::cut_y, ExcInternalError());
            // children:
            //   .-----.
            //   |  1  |
            //   .-----.
            //   |  0  |
            //   .-----.
            // lines:
            //   .---5---.
            //   1       3
            //   .---6---.
            //   0       2
            //   .---4---.
            subcells[0]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[0]->index(),
                                            new_lines[2]->index(),
                                            new_lines[4]->index(),
                                            new_lines[6]->index()));
            subcells[1]->set (internal::Triangulation::
                              TriaObject<2>(new_lines[1]->index(),
                                            new_lines[3]->index(),
                                            new_lines[6]->index(),
                                            new_lines[5]->index()));
          }

        types::subdomain_id subdomainid = cell->subdomain_id();

        for (unsigned int i=0; i<n_children; ++i)
          {
            subcells[i]->set_used_flag();
            subcells[i]->clear_refine_flag();
            subcells[i]->clear_user_flag();
            subcells[i]->clear_user_data();
            subcells[i]->clear_children();
            // inherit material
            // properties
            subcells[i]->set_material_id (cell->material_id());
            subcells[i]->set_manifold_id (cell->manifold_id());
            subcells[i]->set_subdomain_id (subdomainid);

            if (i%2==0)
              subcells[i]->set_parent (cell->index ());
          }



        // set child index for
        // even children children
        // i=0,2 (0)
        for (unsigned int i=0; i<n_children/2; ++i)
          cell->set_children (2*i, subcells[2*i]->index());
        // set the refine case
        cell->set_refinement_case(ref_case);

        // note that the
        // refinement flag was
        // already cleared at the
        // beginning of this function

        if (dim < spacedim)
          for (unsigned int c=0; c<n_children; ++c)
            cell->child(c)->set_direction_flag (cell->direction_flag());

      }



      /**
       * A function that performs the
       * refinement of a triangulation in 1d.
       */
      template <int spacedim>
      static
      typename Triangulation<1,spacedim>::DistortedCellList
      execute_refinement (Triangulation<1,spacedim> &triangulation,
                          const bool /*check_for_distorted_cells*/)
      {
        const unsigned int dim = 1;

        // check whether a new level is needed we have to check for
        // this on the highest level only (on this, all used cells are
        // also active, so we only have to check for this)
        {
          typename Triangulation<dim,spacedim>::raw_cell_iterator
          cell = triangulation.begin_active (triangulation.levels.size()-1),
          endc = triangulation.end();
          for (; cell != endc; ++cell)
            if (cell->used())
              if (cell->refine_flag_set())
                {
                  triangulation.levels
                  .push_back (new internal::Triangulation::TriaLevel<dim>);
                  break;
                }
        }


        // check how much space is needed on every level we need not
        // check the highest level since either - on the highest level
        // no cells are flagged for refinement - there are, but
        // prepare_refinement added another empty level
        unsigned int needed_vertices = 0;
        for (int level=triangulation.levels.size()-2; level>=0; --level)
          {
            // count number of flagged
            // cells on this level
            unsigned int flagged_cells = 0;
            typename Triangulation<dim,spacedim>::active_cell_iterator
            acell = triangulation.begin_active(level),
            aendc = triangulation.begin_active(level+1);
            for (; acell!=aendc; ++acell)
              if (acell->refine_flag_set())
                ++flagged_cells;

            // count number of used cells
            // on the next higher level
            const unsigned int used_cells
              =  std::count_if (triangulation.levels[level+1]->cells.used.begin(),
                                triangulation.levels[level+1]->cells.used.end(),
                                std::bind2nd (std::equal_to<bool>(), true));

            // reserve space for the used_cells cells already existing
            // on the next higher level as well as for the
            // 2*flagged_cells that will be created on that level
            triangulation.levels[level+1]
            ->reserve_space(used_cells+
                            GeometryInfo<1>::max_children_per_cell *
                            flagged_cells,
                            1,
                            spacedim);
            // reserve space for 2*flagged_cells new lines on the next
            // higher level
            triangulation.levels[level+1]->cells
            .reserve_space (GeometryInfo<1>::max_children_per_cell *
                            flagged_cells,
                            0);

            needed_vertices += flagged_cells;
          }

        // add to needed vertices how many
        // vertices are already in use
        needed_vertices += std::count_if (triangulation.vertices_used.begin(),
                                          triangulation.vertices_used.end(),
                                          std::bind2nd (std::equal_to<bool>(),
                                                        true));
        // if we need more vertices: create them, if not: leave the
        // array as is, since shrinking is not really possible because
        // some of the vertices at the end may be in use
        if (needed_vertices > triangulation.vertices.size())
          {
            triangulation.vertices.resize (needed_vertices,
                                           Point<spacedim>());
            triangulation.vertices_used.resize (needed_vertices, false);
          }


        // Do REFINEMENT on every level; exclude highest level as
        // above

        // index of next unused vertex
        unsigned int next_unused_vertex = 0;

        for (int level=triangulation.levels.size()-2; level>=0; --level)
          {
            typename Triangulation<dim,spacedim>::active_cell_iterator
            cell = triangulation.begin_active(level),
            endc = triangulation.begin_active(level+1);

            typename Triangulation<dim,spacedim>::raw_cell_iterator
            next_unused_cell = triangulation.begin_raw (level+1);

            for (; (cell!=endc) && (cell->level()==level); ++cell)
              if (cell->refine_flag_set())
                {
                  // clear refinement flag
                  cell->clear_refine_flag ();

                  // search for next unused
                  // vertex
                  while (triangulation.vertices_used[next_unused_vertex] == true)
                    ++next_unused_vertex;
                  Assert (next_unused_vertex < triangulation.vertices.size(),
                          ExcTooFewVerticesAllocated());

                  // Now we always ask the cell itself where to put
                  // the new point. The cell in turn will query the
                  // manifold object internally.
                  triangulation.vertices[next_unused_vertex] =
                    cell->center(true);

                  triangulation.vertices_used[next_unused_vertex] = true;

                  // search for next two unused cell (++ takes care of
                  // the end of the vector)
                  typename Triangulation<dim,spacedim>::raw_cell_iterator
                  first_child,
                  second_child;
                  while (next_unused_cell->used() == true)
                    ++next_unused_cell;
                  first_child = next_unused_cell;
                  first_child->set_used_flag ();
                  first_child->clear_user_data ();
                  ++next_unused_cell;
                  Assert (next_unused_cell->used() == false,
                          ExcCellShouldBeUnused());
                  second_child = next_unused_cell;
                  second_child->set_used_flag ();
                  second_child->clear_user_data ();

                  types::subdomain_id subdomainid = cell->subdomain_id();

                  // insert first child
                  cell->set_children (0, first_child->index());
                  first_child->clear_children ();
                  first_child->set (internal::Triangulation
                                    ::TriaObject<1> (cell->vertex_index(0),
                                                     next_unused_vertex));
                  first_child->set_material_id (cell->material_id());
                  first_child->set_manifold_id (cell->manifold_id());
                  first_child->set_subdomain_id (subdomainid);
                  first_child->set_direction_flag (cell->direction_flag());

                  first_child->set_parent (cell->index ());

                  // Set manifold id of the right face. Only do this
                  // on the first child.
                  first_child->face(1)->set_manifold_id(cell->manifold_id());

                  // reset neighborship info (refer to
                  // internal::Triangulation::TriaLevel<0> for
                  // details)
                  first_child->set_neighbor (1, second_child);
                  if (cell->neighbor(0).state() != IteratorState::valid)
                    first_child->set_neighbor (0, cell->neighbor(0));
                  else if (cell->neighbor(0)->active())
                    {
                      // since the neighbors level is always <=level,
                      // if the cell is active, then there are no
                      // cells to the left which may want to know
                      // about this new child cell.
                      Assert (cell->neighbor (0)->level () <= cell->level (),
                              ExcInternalError ());
                      first_child->set_neighbor (0, cell->neighbor(0));
                    }
                  else
                    // left neighbor is refined
                    {
                      // set neighbor to cell on same level
                      const unsigned int nbnb = cell->neighbor_of_neighbor (0);
                      first_child->set_neighbor (0, cell->neighbor(0)->child(nbnb));

                      // reset neighbor info of all right descendant
                      // of the left neighbor of cell
                      typename Triangulation<dim,spacedim>::cell_iterator
                      left_neighbor = cell->neighbor(0);
                      while (left_neighbor->has_children())
                        {
                          left_neighbor = left_neighbor->child(nbnb);
                          left_neighbor->set_neighbor (nbnb, first_child);
                        }
                    }

                  // insert second child
                  second_child->clear_children ();
                  second_child->set (internal::Triangulation
                                     ::TriaObject<1>(next_unused_vertex,
                                                     cell->vertex_index(1)));
                  second_child->set_neighbor (0, first_child);
                  second_child->set_material_id (cell->material_id());
                  second_child->set_manifold_id (cell->manifold_id());
                  second_child->set_subdomain_id (subdomainid);
                  second_child->set_direction_flag (cell->direction_flag());

                  if (cell->neighbor(1).state() != IteratorState::valid)
                    second_child->set_neighbor (1, cell->neighbor(1));
                  else if (cell->neighbor(1)->active())
                    {
                      Assert (cell->neighbor (1)->level () <= cell->level (),
                              ExcInternalError ());
                      second_child->set_neighbor (1, cell->neighbor(1));
                    }
                  else
                    // right neighbor is refined same as above
                    {
                      const unsigned int nbnb = cell->neighbor_of_neighbor (1);
                      second_child->set_neighbor (1, cell->neighbor(1)->child(nbnb));

                      typename Triangulation<dim,spacedim>::cell_iterator
                      right_neighbor = cell->neighbor(1);
                      while (right_neighbor->has_children())
                        {
                          right_neighbor = right_neighbor->child(nbnb);
                          right_neighbor->set_neighbor (nbnb, second_child);
                        }
                    }
                }
          }

        // in 1d, we can not have distorted children unless the parent
        // was already distorted (that is because we don't use
        // boundary information for 1d triangulations). so return an
        // empty list
        return typename Triangulation<1,spacedim>::DistortedCellList();
      }


      /**
       * A function that performs the refinement of a triangulation in
       * 2d.
       */
      template <int spacedim>
      static
      typename Triangulation<2,spacedim>::DistortedCellList
      execute_refinement (Triangulation<2,spacedim> &triangulation,
                          const bool check_for_distorted_cells)
      {
        const unsigned int dim = 2;

        // check whether a new level is needed we have to check for
        // this on the highest level only (on this, all used cells are
        // also active, so we only have to check for this)
        if (true)
          {
            typename Triangulation<dim,spacedim>::raw_cell_iterator
            cell = triangulation.begin_active (triangulation.levels.size()-1),
            endc = triangulation.end();
            for (; cell != endc; ++cell)
              if (cell->used())
                if (cell->refine_flag_set())
                  {
                    triangulation.levels.push_back (new internal::Triangulation::TriaLevel<dim>);
                    break;
                  }
          }


        // first clear user flags and pointers of lines; we're going
        // to use them to flag which lines need refinement
        for (typename Triangulation<dim,spacedim>::line_iterator
             line=triangulation.begin_line(); line!=triangulation.end_line(); ++line)
          {
            line->clear_user_flag();
            line->clear_user_data();
          }
        // running over all cells and lines count the number
        // n_single_lines of lines which can be stored as single
        // lines, e.g. inner lines
        unsigned int n_single_lines=0;

        // New lines to be created: number lines which are stored in
        // pairs (the children of lines must be stored in pairs)
        unsigned int n_lines_in_pairs = 0;

        // check how much space is needed on every level we need not
        // check the highest level since either - on the highest level
        // no cells are flagged for refinement - there are, but
        // prepare_refinement added another empty level
        unsigned int needed_vertices = 0;
        for (int level=triangulation.levels.size()-2; level>=0; --level)
          {
            // count number of flagged cells on this level and compute
            // how many new vertices and new lines will be needed
            unsigned int needed_cells = 0;

            typename Triangulation<dim,spacedim>::active_cell_iterator
            cell = triangulation.begin_active(level),
            endc = triangulation.begin_active(level+1);
            for (; cell!=endc; ++cell)
              if (cell->refine_flag_set())
                {
                  if (cell->refine_flag_set()==RefinementCase<dim>::cut_xy)
                    {
                      needed_cells += 4;

                      // new vertex at center of cell is needed in any
                      // case
                      ++needed_vertices;

                      // the four inner lines can be stored as singles
                      n_single_lines += 4;
                    }
                  else // cut_x || cut_y
                    {
                      // set the flag showing that anisotropic
                      // refinement is used for at least one cell
                      triangulation.anisotropic_refinement = true;

                      needed_cells += 2;
                      // no vertex at center

                      // the inner line can be stored as single
                      n_single_lines += 1;

                    }

                  // mark all faces (lines) for refinement; checking
                  // locally whether the neighbor would also like to
                  // refine them is rather difficult for lines so we
                  // only flag them and after visiting all cells, we
                  // decide which lines need refinement;
                  for (unsigned int line_no=0; line_no<GeometryInfo<dim>::faces_per_cell;
                       ++line_no)
                    {
                      if (GeometryInfo<dim>::face_refinement_case(
                            cell->refine_flag_set(), line_no)==RefinementCase<1>::cut_x)
                        {
                          typename Triangulation<dim,spacedim>::line_iterator
                          line = cell->line(line_no);
                          if (line->has_children() == false)
                            {
                              line->set_user_flag ();
//TODO[WB]: we overwrite the user_index here because we later on need
// to find out which boundary object we have to ask to refine this
// line. we can't use the boundary_indicator field because that can
// only be used for lines at the boundary of the domain, but we also
// need a domain description for interior lines in the codim-1 case
                              if (spacedim > dim)
                                {
                                  if (line->at_boundary())
                                    // if possible honor boundary
                                    // indicator
                                    line->set_user_index(line->boundary_indicator());
                                  else
                                    // otherwise take manifold
                                    // description from the adjacent
                                    // cell
                                    line->set_user_index(cell->material_id());
                                }
                            }
                        }
                    }
                }


            // count number of used cells on the next higher level
            const unsigned int used_cells
              = std::count_if (triangulation.levels[level+1]->cells.used.begin(),
                               triangulation.levels[level+1]->cells.used.end(),
                               std::bind2nd (std::equal_to<bool>(), true));


            // reserve space for the used_cells cells already existing
            // on the next higher level as well as for the
            // needed_cells that will be created on that level
            triangulation.levels[level+1]
            ->reserve_space (used_cells+needed_cells, 2, spacedim);

            // reserve space for needed_cells new quads on the next
            // higher level
            triangulation.levels[level+1]->cells.
            reserve_space (needed_cells,0);
          }

        // now count the lines which were flagged for refinement
        for (typename Triangulation<dim,spacedim>::line_iterator
             line=triangulation.begin_line(); line!=triangulation.end_line(); ++line)
          if (line->user_flag_set())
            {
              Assert (line->has_children() == false, ExcInternalError());
              n_lines_in_pairs += 2;
              needed_vertices  += 1;
            }
        // reserve space for n_lines_in_pairs new lines.  note, that
        // we can't reserve space for the single lines here as well,
        // as all the space reserved for lines in pairs would be
        // counted as unused and we would end up with too little space
        // to store all lines. memory reservation for n_single_lines
        // can only be done AFTER we refined the lines of the current
        // cells
        triangulation.faces->lines.
        reserve_space (n_lines_in_pairs, 0);

        // add to needed vertices how many vertices are already in use
        needed_vertices += std::count_if (triangulation.vertices_used.begin(), triangulation.vertices_used.end(),
                                          std::bind2nd (std::equal_to<bool>(), true));
        // if we need more vertices: create them, if not: leave the
        // array as is, since shrinking is not really possible because
        // some of the vertices at the end may be in use
        if (needed_vertices > triangulation.vertices.size())
          {
            triangulation.vertices.resize (needed_vertices, Point<spacedim>());
            triangulation.vertices_used.resize (needed_vertices, false);
          }


        // Do REFINEMENT on every level; exclude highest level as
        // above

        //  index of next unused vertex
        unsigned int next_unused_vertex = 0;

        // first the refinement of lines.  children are stored
        // pairwise
        if (true)
          {
            // only active objects can be refined further
            typename Triangulation<dim,spacedim>::active_line_iterator
            line = triangulation.begin_active_line(),
            endl = triangulation.end_line();
            typename Triangulation<dim,spacedim>::raw_line_iterator
            next_unused_line = triangulation.begin_raw_line ();

            for (; line!=endl; ++line)
              if (line->user_flag_set())
                {
                  // this line needs to be refined

                  // find the next unused vertex and set it
                  // appropriately
                  while (triangulation.vertices_used[next_unused_vertex] == true)
                    ++next_unused_vertex;
                  Assert (next_unused_vertex < triangulation.vertices.size(),
                          ExcTooFewVerticesAllocated());
                  triangulation.vertices_used[next_unused_vertex] = true;

                  if (spacedim == dim)
                    {
                      // for the case of a domain in an
                      // equal-dimensional space we only have to treat
                      // boundary lines differently; for interior
                      // lines we can compute the midpoint as the mean
                      // of the two vertices: if (line->at_boundary())
                      triangulation.vertices[next_unused_vertex]
                        = line->center(true);
                    }
                  else
                    // however, if spacedim>dim, we always have to ask
                    // the boundary object for its answer. We use the
                    // same object of the cell (which was stored in
                    // line->user_index() before) unless a manifold_id
                    // has been set on this very line.
                    if (line->manifold_id() == numbers::invalid_manifold_id)
                      triangulation.vertices[next_unused_vertex]
                        = triangulation.get_manifold(line->user_index()).get_new_point_on_line (line);
                    else
                      triangulation.vertices[next_unused_vertex]
                        = line->center(true);

                  // now that we created the right point, make up the
                  // two child lines.  To this end, find a pair of
                  // unused lines
                  bool pair_found=false;
                  for (; next_unused_line!=endl; ++next_unused_line)
                    if (!next_unused_line->used() &&
                        !(++next_unused_line)->used())
                      {
                        // go back to the first of the two unused
                        // lines
                        --next_unused_line;
                        pair_found=true;
                        break;
                      }
                  Assert (pair_found, ExcInternalError());

                  // there are now two consecutive unused lines, such
                  // that the children of a line will be consecutive.
                  // then set the child pointer of the present line
                  line->set_children (0, next_unused_line->index());

                  // set the two new lines
                  const typename Triangulation<dim,spacedim>::raw_line_iterator
                  children[2] = { next_unused_line,
                                  ++next_unused_line
                                };
                  // some tests; if any of the iterators should be
                  // invalid, then already dereferencing will fail
                  Assert (children[0]->used() == false, ExcCellShouldBeUnused());
                  Assert (children[1]->used() == false, ExcCellShouldBeUnused());

                  children[0]->set (internal::Triangulation
                                    ::TriaObject<1>(line->vertex_index(0),
                                                    next_unused_vertex));
                  children[1]->set (internal::Triangulation
                                    ::TriaObject<1>(next_unused_vertex,
                                                    line->vertex_index(1)));

                  children[0]->set_used_flag();
                  children[1]->set_used_flag();
                  children[0]->clear_children();
                  children[1]->clear_children();
                  children[0]->clear_user_data();
                  children[1]->clear_user_data();
                  children[0]->clear_user_flag();
                  children[1]->clear_user_flag();

                  children[0]->set_boundary_indicator (line->boundary_indicator());
                  children[1]->set_boundary_indicator (line->boundary_indicator());

                  children[0]->set_manifold_id (line->manifold_id());
                  children[1]->set_manifold_id (line->manifold_id());

                  // finally clear flag indicating the need for
                  // refinement
                  line->clear_user_flag ();
                }
          }


        // Now set up the new cells

        // reserve space for inner lines (can be stored as single
        // lines)
        triangulation.faces->lines.
        reserve_space (0,n_single_lines);

        typename Triangulation<2,spacedim>::DistortedCellList
        cells_with_distorted_children;

        // reset next_unused_line, as now also single empty places in
        // the vector can be used
        typename Triangulation<dim,spacedim>::raw_line_iterator
        next_unused_line = triangulation.begin_raw_line ();

        for (int level=0; level<static_cast<int>(triangulation.levels.size())-1; ++level)
          {

            // Remember: as we don't operate on the finest level,
            // begin_*(level+1) is allowed
            typename Triangulation<dim,spacedim>::active_cell_iterator
            cell = triangulation.begin_active(level),
            endc = triangulation.begin_active(level+1);

            typename Triangulation<dim,spacedim>::raw_cell_iterator
            next_unused_cell = triangulation.begin_raw (level+1);

            for (; cell!=endc; ++cell)
              if (cell->refine_flag_set())
                {
                  // set the user flag to indicate, that at least one
                  // line is at the boundary

                  // TODO[Tobias Leicht] find a better place to set
                  // this flag, so that we do not need so much time to
                  // check each cell here
                  if (cell->at_boundary())
                    cell->set_user_flag();

                  // actually set up the children and update neighbor
                  // information
                  create_children (triangulation,
                                   next_unused_vertex,
                                   next_unused_line,
                                   next_unused_cell,
                                   cell);

                  if ((check_for_distorted_cells == true)
                      &&
                      has_distorted_children (cell,
                                              internal::int2type<dim>(),
                                              internal::int2type<spacedim>()))
                    cells_with_distorted_children.distorted_cells.push_back (cell);
                }
          }

        return cells_with_distorted_children;
      }


      /**
       * A function that performs the refinement of a triangulation in
       * 3d.
       */
      template <int spacedim>
      static
      typename Triangulation<3,spacedim>::DistortedCellList
      execute_refinement (Triangulation<3,spacedim> &triangulation,
                          const bool check_for_distorted_cells)
      {
        const unsigned int dim = 3;

        // this function probably also works for spacedim>3 but it
        // isn't tested. it will probably be necessary to pull new
        // vertices onto the manifold just as we do for the other
        // functions above.
        Assert (spacedim == 3, ExcNotImplemented());

        // check whether a new level is needed we have to check for
        // this on the highest level only (on this, all used cells are
        // also active, so we only have to check for this)
        if (true)
          {
            typename Triangulation<dim,spacedim>::raw_cell_iterator
            cell = triangulation.begin_active (triangulation.levels.size()-1),
            endc = triangulation.end();
            for (; cell != endc; ++cell)
              if (cell->used())
                if (cell->refine_flag_set())
                  {
                    triangulation.levels.push_back (new internal::Triangulation::TriaLevel<dim>);
                    break;
                  }
          }


        // first clear user flags for quads and lines; we're going to
        // use them to flag which lines and quads need refinement
        triangulation.faces->quads.clear_user_data();

        for (typename Triangulation<dim,spacedim>::line_iterator
             line=triangulation.begin_line(); line!=triangulation.end_line(); ++line)
          line->clear_user_flag();
        for (typename Triangulation<dim,spacedim>::quad_iterator
             quad=triangulation.begin_quad(); quad!=triangulation.end_quad(); ++quad)
          quad->clear_user_flag();

        // create an array of face refine cases. User indices of faces
        // will be set to values corresponding with indices in this
        // array.
        const RefinementCase<dim-1>  face_refinement_cases[4]=
        {
          RefinementCase<dim-1>::no_refinement,
          RefinementCase<dim-1>::cut_x,
          RefinementCase<dim-1>::cut_y,
          RefinementCase<dim-1>::cut_xy
        };

        // check how much space is needed on every level we need not
        // check the highest level since either
        // - on the highest level no cells are flagged for refinement
        // - there are, but prepare_refinement added another empty
        // level which then is the highest level

        // variables to hold the number of newly to be created
        // vertices, lines and quads. as these are stored globally,
        // declare them outside the loop over al levels. we need lines
        // and quads in pairs for refinement of old ones and lines and
        // quads, that can be stored as single ones, as they are newly
        // created in the inside of an existing cell
        unsigned int needed_vertices = 0;
        unsigned int needed_lines_single  = 0;
        unsigned int needed_quads_single  = 0;
        unsigned int needed_lines_pair  = 0;
        unsigned int needed_quads_pair  = 0;
        for (int level=triangulation.levels.size()-2; level>=0; --level)
          {
            // count number of flagged cells on this level and compute
            // how many new vertices and new lines will be needed
            unsigned int new_cells = 0;

            typename Triangulation<dim,spacedim>::active_cell_iterator
            acell = triangulation.begin_active(level),
            aendc = triangulation.begin_active(level+1);
            for (; acell!=aendc; ++acell)
              if (acell->refine_flag_set())
                {
                  RefinementCase<dim> ref_case=acell->refine_flag_set();

                  // now for interior vertices, lines and quads, which
                  // are needed in any case
                  if (ref_case==RefinementCase<dim>::cut_x ||
                      ref_case==RefinementCase<dim>::cut_y ||
                      ref_case==RefinementCase<dim>::cut_z)
                    {
                      ++needed_quads_single;
                      new_cells+=2;
                      triangulation.anisotropic_refinement=true;
                    }
                  else if (ref_case==RefinementCase<dim>::cut_xy ||
                           ref_case==RefinementCase<dim>::cut_xz ||
                           ref_case==RefinementCase<dim>::cut_yz)
                    {
                      ++needed_lines_single;
                      needed_quads_single += 4;
                      new_cells+=4;
                      triangulation.anisotropic_refinement=true;
                    }
                  else if  (ref_case==RefinementCase<dim>::cut_xyz)
                    {
                      ++needed_vertices;
                      needed_lines_single += 6;
                      needed_quads_single += 12;
                      new_cells+=8;
                    }
                  else
                    {
                      // we should never get here
                      Assert(false, ExcInternalError());
                    }

                  // mark all faces for refinement; checking locally
                  // if and how the neighbor would like to refine
                  // these is difficult so we only flag them and after
                  // visiting all cells, we decide which faces need
                  // which refinement;
                  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell;
                       ++face)
                    {
                      typename Triangulation<dim,spacedim>::face_iterator
                      aface = acell->face(face);
                      // get the RefineCase this faces has for the
                      // given RefineCase of the cell
                      RefinementCase<dim-1> face_ref_case=
                        GeometryInfo<dim>::face_refinement_case(ref_case,
                                                                face,
                                                                acell->face_orientation(face),
                                                                acell->face_flip(face),
                                                                acell->face_rotation(face));
                      // only do something, if this face has to be
                      // refined
                      if (face_ref_case)
                        {
                          if (face_ref_case==RefinementCase<dim-1>::isotropic_refinement)
                            {
                              if (aface->number_of_children()<4)
                                // we use user_flags to denote needed
                                // isotropic refinement
                                aface->set_user_flag();
                            }
                          else if (aface->refinement_case()!=face_ref_case)
                            // we use user_indices to denote needed
                            // anisotropic refinement. note, that we
                            // can have at most one anisotropic
                            // refinement case for this face, as
                            // otherwise prepare_refinement() would
                            // have changed one of the cells to yield
                            // isotropic refinement at this
                            // face. therefore we set the user_index
                            // uniquely
                            {
                              Assert(aface->refinement_case()==RefinementCase<dim-1>::isotropic_refinement ||
                                     aface->refinement_case()==RefinementCase<dim-1>::no_refinement,
                                     ExcInternalError());
                              aface->set_user_index(face_ref_case);
                            }
                        }
                    }// for all faces

                  // flag all lines, that have to be refined
                  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                    if (GeometryInfo<dim>::line_refinement_case(ref_case,line) &&
                        !acell->line(line)->has_children())
                      acell->line(line)->set_user_flag();

                }// if refine_flag set and for all cells on this level


            // count number of used cells on the next higher level
            const unsigned int used_cells
              = std::count_if (triangulation.levels[level+1]->cells.used.begin(),
                               triangulation.levels[level+1]->cells.used.end(),
                               std::bind2nd (std::equal_to<bool>(), true));


            // reserve space for the used_cells cells already existing
            // on the next higher level as well as for the
            // 8*flagged_cells that will be created on that level
            triangulation.levels[level+1]
            ->reserve_space (used_cells+new_cells, 3, spacedim);
            // reserve space for 8*flagged_cells new hexes on the next
            // higher level
            triangulation.levels[level+1]->cells.reserve_space (new_cells);
          }// for all levels
        // now count the quads and lines which were flagged for
        // refinement
        for (typename Triangulation<dim,spacedim>::quad_iterator
             quad=triangulation.begin_quad(); quad!=triangulation.end_quad(); ++quad)
          {
            if (quad->user_flag_set())
              {
                // isotropic refinement: 1 interior vertex, 4 quads
                // and 4 interior lines. we store the interior lines
                // in pairs in case the face is already or will be
                // refined anisotropically
                needed_quads_pair += 4;
                needed_lines_pair += 4;
                needed_vertices += 1;
              }
            if (quad->user_index())
              {
                // anisotropic refinement: 1 interior
                // line and two quads
                needed_quads_pair += 2;
                needed_lines_single += 1;
                // there is a kind of complicated situation here which
                // requires our attention. if the quad is refined
                // isotropcally, two of the interior lines will get a
                // new mother line - the interior line of our
                // anisotropically refined quad. if those two lines
                // are not consecutive, we cannot do so and have to
                // replace them by two lines that are consecutive. we
                // try to avoid that situation, but it may happen
                // nevertheless throug repeated refinement and
                // coarsening. thus we have to check here, as we will
                // need some additional space to store those new lines
                // in case we need them...
                if (quad->has_children())
                  {
                    Assert(quad->refinement_case()==RefinementCase<dim-1>::isotropic_refinement, ExcInternalError());
                    if ((face_refinement_cases[quad->user_index()]==RefinementCase<dim-1>::cut_x
                         && (quad->child(0)->line_index(1)+1!=quad->child(2)->line_index(1))) ||
                        (face_refinement_cases[quad->user_index()]==RefinementCase<dim-1>::cut_y
                         && (quad->child(0)->line_index(3)+1!=quad->child(1)->line_index(3))))
                      needed_lines_pair +=2;
                  }
              }
          }

        for (typename Triangulation<dim,spacedim>::line_iterator
             line=triangulation.begin_line(); line!=triangulation.end_line(); ++line)
          if (line->user_flag_set())
            {
              needed_lines_pair += 2;
              needed_vertices += 1;
            }

        // reserve space for needed_lines new lines stored in pairs
        triangulation.faces->lines.
        reserve_space (needed_lines_pair,needed_lines_single);
        // reserve space for needed_quads new quads stored in pairs
        triangulation.faces->quads.
        reserve_space (needed_quads_pair,needed_quads_single);


        // add to needed vertices how many vertices are already in use
        needed_vertices += std::count_if (triangulation.vertices_used.begin(), triangulation.vertices_used.end(),
                                          std::bind2nd (std::equal_to<bool>(), true));
        // if we need more vertices: create them, if not: leave the
        // array as is, since shrinking is not really possible because
        // some of the vertices at the end may be in use
        if (needed_vertices > triangulation.vertices.size())
          {
            triangulation.vertices.resize (needed_vertices, Point<spacedim>());
            triangulation.vertices_used.resize (needed_vertices, false);
          }


        ///////////////////////////////////////////
        // Before we start with the actual refinement, we do some
        // sanity checks if in debug mode. especially, we try to catch
        // the notorious problem with lines being twice refined,
        // i.e. there are cells adjacent at one line ("around the
        // edge", but not at a face), with two cells differing by more
        // than one refinement level
        //
        // this check is very simple to implement here, since we have
        // all lines flagged if they shall be refined
#ifdef DEBUG
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
             cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell)
          if (!cell->refine_flag_set())
            for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
              if (cell->line(line)->has_children())
                for (unsigned int c=0; c<2; ++c)
                  Assert (cell->line(line)->child(c)->user_flag_set() == false,
                          ExcInternalError());
#endif

        ///////////////////////////////////////////
        // Do refinement on every level
        //
        // To make life a bit easier, we first refine those lines and
        // quads that were flagged for refinement and then compose the
        // newly to be created cells.
        //
        // index of next unused vertex
        unsigned int next_unused_vertex = 0;

        // first for lines
        if (true)
          {
            // only active objects can be refined further
            typename Triangulation<dim,spacedim>::active_line_iterator
            line = triangulation.begin_active_line(),
            endl = triangulation.end_line();
            typename Triangulation<dim,spacedim>::raw_line_iterator
            next_unused_line = triangulation.begin_raw_line ();

            for (; line!=endl; ++line)
              if (line->user_flag_set())
                {
                  // this line needs to be refined

                  // find the next unused vertex and set it
                  // appropriately
                  while (triangulation.vertices_used[next_unused_vertex] == true)
                    ++next_unused_vertex;
                  Assert (next_unused_vertex < triangulation.vertices.size(),
                          ExcTooFewVerticesAllocated());
                  triangulation.vertices_used[next_unused_vertex] = true;

                  triangulation.vertices[next_unused_vertex]
                    = line->center(true);

                  // now that we created the right point, make up the
                  // two child lines (++ takes care of the end of the
                  // vector)
                  next_unused_line=triangulation.faces->lines.next_free_pair_object(triangulation);
                  Assert(next_unused_line.state() == IteratorState::valid,
                         ExcInternalError());

                  // now we found two consecutive unused lines, such
                  // that the children of a line will be consecutive.
                  // then set the child pointer of the present line
                  line->set_children (0, next_unused_line->index());

                  // set the two new lines
                  const typename Triangulation<dim,spacedim>::raw_line_iterator
                  children[2] = { next_unused_line,
                                  ++next_unused_line
                                };

                  // some tests; if any of the iterators should be
                  // invalid, then already dereferencing will fail
                  Assert (children[0]->used() == false, ExcCellShouldBeUnused());
                  Assert (children[1]->used() == false, ExcCellShouldBeUnused());

                  children[0]->set (internal::Triangulation
                                    ::TriaObject<1>(line->vertex_index(0),
                                                    next_unused_vertex));
                  children[1]->set (internal::Triangulation
                                    ::TriaObject<1>(next_unused_vertex,
                                                    line->vertex_index(1)));

                  children[0]->set_used_flag();
                  children[1]->set_used_flag();
                  children[0]->clear_children();
                  children[1]->clear_children();
                  children[0]->clear_user_data();
                  children[1]->clear_user_data();
                  children[0]->clear_user_flag();
                  children[1]->clear_user_flag();

                  children[0]->set_boundary_indicator (line->boundary_indicator());
                  children[1]->set_boundary_indicator (line->boundary_indicator());

                  children[0]->set_manifold_id (line->manifold_id());
                  children[1]->set_manifold_id (line->manifold_id());

                  // finally clear flag
                  // indicating the need
                  // for refinement
                  line->clear_user_flag ();
                }
          }


        ///////////////////////////////////////
        // now refine marked quads
        ///////////////////////////////////////

        // here we encounter several cases:

        // a) the quad is unrefined and shall be refined isotropically

        // b) the quad is unrefined and shall be refined
        // anisotropically

        // c) the quad is unrefined and shall be refined both
        // anisotropically and isotropically (this is reduced to case
        // b) and then case b) for the children again)

        // d) the quad is refined anisotropically and shall be refined
        // isotropically (this is reduced to case b) for the
        // anisotropic children)

        // e) the quad is refined isotropically and shall be refined
        // anisotropically (this is transformed to case c), however we
        // might have to renumber/rename children...)

        // we need a loop in cases c) and d), as the anisotropic
        // children migt have a lower index than the mother quad
        for (unsigned int loop=0; loop<2; ++loop)
          {
            // usually, only active objects can be refined
            // further. however, in cases d) and e) that is not true,
            // so we have to use 'normal' iterators here
            typename Triangulation<dim,spacedim>::quad_iterator
            quad = triangulation.begin_quad(),
            endq = triangulation.end_quad();
            typename Triangulation<dim,spacedim>::raw_line_iterator
            next_unused_line = triangulation.begin_raw_line ();
            typename Triangulation<dim,spacedim>::raw_quad_iterator
            next_unused_quad = triangulation.begin_raw_quad ();

            for (; quad!=endq; ++quad)
              {
                if (quad->user_index())
                  {
                    RefinementCase<dim-1> aniso_quad_ref_case=face_refinement_cases[quad->user_index()];
                    // there is one unlikely event here, where we
                    // already have refind the face: if the face was
                    // refined anisotropically and we want to refine
                    // it isotropically, both children are flagged for
                    // anisotropic refinement. however, if those
                    // children were already flagged for anisotropic
                    // refinement, they might already be processed and
                    // refined.
                    if (aniso_quad_ref_case == quad->refinement_case())
                      continue;

                    Assert(quad->refinement_case()==RefinementCase<dim-1>::cut_xy ||
                           quad->refinement_case()==RefinementCase<dim-1>::no_refinement,
                           ExcInternalError());

                    // this quad needs to be refined anisotropically
                    Assert(quad->user_index() == RefinementCase<dim-1>::cut_x ||
                           quad->user_index() == RefinementCase<dim-1>::cut_y,
                           ExcInternalError());

                    // make the new line interior to the quad
                    typename Triangulation<dim,spacedim>::raw_line_iterator new_line;

                    new_line=triangulation.faces->lines.next_free_single_object(triangulation);
                    Assert (new_line->used() == false,
                            ExcCellShouldBeUnused());

                    // first collect the
                    // indices of the vertices:
                    // *--1--*
                    // |  |  |
                    // |  |  |    cut_x
                    // |  |  |
                    // *--0--*
                    //
                    // *-----*
                    // |     |
                    // 0-----1    cut_y
                    // |     |
                    // *-----*
                    unsigned int vertex_indices[2];
                    if (aniso_quad_ref_case==RefinementCase<dim-1>::cut_x)
                      {
                        vertex_indices[0]=quad->line(2)->child(0)->vertex_index(1);
                        vertex_indices[1]=quad->line(3)->child(0)->vertex_index(1);
                      }
                    else
                      {
                        vertex_indices[0]=quad->line(0)->child(0)->vertex_index(1);
                        vertex_indices[1]=quad->line(1)->child(0)->vertex_index(1);
                      }

                    new_line->set (internal::Triangulation::
                                   TriaObject<1>(vertex_indices[0], vertex_indices[1]));
                    new_line->set_used_flag();
                    new_line->clear_user_flag();
                    new_line->clear_user_data();
                    new_line->clear_children();
                    new_line->set_boundary_indicator(quad->boundary_indicator());
                    new_line->set_manifold_id(quad->manifold_id());

                    // child 0 and 1 of a line are switched if the
                    // line orientation is false. set up a miniature
                    // table, indicating which child to take for line
                    // orientations false and true. first index: child
                    // index in standard orientation, second index:
                    // line orientation
                    const unsigned int index[2][2]=
                    {
                      {1,0},   // child 0, line_orientation=false and true
                      {0,1}
                    };  // child 1, line_orientation=false and true

                    // find some space (consecutive) for the two newly
                    // to be created quads.
                    typename Triangulation<dim,spacedim>::raw_quad_iterator new_quads[2];

                    next_unused_quad=triangulation.faces->quads.next_free_pair_object(triangulation);
                    new_quads[0] = next_unused_quad;
                    Assert (new_quads[0]->used() == false, ExcCellShouldBeUnused());

                    ++next_unused_quad;
                    new_quads[1] = next_unused_quad;
                    Assert (new_quads[1]->used() == false, ExcCellShouldBeUnused());


                    if (aniso_quad_ref_case==RefinementCase<dim-1>::cut_x)
                      {
                        new_quads[0]->set (internal::Triangulation
                                           ::TriaObject<2>(quad->line_index(0),
                                                           new_line->index(),
                                                           quad->line(2)->child(index[0][quad->line_orientation(2)])->index(),
                                                           quad->line(3)->child(index[0][quad->line_orientation(3)])->index()));
                        new_quads[1]->set (internal::Triangulation
                                           ::TriaObject<2>(new_line->index(),
                                                           quad->line_index(1),
                                                           quad->line(2)->child(index[1][quad->line_orientation(2)])->index(),
                                                           quad->line(3)->child(index[1][quad->line_orientation(3)])->index()));
                      }
                    else
                      {
                        new_quads[0]->set (internal::Triangulation
                                           ::TriaObject<2>(quad->line(0)->child(index[0][quad->line_orientation(0)])->index(),
                                                           quad->line(1)->child(index[0][quad->line_orientation(1)])->index(),
                                                           quad->line_index(2),
                                                           new_line->index()));
                        new_quads[1]->set (internal::Triangulation
                                           ::TriaObject<2>(quad->line(0)->child(index[1][quad->line_orientation(0)])->index(),
                                                           quad->line(1)->child(index[1][quad->line_orientation(1)])->index(),
                                                           new_line->index(),
                                                           quad->line_index(3)));
                      }

                    for (unsigned int i=0; i<2; ++i)
                      {
                        new_quads[i]->set_used_flag();
                        new_quads[i]->clear_user_flag();
                        new_quads[i]->clear_user_data();
                        new_quads[i]->clear_children();
                        new_quads[i]->set_boundary_indicator (quad->boundary_indicator());
                        new_quads[i]->set_manifold_id (quad->manifold_id());
                        // set all line orientations to true, change
                        // this after the loop, as we have to consider
                        // different lines for each child
                        for (unsigned int j=0; j<GeometryInfo<dim>::lines_per_face; ++j)
                          new_quads[i]->set_line_orientation(j,true);
                      }
                    // now set the line orientation of children of
                    // outer lines correctly, the lines in the
                    // interior of the refined quad are automatically
                    // oriented conforming to the standard
                    new_quads[0]->set_line_orientation(0,quad->line_orientation(0));
                    new_quads[0]->set_line_orientation(2,quad->line_orientation(2));
                    new_quads[1]->set_line_orientation(1,quad->line_orientation(1));
                    new_quads[1]->set_line_orientation(3,quad->line_orientation(3));
                    if (aniso_quad_ref_case==RefinementCase<dim-1>::cut_x)
                      {
                        new_quads[0]->set_line_orientation(3,quad->line_orientation(3));
                        new_quads[1]->set_line_orientation(2,quad->line_orientation(2));
                      }
                    else
                      {
                        new_quads[0]->set_line_orientation(1,quad->line_orientation(1));
                        new_quads[1]->set_line_orientation(0,quad->line_orientation(0));
                      }

                    // test, whether this face is refined
                    // isotropically already. if so, set the correct
                    // children pointers.
                    if (quad->refinement_case()==RefinementCase<dim-1>::cut_xy)
                      {
                        // we will put a new refinemnt level of
                        // anisotropic refinement between the
                        // unrefined and isotropically refined quad
                        // ending up with the same fine quads but
                        // introducing anisotropically refined ones as
                        // children of the unrefined quad and mother
                        // cells of the original fine ones.

                        // this process includes the creation of a new
                        // middle line which we will assign as the
                        // mother line of two of the existing inner
                        // lines. If those inner lines are not
                        // consecutive in memory, we won't find them
                        // later on, so we have to create new ones
                        // instead and replace all occurrences of the
                        // old ones with those new ones. As this is
                        // kind of ugly, we hope we don't have to do
                        // it often...
                        typename Triangulation<dim,spacedim>::line_iterator old_child[2];
                        if (aniso_quad_ref_case==RefinementCase<dim-1>::cut_x)
                          {
                            old_child[0]=quad->child(0)->line(1);
                            old_child[1]=quad->child(2)->line(1);
                          }
                        else
                          {
                            Assert(aniso_quad_ref_case==RefinementCase<dim-1>::cut_y, ExcInternalError());

                            old_child[0]=quad->child(0)->line(3);
                            old_child[1]=quad->child(1)->line(3);
                          }

                        if (old_child[0]->index()+1 != old_child[1]->index())
                          {
                            // this is exactly the ugly case we taked
                            // about. so, no coimplaining, lets get
                            // two new lines and copy all info
                            typename Triangulation<dim,spacedim>::raw_line_iterator new_child[2];

                            new_child[0]=new_child[1]=triangulation.faces->lines.next_free_pair_object(triangulation);
                            ++new_child[1];

                            new_child[0]->set_used_flag();
                            new_child[1]->set_used_flag();

                            const int old_index_0=old_child[0]->index(),
                                      old_index_1=old_child[1]->index(),
                                      new_index_0=new_child[0]->index(),
                                      new_index_1=new_child[1]->index();

                            // loop over all quads and replace the old
                            // lines
                            for (unsigned int q=0; q<triangulation.faces->quads.cells.size(); ++q)
                              for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
                                {
                                  const int this_index=triangulation.faces->quads.cells[q].face(l);
                                  if (this_index==old_index_0)
                                    triangulation.faces->quads.cells[q].set_face(l,new_index_0);
                                  else if (this_index==old_index_1)
                                    triangulation.faces->quads.cells[q].set_face(l,new_index_1);
                                }
                            // now we have to copy all information of
                            // the two lines
                            for (unsigned int i=0; i<2; ++i)
                              {
                                Assert(!old_child[i]->has_children(), ExcInternalError());

                                new_child[i]->set(internal::Triangulation::TriaObject<1>(old_child[i]->vertex_index(0),
                                                                                         old_child[i]->vertex_index(1)));
                                new_child[i]->set_boundary_indicator(old_child[i]->boundary_indicator());
                                new_child[i]->set_manifold_id(old_child[i]->manifold_id());
                                new_child[i]->set_user_index(old_child[i]->user_index());
                                if (old_child[i]->user_flag_set())
                                  new_child[i]->set_user_flag();
                                else
                                  new_child[i]->clear_user_flag();

                                new_child[i]->clear_children();

                                old_child[i]->clear_user_flag();
                                old_child[i]->clear_user_index();
                                old_child[i]->clear_used_flag();
                              }
                          }
                        // now that we cared about the lines, go on
                        // with the quads themselves, where we might
                        // encounter similar situations...
                        if (aniso_quad_ref_case==RefinementCase<dim-1>::cut_x)
                          {
                            new_line->set_children(0, quad->child(0)->line_index(1));
                            Assert(new_line->child(1)==quad->child(2)->line(1),
                                   ExcInternalError());
                            // now evereything is quite
                            // complicated. we have the children
                            // numbered according to
                            //
                            // *---*---*
                            // |n+2|n+3|
                            // *---*---*
                            // | n |n+1|
                            // *---*---*
                            //
                            // from the original isotropic
                            // refinement. we have to reorder them as
                            //
                            // *---*---*
                            // |n+1|n+3|
                            // *---*---*
                            // | n |n+2|
                            // *---*---*
                            //
                            // such that n and n+1 are consecutive
                            // children of m and n+2 and n+3 are
                            // consecutive children of m+1, where m
                            // and m+1 are given as in
                            //
                            // *---*---*
                            // |   |   |
                            // | m |m+1|
                            // |   |   |
                            // *---*---*
                            //
                            // this is a bit ugly, of course: loop
                            // over all cells on all levels and look
                            // for faces n+1 (switch_1) and n+2
                            // (switch_2).
                            const typename Triangulation<dim,spacedim>::quad_iterator
                            switch_1=quad->child(1),
                            switch_2=quad->child(2);
                            const int switch_1_index=switch_1->index();
                            const int switch_2_index=switch_2->index();
                            for (unsigned int l=0; l<triangulation.levels.size(); ++l)
                              for (unsigned int h=0; h<triangulation.levels[l]->cells.cells.size(); ++h)
                                for (unsigned int q=0; q<GeometryInfo<dim>::faces_per_cell; ++q)
                                  {
                                    const int face_index=triangulation.levels[l]->cells.cells[h].face(q);
                                    if (face_index==switch_1_index)
                                      triangulation.levels[l]->cells.cells[h].set_face(q,switch_2_index);
                                    else if (face_index==switch_2_index)
                                      triangulation.levels[l]->cells.cells[h].set_face(q,switch_1_index);
                                  }
                            // now we have to copy all information of
                            // the two quads
                            const unsigned int switch_1_lines[4]=
                            {
                              switch_1->line_index(0),
                              switch_1->line_index(1),
                              switch_1->line_index(2),
                              switch_1->line_index(3)
                            };
                            const bool switch_1_line_orientations[4]=
                            {
                              switch_1->line_orientation(0),
                              switch_1->line_orientation(1),
                              switch_1->line_orientation(2),
                              switch_1->line_orientation(3)
                            };
                            const types::boundary_id switch_1_boundary_indicator=switch_1->boundary_indicator();
                            const unsigned int switch_1_user_index=switch_1->user_index();
                            const bool switch_1_user_flag=switch_1->user_flag_set();
                            const RefinementCase<dim-1> switch_1_refinement_case=switch_1->refinement_case();
                            const int switch_1_first_child_pair=(switch_1_refinement_case ? switch_1->child_index(0) : -1);
                            const int switch_1_second_child_pair=(switch_1_refinement_case==RefinementCase<dim-1>::cut_xy ? switch_1->child_index(2) : -1);

                            switch_1->set(internal::Triangulation::TriaObject<2>(switch_2->line_index(0),
                                                                                 switch_2->line_index(1),
                                                                                 switch_2->line_index(2),
                                                                                 switch_2->line_index(3)));
                            switch_1->set_line_orientation(0, switch_2->line_orientation(0));
                            switch_1->set_line_orientation(1, switch_2->line_orientation(1));
                            switch_1->set_line_orientation(2, switch_2->line_orientation(2));
                            switch_1->set_line_orientation(3, switch_2->line_orientation(3));
                            switch_1->set_boundary_indicator(switch_2->boundary_indicator());
                            switch_1->set_manifold_id(switch_2->manifold_id());
                            switch_1->set_user_index(switch_2->user_index());
                            if (switch_2->user_flag_set())
                              switch_1->set_user_flag();
                            else
                              switch_1->clear_user_flag();
                            switch_1->clear_refinement_case();
                            switch_1->set_refinement_case(switch_2->refinement_case());
                            switch_1->clear_children();
                            if (switch_2->refinement_case())
                              switch_1->set_children(0, switch_2->child_index(0));
                            if (switch_2->refinement_case()==RefinementCase<dim-1>::cut_xy)
                              switch_1->set_children(2, switch_2->child_index(2));

                            switch_2->set(internal::Triangulation::TriaObject<2>(switch_1_lines[0],
                                                                                 switch_1_lines[1],
                                                                                 switch_1_lines[2],
                                                                                 switch_1_lines[3]));
                            switch_2->set_line_orientation(0, switch_1_line_orientations[0]);
                            switch_2->set_line_orientation(1, switch_1_line_orientations[1]);
                            switch_2->set_line_orientation(2, switch_1_line_orientations[2]);
                            switch_2->set_line_orientation(3, switch_1_line_orientations[3]);
                            switch_2->set_boundary_indicator(switch_1_boundary_indicator);
                            switch_2->set_manifold_id(switch_1->manifold_id());
                            switch_2->set_user_index(switch_1_user_index);
                            if (switch_1_user_flag)
                              switch_2->set_user_flag();
                            else
                              switch_2->clear_user_flag();
                            switch_2->clear_refinement_case();
                            switch_2->set_refinement_case(switch_1_refinement_case);
                            switch_2->clear_children();
                            switch_2->set_children(0, switch_1_first_child_pair);
                            switch_2->set_children(2, switch_1_second_child_pair);

                            new_quads[0]->set_refinement_case(RefinementCase<2>::cut_y);
                            new_quads[0]->set_children(0, quad->child_index(0));
                            new_quads[1]->set_refinement_case(RefinementCase<2>::cut_y);
                            new_quads[1]->set_children(0, quad->child_index(2));
                          }
                        else
                          {
                            new_quads[0]->set_refinement_case(RefinementCase<2>::cut_x);
                            new_quads[0]->set_children(0, quad->child_index(0));
                            new_quads[1]->set_refinement_case(RefinementCase<2>::cut_x);
                            new_quads[1]->set_children(0, quad->child_index(2));
                            new_line->set_children(0, quad->child(0)->line_index(3));
                            Assert(new_line->child(1)==quad->child(1)->line(3),
                                   ExcInternalError());
                          }
                        quad->clear_children();
                      }

                    // note these quads as children to the present one
                    quad->set_children (0, new_quads[0]->index());

                    quad->set_refinement_case(aniso_quad_ref_case);

                    // finally clear flag indicating the need for
                    // refinement
                    quad->clear_user_data ();
                  } // if (anisotropic refinement)

                if (quad->user_flag_set())
                  {
                    // this quad needs to be refined isotropically

                    // first of all: we only get here in the first run
                    // of the loop
                    Assert(loop==0,ExcInternalError());

                    // find the next unused vertex. we'll need this in
                    // any case
                    while (triangulation.vertices_used[next_unused_vertex] == true)
                      ++next_unused_vertex;
                    Assert (next_unused_vertex < triangulation.vertices.size(),
                            ExcTooFewVerticesAllocated());

                    // now: if the quad is refined anisotropically
                    // already, set the anisotropic refinement flag
                    // for both children. Additionally, we have to
                    // refine the inner line, as it is an outer line
                    // of the two (anisotropic) children
                    const RefinementCase<dim-1> quad_ref_case=quad->refinement_case();

                    if (quad_ref_case==RefinementCase<dim-1>::cut_x ||
                        quad_ref_case==RefinementCase<dim-1>::cut_y)
                      {
                        // set the 'opposite' refine case for children
                        quad->child(0)->set_user_index(RefinementCase<dim-1>::cut_xy-quad_ref_case);
                        quad->child(1)->set_user_index(RefinementCase<dim-1>::cut_xy-quad_ref_case);
                        // refine the inner line
                        typename Triangulation<dim,spacedim>::line_iterator middle_line;
                        if (quad_ref_case==RefinementCase<dim-1>::cut_x)
                          middle_line=quad->child(0)->line(1);
                        else
                          middle_line=quad->child(0)->line(3);

                        // if the face has been refined
                        // anisotropically in the last refinement step
                        // it might be, that it is flagged already and
                        // that the middle line is thus refined
                        // already. if not create children.
                        if (!middle_line->has_children())
                          {
                            // set the middle vertex
                            // appropriately. double refinement of
                            // quads can only happen in the interior
                            // of the domain, so we need not care
                            // about boundary quads here
                            triangulation.vertices[next_unused_vertex]
                              = middle_line->center(true);
                            triangulation.vertices_used[next_unused_vertex] = true;

                            // now search a slot for the two
                            // child lines
                            next_unused_line=triangulation.faces->lines.next_free_pair_object(triangulation);

                            // set the child pointer of the present
                            // line
                            middle_line->set_children (0, next_unused_line->index());

                            // set the two new lines
                            const typename Triangulation<dim,spacedim>::raw_line_iterator
                            children[2] = { next_unused_line,
                                            ++next_unused_line
                                          };

                            // some tests; if any of the iterators
                            // should be invalid, then already
                            // dereferencing will fail
                            Assert (children[0]->used() == false, ExcCellShouldBeUnused());
                            Assert (children[1]->used() == false, ExcCellShouldBeUnused());

                            children[0]->set (internal::Triangulation::
                                              TriaObject<1>(middle_line->vertex_index(0),
                                                            next_unused_vertex));
                            children[1]->set (internal::Triangulation::
                                              TriaObject<1>(next_unused_vertex,
                                                            middle_line->vertex_index(1)));

                            children[0]->set_used_flag();
                            children[1]->set_used_flag();
                            children[0]->clear_children();
                            children[1]->clear_children();
                            children[0]->clear_user_data();
                            children[1]->clear_user_data();
                            children[0]->clear_user_flag();
                            children[1]->clear_user_flag();

                            children[0]->set_boundary_indicator (middle_line->boundary_indicator());
                            children[1]->set_boundary_indicator (middle_line->boundary_indicator());

                            children[0]->set_manifold_id (middle_line->manifold_id());
                            children[1]->set_manifold_id (middle_line->manifold_id());
                          }
                        // now remove the flag from the quad and go to
                        // the next quad, the actual refinement of the
                        // quad takes place later on in this pass of
                        // the loop or in the next one
                        quad->clear_user_flag();
                        continue;
                      } // if (several refinement cases)

                    // if we got here, we have an unrefined quad and
                    // have to do the usual work like in an purely
                    // isotropic refinement
                    Assert(quad_ref_case==RefinementCase<dim-1>::no_refinement, ExcInternalError());

                    // set the middle vertex
                    // appropriately
                    if (quad->at_boundary() ||
                        (quad->manifold_id() != numbers::invalid_manifold_id) )
                      triangulation.vertices[next_unused_vertex]
                        = quad->center(true);
                    else
                      {
                        // it might be that the quad itself is not at
                        // the boundary, but that one of its lines
                        // actually is. in this case, the newly
                        // created vertices at the centers of the
                        // lines are not necessarily the mean values
                        // of the adjacent vertices, so do not compute
                        // the new vertex as the mean value of the 4
                        // vertices of the face, but rather as a
                        // weighted mean value of the 8 vertices which
                        // we already have (the four old ones, and the
                        // four ones inserted as middle points for the
                        // four lines). summing up some more points is
                        // generally cheaper than first asking whether
                        // one of the lines is at the boundary
                        //
                        // note that the exact weights are chosen such
                        // as to minimize the distortion of the four
                        // new quads from the optimal shape; their
                        // derivation and values is copied over from
                        // the @p{MappingQ::set_laplace_on_vector}
                        // function
                        triangulation.vertices[next_unused_vertex] =
                          quad->center(true, true);
                      }
                    triangulation.vertices_used[next_unused_vertex] = true;
                    // now that we created the right point, make up
                    // the four lines interior to the quad (++ takes
                    // care of the end of the vector)
                    typename Triangulation<dim,spacedim>::raw_line_iterator new_lines[4];

                    for (unsigned int i=0; i<4; ++i)
                      {
                        if (i%2==0)
                          // search a free pair of lines for 0. and
                          // 2. line, so that two of them end up
                          // together, which is necessary if later on
                          // we want to refine the quad
                          // anisotropically and the two lines end up
                          // as children of new line
                          next_unused_line=triangulation.faces->lines.next_free_pair_object(triangulation);

                        new_lines[i] = next_unused_line;
                        ++next_unused_line;

                        Assert (new_lines[i]->used() == false,
                                ExcCellShouldBeUnused());
                      }

                    // set the data of the four lines.  first collect
                    // the indices of the five vertices:
                    //
                    // *--3--*
                    // |  |  |
                    // 0--4--1
                    // |  |  |
                    // *--2--*
                    //
                    // the lines are numbered as follows:
                    //
                    // *--*--*
                    // |  1  |
                    // *2-*-3*
                    // |  0  |
                    // *--*--*

                    const unsigned int vertex_indices[5]
                      = { quad->line(0)->child(0)->vertex_index(1),
                          quad->line(1)->child(0)->vertex_index(1),
                          quad->line(2)->child(0)->vertex_index(1),
                          quad->line(3)->child(0)->vertex_index(1),
                          next_unused_vertex
                        };

                    new_lines[0]->set (internal::Triangulation::
                                       TriaObject<1>(vertex_indices[2], vertex_indices[4]));
                    new_lines[1]->set (internal::Triangulation::
                                       TriaObject<1>(vertex_indices[4], vertex_indices[3]));
                    new_lines[2]->set (internal::Triangulation::
                                       TriaObject<1>(vertex_indices[0], vertex_indices[4]));
                    new_lines[3]->set (internal::Triangulation::
                                       TriaObject<1>(vertex_indices[4], vertex_indices[1]));

                    for (unsigned int i=0; i<4; ++i)
                      {
                        new_lines[i]->set_used_flag();
                        new_lines[i]->clear_user_flag();
                        new_lines[i]->clear_user_data();
                        new_lines[i]->clear_children();
                        new_lines[i]->set_boundary_indicator(quad->boundary_indicator());
                        new_lines[i]->set_manifold_id(quad->manifold_id());
                      }

                    // now for the quads. again, first collect some
                    // data about the indices of the lines, with the
                    // following numbering:
                    //
                    //   .-6-.-7-.
                    //   1   9   3
                    //   .-10.11-.
                    //   0   8   2
                    //   .-4-.-5-.

                    // child 0 and 1 of a line are switched if the
                    // line orientation is false. set up a miniature
                    // table, indicating which child to take for line
                    // orientations false and true. first index: child
                    // index in standard orientation, second index:
                    // line orientation
                    const unsigned int index[2][2]=
                    {
                      {1,0},   // child 0, line_orientation=false and true
                      {0,1}
                    };  // child 1, line_orientation=false and true

                    const int line_indices[12]
                      = { quad->line(0)->child(index[0][quad->line_orientation(0)])->index(),
                          quad->line(0)->child(index[1][quad->line_orientation(0)])->index(),
                          quad->line(1)->child(index[0][quad->line_orientation(1)])->index(),
                          quad->line(1)->child(index[1][quad->line_orientation(1)])->index(),
                          quad->line(2)->child(index[0][quad->line_orientation(2)])->index(),
                          quad->line(2)->child(index[1][quad->line_orientation(2)])->index(),
                          quad->line(3)->child(index[0][quad->line_orientation(3)])->index(),
                          quad->line(3)->child(index[1][quad->line_orientation(3)])->index(),
                          new_lines[0]->index(),
                          new_lines[1]->index(),
                          new_lines[2]->index(),
                          new_lines[3]->index()
                        };

                    // find some space (consecutive)
                    // for the first two newly to be
                    // created quads.
                    typename Triangulation<dim,spacedim>::raw_quad_iterator new_quads[4];

                    next_unused_quad=triangulation.faces->quads.next_free_pair_object(triangulation);

                    new_quads[0] = next_unused_quad;
                    Assert (new_quads[0]->used() == false, ExcCellShouldBeUnused());

                    ++next_unused_quad;
                    new_quads[1] = next_unused_quad;
                    Assert (new_quads[1]->used() == false, ExcCellShouldBeUnused());

                    next_unused_quad=triangulation.faces->quads.next_free_pair_object(triangulation);
                    new_quads[2] = next_unused_quad;
                    Assert (new_quads[2]->used() == false, ExcCellShouldBeUnused());

                    ++next_unused_quad;
                    new_quads[3] = next_unused_quad;
                    Assert (new_quads[3]->used() == false, ExcCellShouldBeUnused());

                    // note these quads as children to the present one
                    quad->set_children (0, new_quads[0]->index());
                    quad->set_children (2, new_quads[2]->index());
                    new_quads[0]->set (internal::Triangulation
                                       ::TriaObject<2>(line_indices[0],
                                                       line_indices[8],
                                                       line_indices[4],
                                                       line_indices[10]));

                    quad->set_refinement_case(RefinementCase<2>::cut_xy);

                    new_quads[0]->set (internal::Triangulation
                                       ::TriaObject<2>(line_indices[0],
                                                       line_indices[8],
                                                       line_indices[4],
                                                       line_indices[10]));
                    new_quads[1]->set (internal::Triangulation
                                       ::TriaObject<2>(line_indices[8],
                                                       line_indices[2],
                                                       line_indices[5],
                                                       line_indices[11]));
                    new_quads[2]->set (internal::Triangulation
                                       ::TriaObject<2>(line_indices[1],
                                                       line_indices[9],
                                                       line_indices[10],
                                                       line_indices[6]));
                    new_quads[3]->set (internal::Triangulation
                                       ::TriaObject<2>(line_indices[9],
                                                       line_indices[3],
                                                       line_indices[11],
                                                       line_indices[7]));
                    for (unsigned int i=0; i<4; ++i)
                      {
                        new_quads[i]->set_used_flag();
                        new_quads[i]->clear_user_flag();
                        new_quads[i]->clear_user_data();
                        new_quads[i]->clear_children();
                        new_quads[i]->set_boundary_indicator (quad->boundary_indicator());
                        new_quads[i]->set_manifold_id (quad->manifold_id());
                        // set all line orientations to true, change
                        // this after the loop, as we have to consider
                        // different lines for each child
                        for (unsigned int j=0; j<GeometryInfo<dim>::lines_per_face; ++j)
                          new_quads[i]->set_line_orientation(j,true);
                      }
                    // now set the line orientation of children of
                    // outer lines correctly, the lines in the
                    // interior of the refined quad are automatically
                    // oriented conforming to the standard
                    new_quads[0]->set_line_orientation(0,quad->line_orientation(0));
                    new_quads[0]->set_line_orientation(2,quad->line_orientation(2));
                    new_quads[1]->set_line_orientation(1,quad->line_orientation(1));
                    new_quads[1]->set_line_orientation(2,quad->line_orientation(2));
                    new_quads[2]->set_line_orientation(0,quad->line_orientation(0));
                    new_quads[2]->set_line_orientation(3,quad->line_orientation(3));
                    new_quads[3]->set_line_orientation(1,quad->line_orientation(1));
                    new_quads[3]->set_line_orientation(3,quad->line_orientation(3));

                    // finally clear flag indicating the need for
                    // refinement
                    quad->clear_user_flag ();
                  } // if (isotropic refinement)
              } // for all quads
          } // looped two times over all quads, all quads refined now

        ///////////////////////////////////
        // Now, finally, set up the new
        // cells
        ///////////////////////////////////

        typename Triangulation<3,spacedim>::DistortedCellList
        cells_with_distorted_children;

        for (unsigned int level=0; level!=triangulation.levels.size()-1; ++level)
          {
            // only active objects can be refined further; remember
            // that we won't operate on the finest level, so
            // triangulation.begin_*(level+1) is allowed
            typename Triangulation<dim,spacedim>::active_hex_iterator
            hex  = triangulation.begin_active_hex(level),
            endh = triangulation.begin_active_hex(level+1);
            typename Triangulation<dim,spacedim>::raw_hex_iterator
            next_unused_hex  = triangulation.begin_raw_hex (level+1);

            for (; hex!=endh; ++hex)
              if (hex->refine_flag_set())
                {
                  // this hex needs to be refined

                  // clear flag indicating the need for refinement. do
                  // it here already, since we can't do it anymore
                  // once the cell has children
                  const RefinementCase<dim> ref_case=hex->refine_flag_set();
                  hex->clear_refine_flag ();
                  hex->set_refinement_case(ref_case);

                  // depending on the refine case we might have to
                  // create additional vertices, lines and quads
                  // interior of the hex before the actual children
                  // can be set up.

                  // in a first step: reserve the needed space for
                  // lines, quads and hexes and initialize them
                  // correctly

                  unsigned int n_new_lines=0;
                  unsigned int n_new_quads=0;
                  unsigned int n_new_hexes=0;
                  switch (ref_case)
                    {
                    case RefinementCase<dim>::cut_x:
                    case RefinementCase<dim>::cut_y:
                    case RefinementCase<dim>::cut_z:
                      n_new_lines=0;
                      n_new_quads=1;
                      n_new_hexes=2;
                      break;
                    case RefinementCase<dim>::cut_xy:
                    case RefinementCase<dim>::cut_xz:
                    case RefinementCase<dim>::cut_yz:
                      n_new_lines=1;
                      n_new_quads=4;
                      n_new_hexes=4;
                      break;
                    case RefinementCase<dim>::cut_xyz:
                      n_new_lines=6;
                      n_new_quads=12;
                      n_new_hexes=8;
                      break;
                    default:
                      Assert(false, ExcInternalError());
                      break;
                    }

                  // find some space for the newly to be created
                  // interior lines and initialize them.
                  std::vector<typename Triangulation<dim,spacedim>::raw_line_iterator>
                  new_lines(n_new_lines);
                  for (unsigned int i=0; i<n_new_lines; ++i)
                    {
                      new_lines[i] = triangulation.faces->lines.next_free_single_object(triangulation);

                      Assert (new_lines[i]->used() == false,
                              ExcCellShouldBeUnused());
                      new_lines[i]->set_used_flag();
                      new_lines[i]->clear_user_flag();
                      new_lines[i]->clear_user_data();
                      new_lines[i]->clear_children();
                      // interior line
                      new_lines[i]->set_boundary_indicator(numbers::internal_face_boundary_id);
                      // they inherit geometry description of the hex they belong to
                      new_lines[i]->set_manifold_id(hex->manifold_id());
                    }

                  // find some space for the newly to be created
                  // interior quads and initialize them.
                  std::vector<typename Triangulation<dim,spacedim>::raw_quad_iterator>
                  new_quads(n_new_quads);
                  for (unsigned int i=0; i<n_new_quads; ++i)
                    {
                      new_quads[i] = triangulation.faces->quads.next_free_single_object(triangulation);

                      Assert (new_quads[i]->used() == false,
                              ExcCellShouldBeUnused());
                      new_quads[i]->set_used_flag();
                      new_quads[i]->clear_user_flag();
                      new_quads[i]->clear_user_data();
                      new_quads[i]->clear_children();
                      // interior quad
                      new_quads[i]->set_boundary_indicator (numbers::internal_face_boundary_id);
                      // they inherit geometry description of the hex they belong to
                      new_quads[i]->set_manifold_id (hex->manifold_id());
                      // set all line orientation flags to true by
                      // default, change this afterwards, if necessary
                      for (unsigned int j=0; j<GeometryInfo<dim>::lines_per_face; ++j)
                        new_quads[i]->set_line_orientation(j,true);
                    }

                  types::subdomain_id subdomainid = hex->subdomain_id();

                  // find some space for the newly to be created hexes
                  // and initialize them.
                  std::vector<typename Triangulation<dim,spacedim>::raw_hex_iterator>
                  new_hexes(n_new_hexes);
                  for (unsigned int i=0; i<n_new_hexes; ++i)
                    {
                      if (i%2==0)
                        next_unused_hex=triangulation.levels[level+1]->cells.next_free_hex(triangulation,level+1);
                      else
                        ++next_unused_hex;

                      new_hexes[i]=next_unused_hex;

                      Assert (new_hexes[i]->used() == false,
                              ExcCellShouldBeUnused());
                      new_hexes[i]->set_used_flag();
                      new_hexes[i]->clear_user_flag();
                      new_hexes[i]->clear_user_data();
                      new_hexes[i]->clear_children();
                      // inherit material
                      // properties
                      new_hexes[i]->set_material_id (hex->material_id());
                      new_hexes[i]->set_manifold_id (hex->manifold_id());
                      new_hexes[i]->set_subdomain_id (subdomainid);

                      if (i%2)
                        new_hexes[i]->set_parent (hex->index ());
                      // set the face_orientation flag to true for all
                      // faces initially, as this is the default value
                      // which is true for all faces interior to the
                      // hex. later on go the other way round and
                      // reset faces that are at the boundary of the
                      // mother cube
                      //
                      // the same is true for the face_flip and
                      // face_rotation flags. however, the latter two
                      // are set to false by default as this is the
                      // standard value
                      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                        {
                          new_hexes[i]->set_face_orientation(f, true);
                          new_hexes[i]->set_face_flip(f, false);
                          new_hexes[i]->set_face_rotation(f, false);
                        }
                    }
                  // note these hexes as children to the present cell
                  for (unsigned int i=0; i<n_new_hexes/2; ++i)
                    hex->set_children (2*i, new_hexes[2*i]->index());

                  // we have to take into account whether the
                  // different faces are oriented correctly or in the
                  // opposite direction, so store that up front

                  // face_orientation
                  const bool f_or[6]
                    = { hex->face_orientation (0),
                        hex->face_orientation (1),
                        hex->face_orientation (2),
                        hex->face_orientation (3),
                        hex->face_orientation (4),
                        hex->face_orientation (5)
                      };

                  // face_flip
                  const bool f_fl[6]
                    = { hex->face_flip (0),
                        hex->face_flip (1),
                        hex->face_flip (2),
                        hex->face_flip (3),
                        hex->face_flip (4),
                        hex->face_flip (5)
                      };

                  // face_rotation
                  const bool f_ro[6]
                    = { hex->face_rotation (0),
                        hex->face_rotation (1),
                        hex->face_rotation (2),
                        hex->face_rotation (3),
                        hex->face_rotation (4),
                        hex->face_rotation (5)
                      };

                  // some commonly used fields which
                  // have varying size
                  const unsigned int *vertex_indices=0;
                  const typename Triangulation<dim,spacedim>::raw_line_iterator
                  *lines=0;
                  const unsigned int *line_indices=0;
                  const bool *line_orientation=0;
                  const int *quad_indices=0;

                  // little helper table, indicating, whether the
                  // child with index 0 or with index 1 can be found
                  // at the standard origin of an anisotropically
                  // refined quads in real orientation index 1:
                  // (RefineCase - 1) index 2: face_flip

                  // index 3: face rotation
                  // note: face orientation has no influence
                  const unsigned int child_at_origin[2][2][2]=
                  {
                    { { 0, 0 },  // RefinementCase<dim>::cut_x, face_flip=false, face_rotation=false and true
                      { 1, 1 }
                    }, // RefinementCase<dim>::cut_x, face_flip=true,  face_rotation=false and true
                    { { 0, 1 },  // RefinementCase<dim>::cut_y, face_flip=false, face_rotation=false and true
                      { 1, 0 }
                    }
                  };// RefinementCase<dim>::cut_y, face_flip=true,  face_rotation=false and true

                  ///////////////////////////////////////
                  //
                  // in the following we will do the same thing for
                  // each refinement case: create a new vertex (if
                  // needed), create new interior lines (if needed),
                  // create new interior quads and afterwards build
                  // the children hexes out of these and the existing
                  // subfaces of the outer quads (which have been
                  // created above). However, even if the steps are
                  // quite similar, the actual work strongly depends
                  // on the actual refinement case. therefore, we use
                  // separate blocks of code for each of these cases,
                  // which hopefully increases the readability to some
                  // extend.

                  switch (ref_case)
                    {
                    case RefinementCase<dim>::cut_x:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_x
                      //
                      // the refined cube will look
                      // like this:
                      //
                      //        *----*----*
                      //       /    /    /|
                      //      /    /    / |
                      //     /    /    /  |
                      //    *----*----*   |
                      //    |    |    |   |
                      //    |    |    |   *
                      //    |    |    |  /
                      //    |    |    | /
                      //    |    |    |/
                      //    *----*----*
                      //
                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering:

                      // face 2: front plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   |   |
                      //   |   0   |
                      //   |   |   |
                      //   *---*---*
                      //       m0
                      // face 3: back plane
                      //   (note: x,y exchanged)
                      //       m1
                      //   *---*---*
                      //   |   |   |
                      //   |   1   |
                      //   |   |   |
                      //   *---*---*
                      // face 4: bottom plane
                      //       *---*---*
                      //      /   /   /
                      //     /   2   /
                      //    /   /   /
                      //   *---*---*
                      //       m0
                      // face 5: top plane
                      //           m1
                      //       *---*---*
                      //      /   /   /
                      //     /   3   /
                      //    /   /   /
                      //   *---*---*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_x[4]
                      =
                      {
                        hex->face(2)->child(0)
                        ->line((hex->face(2)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //0
                        hex->face(3)->child(0)
                        ->line((hex->face(3)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //1
                        hex->face(4)->child(0)
                        ->line((hex->face(4)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //2
                        hex->face(5)->child(0)
                        ->line((hex->face(5)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3)         //3
                      };

                      lines=&lines_x[0];

                      unsigned int line_indices_x[4];

                      for (unsigned int i=0; i<4; ++i)
                        line_indices_x[i]=lines[i]->index();
                      line_indices=&line_indices_x[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_x[4];

                      // the middle vertice marked as m0 above is the
                      // start vertex for lines 0 and 2 in standard
                      // orientation, whereas m1 is the end vertex of
                      // lines 1 and 3 in standard orientation
                      const unsigned int middle_vertices[2]=
                      {
                        hex->line(2)->child(0)->vertex_index(1),
                        hex->line(7)->child(0)->vertex_index(1)
                      };

                      for (unsigned int i=0; i<4; ++i)
                        if (lines[i]->vertex_index(i%2)==middle_vertices[i%2])
                          line_orientation_x[i]=true;
                        else
                          {
                            // it must be the other
                            // way round then
                            Assert(lines[i]->vertex_index((i+1)%2)==middle_vertices[i%2],
                                   ExcInternalError());
                            line_orientation_x[i]=false;
                          }

                      line_orientation=&line_orientation_x[0];

                      // set up the new quad, line numbering is as
                      // indicated above
                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[0],
                                                         line_indices[1],
                                                         line_indices[2],
                                                         line_indices[3]));

                      new_quads[0]->set_line_orientation(0,line_orientation[0]);
                      new_quads[0]->set_line_orientation(1,line_orientation[1]);
                      new_quads[0]->set_line_orientation(2,line_orientation[2]);
                      new_quads[0]->set_line_orientation(3,line_orientation[3]);

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    / | x
                      //   /  | *-------*      *---------*
                      //  *   | |       |     /         /
                      //  | 0 | |       |    /         /
                      //  |   * |       |   /         /
                      //  |  /  *-------*y *---------*x
                      //  | /
                      //  |/
                      //  *
                      //
                      // children of the faces of the old hex
                      //
                      //      *---*---*        *---*---*
                      //     /|   |   |       /   /   /|
                      //    / |   |   |      / 9 / 10/ |
                      //   /  | 5 | 6 |     /   /   /  |
                      //  *   |   |   |    *---*---*   |
                      //  | 1 *---*---*    |   |   | 2 *
                      //  |  /   /   /     |   |   |  /
                      //  | / 7 / 8 /      | 3 | 4 | /
                      //  |/   /   /       |   |   |/
                      //  *---*---*        *---*---*
                      //
                      // note that we have to take care of the
                      // orientation of faces.
                      const int quad_indices_x[11]
                      =
                      {
                        new_quads[0]->index(),     //0

                        hex->face(0)->index(),     //1

                        hex->face(1)->index(),     //2

                        hex->face(2)->child_index(  child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),  //3
                        hex->face(2)->child_index(1-child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),

                        hex->face(3)->child_index(  child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),  //5
                        hex->face(3)->child_index(1-child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),

                        hex->face(4)->child_index(  child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),  //7
                        hex->face(4)->child_index(1-child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),

                        hex->face(5)->child_index(  child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]]),  //9
                        hex->face(5)->child_index(1-child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]])

                      };
                      quad_indices=&quad_indices_x[0];

                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[1],
                                                         quad_indices[0],
                                                         quad_indices[3],
                                                         quad_indices[5],
                                                         quad_indices[7],
                                                         quad_indices[9]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[0],
                                                         quad_indices[2],
                                                         quad_indices[4],
                                                         quad_indices[6],
                                                         quad_indices[8],
                                                         quad_indices[10]));
                      break;
                    }
                    case RefinementCase<dim>::cut_y:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_y
                      //
                      // the refined cube will look like this:
                      //
                      //        *---------*
                      //       /         /|
                      //      *---------* |
                      //     /         /| |
                      //    *---------* | |
                      //    |         | | |
                      //    |         | | *
                      //    |         | |/
                      //    |         | *
                      //    |         |/
                      //    *---------*
                      //
                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering:

                      // face 0: left plane
                      //       *
                      //      /|
                      //     * |
                      //    /| |
                      //   * | |
                      //   | 0 |
                      //   | | *
                      //   | |/
                      //   | *m0
                      //   |/
                      //   *
                      // face 1: right plane
                      //       *
                      //      /|
                      //   m1* |
                      //    /| |
                      //   * | |
                      //   | 1 |
                      //   | | *
                      //   | |/
                      //   | *
                      //   |/
                      //   *
                      // face 4: bottom plane
                      //       *-------*
                      //      /       /
                      //   m0*---2---*
                      //    /       /
                      //   *-------*
                      // face 5: top plane
                      //       *-------*
                      //      /       /
                      //     *---3---*m1
                      //    /       /
                      //   *-------*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_y[4]
                      =
                      {
                        hex->face(0)->child(0)
                        ->line((hex->face(0)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //0
                        hex->face(1)->child(0)
                        ->line((hex->face(1)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //1
                        hex->face(4)->child(0)
                        ->line((hex->face(4)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //2
                        hex->face(5)->child(0)
                        ->line((hex->face(5)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3)         //3
                      };

                      lines=&lines_y[0];

                      unsigned int line_indices_y[4];

                      for (unsigned int i=0; i<4; ++i)
                        line_indices_y[i]=lines[i]->index();
                      line_indices=&line_indices_y[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_y[4];

                      // the middle vertice marked as m0 above is the
                      // start vertex for lines 0 and 2 in standard
                      // orientation, whereas m1 is the end vertex of
                      // lines 1 and 3 in standard orientation
                      const unsigned int middle_vertices[2]=
                      {
                        hex->line(0)->child(0)->vertex_index(1),
                        hex->line(5)->child(0)->vertex_index(1)
                      };

                      for (unsigned int i=0; i<4; ++i)
                        if (lines[i]->vertex_index(i%2)==middle_vertices[i%2])
                          line_orientation_y[i]=true;
                        else
                          {
                            // it must be the other way round then
                            Assert(lines[i]->vertex_index((i+1)%2)==middle_vertices[i%2],
                                   ExcInternalError());
                            line_orientation_y[i]=false;
                          }

                      line_orientation=&line_orientation_y[0];

                      // set up the new quad, line numbering is as
                      // indicated above
                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[2],
                                                         line_indices[3],
                                                         line_indices[0],
                                                         line_indices[1]));

                      new_quads[0]->set_line_orientation(0,line_orientation[2]);
                      new_quads[0]->set_line_orientation(1,line_orientation[3]);
                      new_quads[0]->set_line_orientation(2,line_orientation[0]);
                      new_quads[0]->set_line_orientation(3,line_orientation[1]);

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    / | x
                      //   /  | *-------*      *---------*
                      //  *   | |       |     /         /
                      //  |   | |   0   |    /         /
                      //  |   * |       |   /         /
                      //  |  /  *-------*y *---------*x
                      //  | /
                      //  |/
                      //  *
                      //
                      // children of the faces of the old hex
                      //
                      //      *-------*        *-------*
                      //     /|       |       /   10  /|
                      //    * |       |      *-------* |
                      //   /| |   6   |     /   9   /| |
                      //  * |2|       |    *-------* |4|
                      //  | | *-------*    |       | | *
                      //  |1|/   8   /     |       |3|/
                      //  | *-------*      |   5   | *
                      //  |/   7   /       |       |/
                      //  *-------*        *-------*
                      //
                      // note that we have to take care of the
                      // orientation of faces.
                      const int quad_indices_y[11]
                      =
                      {
                        new_quads[0]->index(),     //0

                        hex->face(0)->child_index(  child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),  //1
                        hex->face(0)->child_index(1-child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),

                        hex->face(1)->child_index(  child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),  //3
                        hex->face(1)->child_index(1-child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),

                        hex->face(2)->index(),     //5

                        hex->face(3)->index(),     //6

                        hex->face(4)->child_index(  child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),  //7
                        hex->face(4)->child_index(1-child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),

                        hex->face(5)->child_index(  child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]]),  //9
                        hex->face(5)->child_index(1-child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]])

                      };
                      quad_indices=&quad_indices_y[0];

                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[1],
                                                         quad_indices[3],
                                                         quad_indices[5],
                                                         quad_indices[0],
                                                         quad_indices[7],
                                                         quad_indices[9]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[2],
                                                         quad_indices[4],
                                                         quad_indices[0],
                                                         quad_indices[6],
                                                         quad_indices[8],
                                                         quad_indices[10]));
                      break;
                    }
                    case RefinementCase<dim>::cut_z:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_z
                      //
                      // the refined cube will look like this:
                      //
                      //        *---------*
                      //       /         /|
                      //      /         / |
                      //     /         /  *
                      //    *---------*  /|
                      //    |         | / |
                      //    |         |/  *
                      //    *---------*  /
                      //    |         | /
                      //    |         |/
                      //    *---------*
                      //
                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering:

                      // face 0: left plane
                      //       *
                      //      /|
                      //     / |
                      //    /  *
                      //   *  /|
                      //   | 0 |
                      //   |/  *
                      // m0*  /
                      //   | /
                      //   |/
                      //   *
                      // face 1: right plane
                      //       *
                      //      /|
                      //     / |
                      //    /  *m1
                      //   *  /|
                      //   | 1 |
                      //   |/  *
                      //   *  /
                      //   | /
                      //   |/
                      //   *
                      // face 2: front plane
                      //   (note: x,y exchanged)
                      //   *-------*
                      //   |       |
                      // m0*---2---*
                      //   |       |
                      //   *-------*
                      // face 3: back plane
                      //   (note: x,y exchanged)
                      //   *-------*
                      //   |       |
                      //   *---3---*m1
                      //   |       |
                      //   *-------*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_z[4]
                      =
                      {
                        hex->face(0)->child(0)
                        ->line((hex->face(0)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //0
                        hex->face(1)->child(0)
                        ->line((hex->face(1)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //1
                        hex->face(2)->child(0)
                        ->line((hex->face(2)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //2
                        hex->face(3)->child(0)
                        ->line((hex->face(3)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3)         //3
                      };

                      lines=&lines_z[0];

                      unsigned int line_indices_z[4];

                      for (unsigned int i=0; i<4; ++i)
                        line_indices_z[i]=lines[i]->index();
                      line_indices=&line_indices_z[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_z[4];

                      // the middle vertex marked as m0 above is the
                      // start vertex for lines 0 and 2 in standard
                      // orientation, whereas m1 is the end vertex of
                      // lines 1 and 3 in standard orientation
                      const unsigned int middle_vertices[2]=
                      {
                        middle_vertex_index<dim,spacedim>(hex->line(8)),
                        middle_vertex_index<dim,spacedim>(hex->line(11))
                      };

                      for (unsigned int i=0; i<4; ++i)
                        if (lines[i]->vertex_index(i%2)==middle_vertices[i%2])
                          line_orientation_z[i]=true;
                        else
                          {
                            // it must be the other way round then
                            Assert(lines[i]->vertex_index((i+1)%2)==middle_vertices[i%2],
                                   ExcInternalError());
                            line_orientation_z[i]=false;
                          }

                      line_orientation=&line_orientation_z[0];

                      // set up the new quad, line numbering is as
                      // indicated above
                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[0],
                                                         line_indices[1],
                                                         line_indices[2],
                                                         line_indices[3]));

                      new_quads[0]->set_line_orientation(0,line_orientation[0]);
                      new_quads[0]->set_line_orientation(1,line_orientation[1]);
                      new_quads[0]->set_line_orientation(2,line_orientation[2]);
                      new_quads[0]->set_line_orientation(3,line_orientation[3]);

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    / | x
                      //   /  | *-------*      *---------*
                      //  *   | |       |     /         /
                      //  |   | |       |    /    0    /
                      //  |   * |       |   /         /
                      //  |  /  *-------*y *---------*x
                      //  | /
                      //  |/
                      //  *
                      //
                      // children of the faces of the old hex
                      //
                      //      *---*---*        *-------*
                      //     /|   8   |       /       /|
                      //    / |       |      /   10  / |
                      //   /  *-------*     /       /  *
                      //  * 2/|       |    *-------* 4/|
                      //  | / |   7   |    |   6   | / |
                      //  |/1 *-------*    |       |/3 *
                      //  *  /       /     *-------*  /
                      //  | /   9   /      |       | /
                      //  |/       /       |   5   |/
                      //  *-------*        *---*---*
                      //
                      // note that we have to take care of the
                      // orientation of faces.
                      const int quad_indices_z[11]
                      =
                      {
                        new_quads[0]->index(),     //0

                        hex->face(0)->child_index(  child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),  //1
                        hex->face(0)->child_index(1-child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),

                        hex->face(1)->child_index(  child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),  //3
                        hex->face(1)->child_index(1-child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),

                        hex->face(2)->child_index(  child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),  //5
                        hex->face(2)->child_index(1-child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),

                        hex->face(3)->child_index(  child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),  //7
                        hex->face(3)->child_index(1-child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),

                        hex->face(4)->index(),     //9

                        hex->face(5)->index()      //10
                      };
                      quad_indices=&quad_indices_z[0];

                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[1],
                                                         quad_indices[3],
                                                         quad_indices[5],
                                                         quad_indices[7],
                                                         quad_indices[9],
                                                         quad_indices[0]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[2],
                                                         quad_indices[4],
                                                         quad_indices[6],
                                                         quad_indices[8],
                                                         quad_indices[0],
                                                         quad_indices[10]));
                      break;
                    }
                    case RefinementCase<dim>::cut_xy:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_xy
                      //
                      // the refined cube will look like this:
                      //
                      //        *----*----*
                      //       /    /    /|
                      //      *----*----* |
                      //     /    /    /| |
                      //    *----*----* | |
                      //    |    |    | | |
                      //    |    |    | | *
                      //    |    |    | |/
                      //    |    |    | *
                      //    |    |    |/
                      //    *----*----*
                      //

                      // first, create the new internal line
                      new_lines[0]->set (internal::Triangulation::
                                         TriaObject<1>(middle_vertex_index<dim,spacedim>(hex->face(4)),
                                                       middle_vertex_index<dim,spacedim>(hex->face(5))));

                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering:

                      // face 0: left plane
                      //       *
                      //      /|
                      //     * |
                      //    /| |
                      //   * | |
                      //   | 0 |
                      //   | | *
                      //   | |/
                      //   | *
                      //   |/
                      //   *
                      // face 1: right plane
                      //       *
                      //      /|
                      //     * |
                      //    /| |
                      //   * | |
                      //   | 1 |
                      //   | | *
                      //   | |/
                      //   | *
                      //   |/
                      //   *
                      // face 2: front plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   |   |
                      //   |   2   |
                      //   |   |   |
                      //   *-------*
                      // face 3: back plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   |   |
                      //   |   3   |
                      //   |   |   |
                      //   *---*---*
                      // face 4: bottom plane
                      //       *---*---*
                      //      /   5   /
                      //     *-6-*-7-*
                      //    /   4   /
                      //   *---*---*
                      // face 5: top plane
                      //       *---*---*
                      //      /   9   /
                      //     *10-*-11*
                      //    /   8   /
                      //   *---*---*
                      // middle planes
                      //     *-------*   *---*---*
                      //    /       /    |   |   |
                      //   /       /     |   12  |
                      //  /       /      |   |   |
                      // *-------*       *---*---*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_xy[13]
                      =
                      {
                        hex->face(0)->child(0)
                        ->line((hex->face(0)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //0
                        hex->face(1)->child(0)
                        ->line((hex->face(1)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //1
                        hex->face(2)->child(0)
                        ->line((hex->face(2)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //2
                        hex->face(3)->child(0)
                        ->line((hex->face(3)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //3

                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[4],f_fl[4],f_ro[4])),        //4
                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[4],f_fl[4],f_ro[4])),        //5
                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[4],f_fl[4],f_ro[4])),        //6
                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[4],f_fl[4],f_ro[4])),        //7

                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[5],f_fl[5],f_ro[5])),        //8
                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[5],f_fl[5],f_ro[5])),        //9
                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[5],f_fl[5],f_ro[5])),        //10
                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[5],f_fl[5],f_ro[5])),        //11

                        new_lines[0]                        //12
                      };

                      lines=&lines_xy[0];

                      unsigned int line_indices_xy[13];

                      for (unsigned int i=0; i<13; ++i)
                        line_indices_xy[i]=lines[i]->index();
                      line_indices=&line_indices_xy[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_xy[13];

                      // the middle vertices of the lines of our
                      // bottom face
                      const unsigned int middle_vertices[4]=
                      {
                        hex->line(0)->child(0)->vertex_index(1),
                        hex->line(1)->child(0)->vertex_index(1),
                        hex->line(2)->child(0)->vertex_index(1),
                        hex->line(3)->child(0)->vertex_index(1),
                      };

                      // note: for lines 0 to 3 the orientation of the
                      // line is 'true', if vertex 0 is on the bottom
                      // face
                      for (unsigned int i=0; i<4; ++i)
                        if (lines[i]->vertex_index(0)==middle_vertices[i])
                          line_orientation_xy[i]=true;
                        else
                          {
                            // it must be the other way round then
                            Assert(lines[i]->vertex_index(1)==middle_vertices[i],
                                   ExcInternalError());
                            line_orientation_xy[i]=false;
                          }

                      // note: for lines 4 to 11 (inner lines of the
                      // outer quads) the following holds: the second
                      // vertex of the even lines in standard
                      // orientation is the vertex in the middle of
                      // the quad, whereas for odd lines the first
                      // vertex is the same middle vertex.
                      for (unsigned int i=4; i<12; ++i)
                        if (lines[i]->vertex_index((i+1)%2) ==
                            middle_vertex_index<dim,spacedim>(hex->face(3+i/4)))
                          line_orientation_xy[i]=true;
                        else
                          {
                            // it must be the other way
                            // round then
                            Assert(lines[i]->vertex_index(i%2) ==
                                   (middle_vertex_index<dim,spacedim>(hex->face(3+i/4))),
                                   ExcInternalError());
                            line_orientation_xy[i]=false;
                          }
                      // for the last line the line orientation is
                      // always true, since it was just constructed
                      // that way

                      line_orientation_xy[12]=true;
                      line_orientation=&line_orientation_xy[0];

                      // set up the 4 quads, numbered as follows (left
                      // quad numbering, right line numbering
                      // extracted from above)
                      //
                      //      *          *
                      //     /|         9|
                      //    * |        * |
                      //  y/| |       8| 3
                      //  * |1|      * | |
                      //  | | |x     | 12|
                      //  |0| *      | | *
                      //  | |/       2 |5
                      //  | *        | *
                      //  |/         |4
                      //  *          *
                      //
                      //  x
                      //  *---*---*      *10-*-11*
                      //  |   |   |      |   |   |
                      //  | 2 | 3 |      0   12  1
                      //  |   |   |      |   |   |
                      //  *---*---*y     *-6-*-7-*

                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[2],
                                                         line_indices[12],
                                                         line_indices[4],
                                                         line_indices[8]));
                      new_quads[1]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[12],
                                                         line_indices[3],
                                                         line_indices[5],
                                                         line_indices[9]));
                      new_quads[2]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[6],
                                                         line_indices[10],
                                                         line_indices[0],
                                                         line_indices[12]));
                      new_quads[3]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[7],
                                                         line_indices[11],
                                                         line_indices[12],
                                                         line_indices[1]));

                      new_quads[0]->set_line_orientation(0,line_orientation[2]);
                      new_quads[0]->set_line_orientation(2,line_orientation[4]);
                      new_quads[0]->set_line_orientation(3,line_orientation[8]);

                      new_quads[1]->set_line_orientation(1,line_orientation[3]);
                      new_quads[1]->set_line_orientation(2,line_orientation[5]);
                      new_quads[1]->set_line_orientation(3,line_orientation[9]);

                      new_quads[2]->set_line_orientation(0,line_orientation[6]);
                      new_quads[2]->set_line_orientation(1,line_orientation[10]);
                      new_quads[2]->set_line_orientation(2,line_orientation[0]);

                      new_quads[3]->set_line_orientation(0,line_orientation[7]);
                      new_quads[3]->set_line_orientation(1,line_orientation[11]);
                      new_quads[3]->set_line_orientation(3,line_orientation[1]);

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    * | x
                      //   /| | *---*---*      *---------*
                      //  * |1| |   |   |     /         /
                      //  | | | | 2 | 3 |    /         /
                      //  |0| * |   |   |   /         /
                      //  | |/  *---*---*y *---------*x
                      //  | *
                      //  |/
                      //  *
                      //
                      // children of the faces of the old hex
                      //
                      //      *---*---*        *---*---*
                      //     /|   |   |       /18 / 19/|
                      //    * |10 | 11|      /---/---* |
                      //   /| |   |   |     /16 / 17/| |
                      //  * |5|   |   |    *---*---* |7|
                      //  | | *---*---*    |   |   | | *
                      //  |4|/14 / 15/     |   |   |6|/
                      //  | *---/---/      | 8 | 9 | *
                      //  |/12 / 13/       |   |   |/
                      //  *---*---*        *---*---*
                      //
                      // note that we have to take care of the
                      // orientation of faces.
                      const int quad_indices_xy[20]
                      =
                      {
                        new_quads[0]->index(),     //0
                        new_quads[1]->index(),
                        new_quads[2]->index(),
                        new_quads[3]->index(),

                        hex->face(0)->child_index(  child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),  //4
                        hex->face(0)->child_index(1-child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),

                        hex->face(1)->child_index(  child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),  //6
                        hex->face(1)->child_index(1-child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),

                        hex->face(2)->child_index(  child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),  //8
                        hex->face(2)->child_index(1-child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),

                        hex->face(3)->child_index(  child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),  //10
                        hex->face(3)->child_index(1-child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),

                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4])),  //12
                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[4],f_fl[4],f_ro[4])),
                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[4],f_fl[4],f_ro[4])),
                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4])),

                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5])),  //16
                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[5],f_fl[5],f_ro[5])),
                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[5],f_fl[5],f_ro[5])),
                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                      };
                      quad_indices=&quad_indices_xy[0];

                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[4],
                                                         quad_indices[0],
                                                         quad_indices[8],
                                                         quad_indices[2],
                                                         quad_indices[12],
                                                         quad_indices[16]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[0],
                                                         quad_indices[6],
                                                         quad_indices[9],
                                                         quad_indices[3],
                                                         quad_indices[13],
                                                         quad_indices[17]));
                      new_hexes[2]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[5],
                                                         quad_indices[1],
                                                         quad_indices[2],
                                                         quad_indices[10],
                                                         quad_indices[14],
                                                         quad_indices[18]));
                      new_hexes[3]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[1],
                                                         quad_indices[7],
                                                         quad_indices[3],
                                                         quad_indices[11],
                                                         quad_indices[15],
                                                         quad_indices[19]));
                      break;
                    }
                    case RefinementCase<dim>::cut_xz:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_xz
                      //
                      // the refined cube will look like this:
                      //
                      //        *----*----*
                      //       /    /    /|
                      //      /    /    / |
                      //     /    /    /  *
                      //    *----*----*  /|
                      //    |    |    | / |
                      //    |    |    |/  *
                      //    *----*----*  /
                      //    |    |    | /
                      //    |    |    |/
                      //    *----*----*
                      //

                      // first, create the new internal line
                      new_lines[0]->set (internal::Triangulation::
                                         TriaObject<1>(middle_vertex_index<dim,spacedim>(hex->face(2)),
                                                       middle_vertex_index<dim,spacedim>(hex->face(3))));

                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering:

                      // face 0: left plane
                      //       *
                      //      /|
                      //     / |
                      //    /  *
                      //   *  /|
                      //   | 0 |
                      //   |/  *
                      //   *  /
                      //   | /
                      //   |/
                      //   *
                      // face 1: right plane
                      //       *
                      //      /|
                      //     / |
                      //    /  *
                      //   *  /|
                      //   | 1 |
                      //   |/  *
                      //   *  /
                      //   | /
                      //   |/
                      //   *
                      // face 2: front plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   5   |
                      //   *-6-*-7-*
                      //   |   4   |
                      //   *---*---*
                      // face 3: back plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   9   |
                      //   *10-*-11*
                      //   |   8   |
                      //   *---*---*
                      // face 4: bottom plane
                      //       *---*---*
                      //      /   /   /
                      //     /   2   /
                      //    /   /   /
                      //   *---*---*
                      // face 5: top plane
                      //       *---*---*
                      //      /   /   /
                      //     /   3   /
                      //    /   /   /
                      //   *---*---*
                      // middle planes
                      //     *---*---*   *-------*
                      //    /   /   /    |       |
                      //   /   12  /     |       |
                      //  /   /   /      |       |
                      // *---*---*       *-------*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_xz[13]
                      =
                      {
                        hex->face(0)->child(0)
                        ->line((hex->face(0)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //0
                        hex->face(1)->child(0)
                        ->line((hex->face(1)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //1
                        hex->face(4)->child(0)
                        ->line((hex->face(4)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //2
                        hex->face(5)->child(0)
                        ->line((hex->face(5)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //3

                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[2],f_fl[2],f_ro[2])),        //4
                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[2],f_fl[2],f_ro[2])),        //5
                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[2],f_fl[2],f_ro[2])),        //6
                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[2],f_fl[2],f_ro[2])),        //7

                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[3],f_fl[3],f_ro[3])),        //8
                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[3],f_fl[3],f_ro[3])),        //9
                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[3],f_fl[3],f_ro[3])),        //10
                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[3],f_fl[3],f_ro[3])),        //11

                        new_lines[0]                        //12
                      };

                      lines=&lines_xz[0];

                      unsigned int line_indices_xz[13];

                      for (unsigned int i=0; i<13; ++i)
                        line_indices_xz[i]=lines[i]->index();
                      line_indices=&line_indices_xz[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_xz[13];

                      // the middle vertices of the
                      // lines of our front face
                      const unsigned int middle_vertices[4]=
                      {
                        hex->line(8)->child(0)->vertex_index(1),
                        hex->line(9)->child(0)->vertex_index(1),
                        hex->line(2)->child(0)->vertex_index(1),
                        hex->line(6)->child(0)->vertex_index(1),
                      };

                      // note: for lines 0 to 3 the orientation of the
                      // line is 'true', if vertex 0 is on the front
                      for (unsigned int i=0; i<4; ++i)
                        if (lines[i]->vertex_index(0)==middle_vertices[i])
                          line_orientation_xz[i]=true;
                        else
                          {
                            // it must be the other way round then
                            Assert(lines[i]->vertex_index(1)==middle_vertices[i],
                                   ExcInternalError());
                            line_orientation_xz[i]=false;
                          }

                      // note: for lines 4 to 11 (inner lines of the
                      // outer quads) the following holds: the second
                      // vertex of the even lines in standard
                      // orientation is the vertex in the middle of
                      // the quad, whereas for odd lines the first
                      // vertex is the same middle vertex.
                      for (unsigned int i=4; i<12; ++i)
                        if (lines[i]->vertex_index((i+1)%2) ==
                            middle_vertex_index<dim,spacedim>(hex->face(1+i/4)))
                          line_orientation_xz[i]=true;
                        else
                          {
                            // it must be the other way
                            // round then
                            Assert(lines[i]->vertex_index(i%2) ==
                                   (middle_vertex_index<dim,spacedim>(hex->face(1+i/4))),
                                   ExcInternalError());
                            line_orientation_xz[i]=false;
                          }
                      // for the last line the line orientation is
                      // always true, since it was just constructed
                      // that way

                      line_orientation_xz[12]=true;
                      line_orientation=&line_orientation_xz[0];

                      // set up the 4 quads, numbered as follows (left
                      // quad numbering, right line numbering
                      // extracted from above), the drawings denote
                      // middle planes
                      //
                      //      *          *
                      //     /|         /|
                      //    / |        3 9
                      //  y/  *       /  *
                      //  * 3/|      *  /|
                      //  | / |x     5 12|8
                      //  |/  *      |/  *
                      //  * 2/       *  /
                      //  | /        4 2
                      //  |/         |/
                      //  *          *
                      //
                      //       y
                      //      *----*----*      *-10-*-11-*
                      //     /    /    /      /    /    /
                      //    / 0  /  1 /      0    12   1
                      //   /    /    /      /    /    /
                      //  *----*----*x     *--6-*--7-*

                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[0],
                                                         line_indices[12],
                                                         line_indices[6],
                                                         line_indices[10]));
                      new_quads[1]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[12],
                                                         line_indices[1],
                                                         line_indices[7],
                                                         line_indices[11]));
                      new_quads[2]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[4],
                                                         line_indices[8],
                                                         line_indices[2],
                                                         line_indices[12]));
                      new_quads[3]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[5],
                                                         line_indices[9],
                                                         line_indices[12],
                                                         line_indices[3]));

                      new_quads[0]->set_line_orientation(0,line_orientation[0]);
                      new_quads[0]->set_line_orientation(2,line_orientation[6]);
                      new_quads[0]->set_line_orientation(3,line_orientation[10]);

                      new_quads[1]->set_line_orientation(1,line_orientation[1]);
                      new_quads[1]->set_line_orientation(2,line_orientation[7]);
                      new_quads[1]->set_line_orientation(3,line_orientation[11]);

                      new_quads[2]->set_line_orientation(0,line_orientation[4]);
                      new_quads[2]->set_line_orientation(1,line_orientation[8]);
                      new_quads[2]->set_line_orientation(2,line_orientation[2]);

                      new_quads[3]->set_line_orientation(0,line_orientation[5]);
                      new_quads[3]->set_line_orientation(1,line_orientation[9]);
                      new_quads[3]->set_line_orientation(3,line_orientation[3]);

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    / | x
                      //   /3 * *-------*      *----*----*
                      //  *  /| |       |     /    /    /
                      //  | / | |       |    /  0 /  1 /
                      //  |/  * |       |   /    /    /
                      //  * 2/  *-------*y *----*----*x
                      //  | /
                      //  |/
                      //  *
                      //
                      // children of the faces
                      // of the old hex
                      //      *---*---*        *---*---*
                      //     /|13 | 15|       /   /   /|
                      //    / |   |   |      /18 / 19/ |
                      //   /  *---*---*     /   /   /  *
                      //  * 5/|   |   |    *---*---* 7/|
                      //  | / |12 | 14|    | 9 | 11| / |
                      //  |/4 *---*---*    |   |   |/6 *
                      //  *  /   /   /     *---*---*  /
                      //  | /16 / 17/      |   |   | /
                      //  |/   /   /       | 8 | 10|/
                      //  *---*---*        *---*---*
                      //
                      // note that we have to take care of the
                      // orientation of faces.
                      const int quad_indices_xz[20]
                      =
                      {
                        new_quads[0]->index(),     //0
                        new_quads[1]->index(),
                        new_quads[2]->index(),
                        new_quads[3]->index(),

                        hex->face(0)->child_index(  child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),  //4
                        hex->face(0)->child_index(1-child_at_origin[hex->face(0)->refinement_case()-1][f_fl[0]][f_ro[0]]),

                        hex->face(1)->child_index(  child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),  //6
                        hex->face(1)->child_index(1-child_at_origin[hex->face(1)->refinement_case()-1][f_fl[1]][f_ro[1]]),

                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2])),  //8
                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[2],f_fl[2],f_ro[2])),
                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[2],f_fl[2],f_ro[2])),
                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2])),

                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3])),  //12
                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[3],f_fl[3],f_ro[3])),
                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[3],f_fl[3],f_ro[3])),
                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3])),

                        hex->face(4)->child_index(  child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),  //16
                        hex->face(4)->child_index(1-child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),

                        hex->face(5)->child_index(  child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]]),  //18
                        hex->face(5)->child_index(1-child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]])
                      };
                      quad_indices=&quad_indices_xz[0];

                      // due to the exchange of x and y for the front
                      // and back face, we order the children
                      // according to
                      //
                      // *---*---*
                      // | 1 | 3 |
                      // *---*---*
                      // | 0 | 2 |
                      // *---*---*
                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[4],
                                                         quad_indices[2],
                                                         quad_indices[8],
                                                         quad_indices[12],
                                                         quad_indices[16],
                                                         quad_indices[0]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[5],
                                                         quad_indices[3],
                                                         quad_indices[9],
                                                         quad_indices[13],
                                                         quad_indices[0],
                                                         quad_indices[18]));
                      new_hexes[2]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[2],
                                                         quad_indices[6],
                                                         quad_indices[10],
                                                         quad_indices[14],
                                                         quad_indices[17],
                                                         quad_indices[1]));
                      new_hexes[3]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[3],
                                                         quad_indices[7],
                                                         quad_indices[11],
                                                         quad_indices[15],
                                                         quad_indices[1],
                                                         quad_indices[19]));
                      break;
                    }
                    case RefinementCase<dim>::cut_yz:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_yz
                      //
                      // the refined cube will look like this:
                      //
                      //        *---------*
                      //       /         /|
                      //      *---------* |
                      //     /         /| |
                      //    *---------* |/|
                      //    |         | * |
                      //    |         |/| *
                      //    *---------* |/
                      //    |         | *
                      //    |         |/
                      //    *---------*
                      //

                      // first, create the new
                      // internal line
                      new_lines[0]->set (internal::Triangulation::
                                         TriaObject<1>(middle_vertex_index<dim,spacedim>(hex->face(0)),
                                                       middle_vertex_index<dim,spacedim>(hex->face(1))));

                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering: (note that face 0 and 1 each are
                      // shown twice for better readability)

                      // face 0: left plane
                      //       *            *
                      //      /|           /|
                      //     * |          * |
                      //    /| *         /| *
                      //   * 5/|        * |7|
                      //   | * |        | * |
                      //   |/| *        |6| *
                      //   * 4/         * |/
                      //   | *          | *
                      //   |/           |/
                      //   *            *
                      // face 1: right plane
                      //       *            *
                      //      /|           /|
                      //     * |          * |
                      //    /| *         /| *
                      //   * 9/|        * |11
                      //   | * |        | * |
                      //   |/| *        |10 *
                      //   * 8/         * |/
                      //   | *          | *
                      //   |/           |/
                      //   *            *
                      // face 2: front plane
                      //   (note: x,y exchanged)
                      //   *-------*
                      //   |       |
                      //   *---0---*
                      //   |       |
                      //   *-------*
                      // face 3: back plane
                      //   (note: x,y exchanged)
                      //   *-------*
                      //   |       |
                      //   *---1---*
                      //   |       |
                      //   *-------*
                      // face 4: bottom plane
                      //       *-------*
                      //      /       /
                      //     *---2---*
                      //    /       /
                      //   *-------*
                      // face 5: top plane
                      //       *-------*
                      //      /       /
                      //     *---3---*
                      //    /       /
                      //   *-------*
                      // middle planes
                      //     *-------*   *-------*
                      //    /       /    |       |
                      //   *---12--*     |       |
                      //  /       /      |       |
                      // *-------*       *-------*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_yz[13]
                      =
                      {
                        hex->face(2)->child(0)
                        ->line((hex->face(2)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //0
                        hex->face(3)->child(0)
                        ->line((hex->face(3)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //1
                        hex->face(4)->child(0)
                        ->line((hex->face(4)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //2
                        hex->face(5)->child(0)
                        ->line((hex->face(5)->refinement_case() == RefinementCase<2>::cut_x) ? 1 : 3),        //3

                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[0],f_fl[0],f_ro[0])),        //4
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[0],f_fl[0],f_ro[0])),        //5
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[0],f_fl[0],f_ro[0])),        //6
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[0],f_fl[0],f_ro[0])),        //7

                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[1],f_fl[1],f_ro[1])),        //8
                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[1],f_fl[1],f_ro[1])),        //9
                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[1],f_fl[1],f_ro[1])),        //10
                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[1],f_fl[1],f_ro[1])),        //11

                        new_lines[0]                        //12
                      };

                      lines=&lines_yz[0];

                      unsigned int line_indices_yz[13];

                      for (unsigned int i=0; i<13; ++i)
                        line_indices_yz[i]=lines[i]->index();
                      line_indices=&line_indices_yz[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_yz[13];

                      // the middle vertices of the lines of our front
                      // face
                      const unsigned int middle_vertices[4]=
                      {
                        hex->line(8)->child(0)->vertex_index(1),
                        hex->line(10)->child(0)->vertex_index(1),
                        hex->line(0)->child(0)->vertex_index(1),
                        hex->line(4)->child(0)->vertex_index(1),
                      };

                      // note: for lines 0 to 3 the orientation of the
                      // line is 'true', if vertex 0 is on the front
                      for (unsigned int i=0; i<4; ++i)
                        if (lines[i]->vertex_index(0)==middle_vertices[i])
                          line_orientation_yz[i]=true;
                        else
                          {
                            // it must be the other way round then
                            Assert(lines[i]->vertex_index(1)==middle_vertices[i],
                                   ExcInternalError());
                            line_orientation_yz[i]=false;
                          }

                      // note: for lines 4 to 11 (inner lines of the
                      // outer quads) the following holds: the second
                      // vertex of the even lines in standard
                      // orientation is the vertex in the middle of
                      // the quad, whereas for odd lines the first
                      // vertex is the same middle vertex.
                      for (unsigned int i=4; i<12; ++i)
                        if (lines[i]->vertex_index((i+1)%2) ==
                            middle_vertex_index<dim,spacedim>(hex->face(i/4-1)))
                          line_orientation_yz[i]=true;
                        else
                          {
                            // it must be the other way
                            // round then
                            Assert(lines[i]->vertex_index(i%2) ==
                                   (middle_vertex_index<dim,spacedim>(hex->face(i/4-1))),
                                   ExcInternalError());
                            line_orientation_yz[i]=false;
                          }
                      // for the last line the line orientation is
                      // always true, since it was just constructed
                      // that way

                      line_orientation_yz[12]=true;
                      line_orientation=&line_orientation_yz[0];

                      // set up the 4 quads, numbered as follows (left
                      // quad numbering, right line numbering
                      // extracted from above)
                      //
                      //  x
                      //  *-------*      *---3---*
                      //  |   3   |      5       9
                      //  *-------*      *---12--*
                      //  |   2   |      4       8
                      //  *-------*y     *---2---*
                      //
                      //       y
                      //      *---------*      *----1----*
                      //     /    1    /      7         11
                      //    *---------*      *----12---*
                      //   /    0    /      6         10
                      //  *---------*x     *----0----*

                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[6],
                                                         line_indices[10],
                                                         line_indices[0],
                                                         line_indices[12]));
                      new_quads[1]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[7],
                                                         line_indices[11],
                                                         line_indices[12],
                                                         line_indices[1]));
                      new_quads[2]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[2],
                                                         line_indices[12],
                                                         line_indices[4],
                                                         line_indices[8]));
                      new_quads[3]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[12],
                                                         line_indices[3],
                                                         line_indices[5],
                                                         line_indices[9]));

                      new_quads[0]->set_line_orientation(0,line_orientation[6]);
                      new_quads[0]->set_line_orientation(1,line_orientation[10]);
                      new_quads[0]->set_line_orientation(2,line_orientation[0]);

                      new_quads[1]->set_line_orientation(0,line_orientation[7]);
                      new_quads[1]->set_line_orientation(1,line_orientation[11]);
                      new_quads[1]->set_line_orientation(3,line_orientation[1]);

                      new_quads[2]->set_line_orientation(0,line_orientation[2]);
                      new_quads[2]->set_line_orientation(2,line_orientation[4]);
                      new_quads[2]->set_line_orientation(3,line_orientation[8]);

                      new_quads[3]->set_line_orientation(1,line_orientation[3]);
                      new_quads[3]->set_line_orientation(2,line_orientation[5]);
                      new_quads[3]->set_line_orientation(3,line_orientation[9]);

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    / | x
                      //   /  | *-------*      *---------*
                      //  *   | |   3   |     /    1    /
                      //  |   | *-------*    *---------*
                      //  |   * |   2   |   /    0    /
                      //  |  /  *-------*y *---------*x
                      //  | /
                      //  |/
                      //  *
                      //
                      // children of the faces
                      // of the old hex
                      //      *-------*        *-------*
                      //     /|       |       /  19   /|
                      //    * |  15   |      *-------* |
                      //   /|7*-------*     /  18   /|11
                      //  * |/|       |    *-------* |/|
                      //  |6* |  14   |    |       10* |
                      //  |/|5*-------*    |  13   |/|9*
                      //  * |/  17   /     *-------* |/
                      //  |4*-------*      |       |8*
                      //  |/  16   /       |  12   |/
                      //  *-------*        *-------*
                      //
                      // note that we have to take care of the
                      // orientation of faces.
                      const int quad_indices_yz[20]
                      =
                      {
                        new_quads[0]->index(),     //0
                        new_quads[1]->index(),
                        new_quads[2]->index(),
                        new_quads[3]->index(),

                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0])),  //4
                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[0],f_fl[0],f_ro[0])),
                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[0],f_fl[0],f_ro[0])),
                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0])),

                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1])),  //8
                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[1],f_fl[1],f_ro[1])),
                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[1],f_fl[1],f_ro[1])),
                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1])),

                        hex->face(2)->child_index(  child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),  //12
                        hex->face(2)->child_index(1-child_at_origin[hex->face(2)->refinement_case()-1][f_fl[2]][f_ro[2]]),

                        hex->face(3)->child_index(  child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),  //14
                        hex->face(3)->child_index(1-child_at_origin[hex->face(3)->refinement_case()-1][f_fl[3]][f_ro[3]]),

                        hex->face(4)->child_index(  child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),  //16
                        hex->face(4)->child_index(1-child_at_origin[hex->face(4)->refinement_case()-1][f_fl[4]][f_ro[4]]),

                        hex->face(5)->child_index(  child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]]),  //18
                        hex->face(5)->child_index(1-child_at_origin[hex->face(5)->refinement_case()-1][f_fl[5]][f_ro[5]])
                      };
                      quad_indices=&quad_indices_yz[0];

                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[4],
                                                         quad_indices[8],
                                                         quad_indices[12],
                                                         quad_indices[2],
                                                         quad_indices[16],
                                                         quad_indices[0]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[5],
                                                         quad_indices[9],
                                                         quad_indices[2],
                                                         quad_indices[14],
                                                         quad_indices[17],
                                                         quad_indices[1]));
                      new_hexes[2]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[6],
                                                         quad_indices[10],
                                                         quad_indices[13],
                                                         quad_indices[3],
                                                         quad_indices[0],
                                                         quad_indices[18]));
                      new_hexes[3]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[7],
                                                         quad_indices[11],
                                                         quad_indices[3],
                                                         quad_indices[15],
                                                         quad_indices[1],
                                                         quad_indices[19]));
                      break;
                    }
                    case RefinementCase<dim>::cut_xyz:
                    {
                      //////////////////////////////
                      //
                      //     RefinementCase<dim>::cut_xyz
                      //     isotropic refinement
                      //
                      // the refined cube will look
                      // like this:
                      //
                      //        *----*----*
                      //       /    /    /|
                      //      *----*----* |
                      //     /    /    /| *
                      //    *----*----* |/|
                      //    |    |    | * |
                      //    |    |    |/| *
                      //    *----*----* |/
                      //    |    |    | *
                      //    |    |    |/
                      //    *----*----*
                      //

                      // find the next unused vertex and set it
                      // appropriately
                      while (triangulation.vertices_used[next_unused_vertex] == true)
                        ++next_unused_vertex;
                      Assert (next_unused_vertex < triangulation.vertices.size(),
                              ExcTooFewVerticesAllocated());
                      triangulation.vertices_used[next_unused_vertex] = true;

                      // the new vertex is definitely in the interior,
                      // so we need not worry about the
                      // boundary. However we need to worry about
                      // Manifolds. Let the cell compute its own
                      // center, by querying the underlying manifold
                      // object.
                      triangulation.vertices[next_unused_vertex] =
                        hex->center(true, true);

                      // set the data of the six lines.  first collect
                      // the indices of the seven vertices (consider
                      // the two planes to be crossed to form the
                      // planes cutting the hex in two vertically and
                      // horizontally)
                      //
                      //     *--3--*   *--5--*
                      //    /  /  /    |  |  |
                      //   0--6--1     0--6--1
                      //  /  /  /      |  |  |
                      // *--2--*       *--4--*
                      // the lines are numbered
                      // as follows:
                      //     *--*--*   *--*--*
                      //    /  1  /    |  5  |
                      //   *2-*-3*     *2-*-3*
                      //  /  0  /      |  4  |
                      // *--*--*       *--*--*
                      //
                      const unsigned int vertex_indices_xyz[7]
                        = { middle_vertex_index<dim,spacedim>(hex->face(0)),
                            middle_vertex_index<dim,spacedim>(hex->face(1)),
                            middle_vertex_index<dim,spacedim>(hex->face(2)),
                            middle_vertex_index<dim,spacedim>(hex->face(3)),
                            middle_vertex_index<dim,spacedim>(hex->face(4)),
                            middle_vertex_index<dim,spacedim>(hex->face(5)),
                            next_unused_vertex
                          };
                      vertex_indices=&vertex_indices_xyz[0];

                      new_lines[0]->set (internal::Triangulation::
                                         TriaObject<1>(vertex_indices[2], vertex_indices[6]));
                      new_lines[1]->set (internal::Triangulation::
                                         TriaObject<1>(vertex_indices[6], vertex_indices[3]));
                      new_lines[2]->set (internal::Triangulation::
                                         TriaObject<1>(vertex_indices[0], vertex_indices[6]));
                      new_lines[3]->set (internal::Triangulation::
                                         TriaObject<1>(vertex_indices[6], vertex_indices[1]));
                      new_lines[4]->set (internal::Triangulation::
                                         TriaObject<1>(vertex_indices[4], vertex_indices[6]));
                      new_lines[5]->set (internal::Triangulation::
                                         TriaObject<1>(vertex_indices[6], vertex_indices[5]));

                      // again, first collect some data about the
                      // indices of the lines, with the following
                      // numbering: (note that face 0 and 1 each are
                      // shown twice for better readability)

                      // face 0: left plane
                      //       *            *
                      //      /|           /|
                      //     * |          * |
                      //    /| *         /| *
                      //   * 1/|        * |3|
                      //   | * |        | * |
                      //   |/| *        |2| *
                      //   * 0/         * |/
                      //   | *          | *
                      //   |/           |/
                      //   *            *
                      // face 1: right plane
                      //       *            *
                      //      /|           /|
                      //     * |          * |
                      //    /| *         /| *
                      //   * 5/|        * |7|
                      //   | * |        | * |
                      //   |/| *        |6| *
                      //   * 4/         * |/
                      //   | *          | *
                      //   |/           |/
                      //   *            *
                      // face 2: front plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   11  |
                      //   *-8-*-9-*
                      //   |   10  |
                      //   *---*---*
                      // face 3: back plane
                      //   (note: x,y exchanged)
                      //   *---*---*
                      //   |   15  |
                      //   *12-*-13*
                      //   |   14  |
                      //   *---*---*
                      // face 4: bottom plane
                      //       *---*---*
                      //      /  17   /
                      //     *18-*-19*
                      //    /   16  /
                      //   *---*---*
                      // face 5: top plane
                      //       *---*---*
                      //      /  21   /
                      //     *22-*-23*
                      //    /   20  /
                      //   *---*---*
                      // middle planes
                      //     *---*---*   *---*---*
                      //    /  25   /    |   29  |
                      //   *26-*-27*     *26-*-27*
                      //  /   24  /      |   28  |
                      // *---*---*       *---*---*

                      // set up a list of line iterators first. from
                      // this, construct lists of line_indices and
                      // line orientations later on
                      const typename Triangulation<dim,spacedim>::raw_line_iterator
                      lines_xyz[30]
                      =
                      {
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[0],f_fl[0],f_ro[0])),        //0
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[0],f_fl[0],f_ro[0])),        //1
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[0],f_fl[0],f_ro[0])),        //2
                        hex->face(0)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[0],f_fl[0],f_ro[0])),        //3

                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[1],f_fl[1],f_ro[1])),        //4
                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[1],f_fl[1],f_ro[1])),        //5
                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[1],f_fl[1],f_ro[1])),        //6
                        hex->face(1)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[1],f_fl[1],f_ro[1])),        //7

                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[2],f_fl[2],f_ro[2])),        //8
                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[2],f_fl[2],f_ro[2])),        //9
                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[2],f_fl[2],f_ro[2])),        //10
                        hex->face(2)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[2],f_fl[2],f_ro[2])),        //11

                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[3],f_fl[3],f_ro[3])),        //12
                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[3],f_fl[3],f_ro[3])),        //13
                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[3],f_fl[3],f_ro[3])),        //14
                        hex->face(3)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[3],f_fl[3],f_ro[3])),        //15

                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[4],f_fl[4],f_ro[4])),        //16
                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[4],f_fl[4],f_ro[4])),        //17
                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[4],f_fl[4],f_ro[4])),        //18
                        hex->face(4)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[4],f_fl[4],f_ro[4])),        //19

                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[5],f_fl[5],f_ro[5])),        //20
                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[5],f_fl[5],f_ro[5])),        //21
                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[5],f_fl[5],f_ro[5])),        //22
                        hex->face(5)->isotropic_child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                        ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[5],f_fl[5],f_ro[5])),        //23

                        new_lines[0],                       //24
                        new_lines[1],                       //25
                        new_lines[2],                       //26
                        new_lines[3],                       //27
                        new_lines[4],                       //28
                        new_lines[5]                        //29
                      };

                      lines=&lines_xyz[0];

                      unsigned int line_indices_xyz[30];
                      for (unsigned int i=0; i<30; ++i)
                        line_indices_xyz[i]=lines[i]->index();
                      line_indices=&line_indices_xyz[0];

                      // the orientation of lines for the inner quads
                      // is quite tricky. as these lines are newly
                      // created ones and thus have no parents, they
                      // cannot inherit this property. set up an array
                      // and fill it with the respective values
                      bool line_orientation_xyz[30];

                      // note: for the first 24 lines (inner lines of
                      // the outer quads) the following holds: the
                      // second vertex of the even lines in standard
                      // orientation is the vertex in the middle of
                      // the quad, whereas for odd lines the first
                      // vertex is the same middle vertex.
                      for (unsigned int i=0; i<24; ++i)
                        if (lines[i]->vertex_index((i+1)%2)==vertex_indices[i/4])
                          line_orientation_xyz[i]=true;
                        else
                          {
                            // it must be the other way
                            // round then
                            Assert(lines[i]->vertex_index(i%2)==vertex_indices[i/4],
                                   ExcInternalError());
                            line_orientation_xyz[i]=false;
                          }
                      // for the last 6 lines the line orientation is
                      // always true, since they were just constructed
                      // that way
                      for (unsigned int i=24; i<30; ++i)
                        line_orientation_xyz[i]=true;
                      line_orientation=&line_orientation_xyz[0];

                      // set up the 12 quads, numbered as follows
                      // (left quad numbering, right line numbering
                      // extracted from above)
                      //
                      //      *          *
                      //     /|        21|
                      //    * |        * 15
                      //  y/|3*      20| *
                      //  * |/|      * |/|
                      //  |2* |x    11 * 14
                      //  |/|1*      |/| *
                      //  * |/       * |17
                      //  |0*       10 *
                      //  |/         |16
                      //  *          *
                      //
                      //  x
                      //  *---*---*      *22-*-23*
                      //  | 5 | 7 |      1  29   5
                      //  *---*---*      *26-*-27*
                      //  | 4 | 6 |      0  28   4
                      //  *---*---*y     *18-*-19*
                      //
                      //       y
                      //      *----*----*      *-12-*-13-*
                      //     / 10 / 11 /      3   25    7
                      //    *----*----*      *-26-*-27-*
                      //   / 8  / 9  /      2   24    6
                      //  *----*----*x     *--8-*--9-*

                      new_quads[0]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[10],
                                                         line_indices[28],
                                                         line_indices[16],
                                                         line_indices[24]));
                      new_quads[1]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[28],
                                                         line_indices[14],
                                                         line_indices[17],
                                                         line_indices[25]));
                      new_quads[2]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[11],
                                                         line_indices[29],
                                                         line_indices[24],
                                                         line_indices[20]));
                      new_quads[3]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[29],
                                                         line_indices[15],
                                                         line_indices[25],
                                                         line_indices[21]));
                      new_quads[4]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[18],
                                                         line_indices[26],
                                                         line_indices[0],
                                                         line_indices[28]));
                      new_quads[5]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[26],
                                                         line_indices[22],
                                                         line_indices[1],
                                                         line_indices[29]));
                      new_quads[6]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[19],
                                                         line_indices[27],
                                                         line_indices[28],
                                                         line_indices[4]));
                      new_quads[7]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[27],
                                                         line_indices[23],
                                                         line_indices[29],
                                                         line_indices[5]));
                      new_quads[8]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[2],
                                                         line_indices[24],
                                                         line_indices[8],
                                                         line_indices[26]));
                      new_quads[9]->set (internal::Triangulation
                                         ::TriaObject<2>(line_indices[24],
                                                         line_indices[6],
                                                         line_indices[9],
                                                         line_indices[27]));
                      new_quads[10]->set (internal::Triangulation
                                          ::TriaObject<2>(line_indices[3],
                                                          line_indices[25],
                                                          line_indices[26],
                                                          line_indices[12]));
                      new_quads[11]->set (internal::Triangulation
                                          ::TriaObject<2>(line_indices[25],
                                                          line_indices[7],
                                                          line_indices[27],
                                                          line_indices[13]));

                      // now reset the line_orientation flags of outer
                      // lines as they cannot be set in a loop (at
                      // least not easily)
                      new_quads[0]->set_line_orientation(0,line_orientation[10]);
                      new_quads[0]->set_line_orientation(2,line_orientation[16]);

                      new_quads[1]->set_line_orientation(1,line_orientation[14]);
                      new_quads[1]->set_line_orientation(2,line_orientation[17]);

                      new_quads[2]->set_line_orientation(0,line_orientation[11]);
                      new_quads[2]->set_line_orientation(3,line_orientation[20]);

                      new_quads[3]->set_line_orientation(1,line_orientation[15]);
                      new_quads[3]->set_line_orientation(3,line_orientation[21]);

                      new_quads[4]->set_line_orientation(0,line_orientation[18]);
                      new_quads[4]->set_line_orientation(2,line_orientation[0]);

                      new_quads[5]->set_line_orientation(1,line_orientation[22]);
                      new_quads[5]->set_line_orientation(2,line_orientation[1]);

                      new_quads[6]->set_line_orientation(0,line_orientation[19]);
                      new_quads[6]->set_line_orientation(3,line_orientation[4]);

                      new_quads[7]->set_line_orientation(1,line_orientation[23]);
                      new_quads[7]->set_line_orientation(3,line_orientation[5]);

                      new_quads[8]->set_line_orientation(0,line_orientation[2]);
                      new_quads[8]->set_line_orientation(2,line_orientation[8]);

                      new_quads[9]->set_line_orientation(1,line_orientation[6]);
                      new_quads[9]->set_line_orientation(2,line_orientation[9]);

                      new_quads[10]->set_line_orientation(0,line_orientation[3]);
                      new_quads[10]->set_line_orientation(3,line_orientation[12]);

                      new_quads[11]->set_line_orientation(1,line_orientation[7]);
                      new_quads[11]->set_line_orientation(3,line_orientation[13]);

                      /////////////////////////////////
                      // create the eight new hexes
                      //
                      // again first collect some data.  here, we need
                      // the indices of a whole lotta quads.

                      // the quads are numbered as follows:
                      //
                      // planes in the interior of the old hex:
                      //
                      //      *
                      //     /|
                      //    * |
                      //   /|3*  *---*---*      *----*----*
                      //  * |/|  | 5 | 7 |     / 10 / 11 /
                      //  |2* |  *---*---*    *----*----*
                      //  |/|1*  | 4 | 6 |   / 8  / 9  /
                      //  * |/   *---*---*y *----*----*x
                      //  |0*
                      //  |/
                      //  *
                      //
                      // children of the faces
                      // of the old hex
                      //      *-------*        *-------*
                      //     /|25   27|       /34   35/|
                      //    15|       |      /       /19
                      //   /  |       |     /32   33/  |
                      //  *   |24   26|    *-------*18 |
                      //  1413*-------*    |21   23| 17*
                      //  |  /30   31/     |       |  /
                      //  12/       /      |       |16
                      //  |/28   29/       |20   22|/
                      //  *-------*        *-------*
                      //
                      // note that we have to
                      // take care of the
                      // orientation of
                      // faces.
                      const int quad_indices_xyz[36]
                      =
                      {
                        new_quads[0]->index(),     //0
                        new_quads[1]->index(),
                        new_quads[2]->index(),
                        new_quads[3]->index(),
                        new_quads[4]->index(),
                        new_quads[5]->index(),
                        new_quads[6]->index(),
                        new_quads[7]->index(),
                        new_quads[8]->index(),
                        new_quads[9]->index(),
                        new_quads[10]->index(),
                        new_quads[11]->index(),    //11

                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0])),  //12
                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[0],f_fl[0],f_ro[0])),
                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[0],f_fl[0],f_ro[0])),
                        hex->face(0)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0])),

                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1])),  //16
                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[1],f_fl[1],f_ro[1])),
                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[1],f_fl[1],f_ro[1])),
                        hex->face(1)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1])),

                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2])),  //20
                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[2],f_fl[2],f_ro[2])),
                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[2],f_fl[2],f_ro[2])),
                        hex->face(2)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2])),

                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3])),  //24
                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[3],f_fl[3],f_ro[3])),
                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[3],f_fl[3],f_ro[3])),
                        hex->face(3)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3])),

                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4])),  //28
                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[4],f_fl[4],f_ro[4])),
                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[4],f_fl[4],f_ro[4])),
                        hex->face(4)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4])),

                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5])),  //32
                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[5],f_fl[5],f_ro[5])),
                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[5],f_fl[5],f_ro[5])),
                        hex->face(5)->isotropic_child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                      };
                      quad_indices=&quad_indices_xyz[0];

                      // bottom children
                      new_hexes[0]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[12],
                                                         quad_indices[0],
                                                         quad_indices[20],
                                                         quad_indices[4],
                                                         quad_indices[28],
                                                         quad_indices[8]));
                      new_hexes[1]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[0],
                                                         quad_indices[16],
                                                         quad_indices[22],
                                                         quad_indices[6],
                                                         quad_indices[29],
                                                         quad_indices[9]));
                      new_hexes[2]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[13],
                                                         quad_indices[1],
                                                         quad_indices[4],
                                                         quad_indices[24],
                                                         quad_indices[30],
                                                         quad_indices[10]));
                      new_hexes[3]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[1],
                                                         quad_indices[17],
                                                         quad_indices[6],
                                                         quad_indices[26],
                                                         quad_indices[31],
                                                         quad_indices[11]));

                      // top children
                      new_hexes[4]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[14],
                                                         quad_indices[2],
                                                         quad_indices[21],
                                                         quad_indices[5],
                                                         quad_indices[8],
                                                         quad_indices[32]));
                      new_hexes[5]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[2],
                                                         quad_indices[18],
                                                         quad_indices[23],
                                                         quad_indices[7],
                                                         quad_indices[9],
                                                         quad_indices[33]));
                      new_hexes[6]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[15],
                                                         quad_indices[3],
                                                         quad_indices[5],
                                                         quad_indices[25],
                                                         quad_indices[10],
                                                         quad_indices[34]));
                      new_hexes[7]->set (internal::Triangulation
                                         ::TriaObject<3>(quad_indices[3],
                                                         quad_indices[19],
                                                         quad_indices[7],
                                                         quad_indices[27],
                                                         quad_indices[11],
                                                         quad_indices[35]));
                      break;
                    }
                    default:
                      // all refinement cases have been treated, there
                      // only remains
                      // RefinementCase<dim>::no_refinement as
                      // untreated enumeration value. However, in that
                      // case we should have aborted much
                      // earlier. thus we should never get here
                      Assert(false, ExcInternalError());
                      break;
                    }//switch (ref_case)

                  // and set face orientation flags. note that new
                  // faces in the interior of the mother cell always
                  // have a correctly oriented face, but the ones on
                  // the outer faces will inherit this flag
                  //
                  // the flag have been set to true for all faces
                  // initially, now go the other way round and reset
                  // faces that are at the boundary of the mother cube
                  //
                  // the same is true for the face_flip and
                  // face_rotation flags. however, the latter two are
                  // set to false by default as this is the standard
                  // value

                  // loop over all faces and all (relevant) subfaces
                  // of that in order to set the correct values for
                  // face_orientation, face_flip and face_rotation,
                  // which are inherited from the corresponding face
                  // of the mother cube
                  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                    for (unsigned int s=0;
                         s<std::max(GeometryInfo<dim-1>::n_children(GeometryInfo<dim>::face_refinement_case(ref_case,f)),
                                    1U);
                         ++s)
                      {
                        const unsigned int current_child
                          =GeometryInfo<dim>::child_cell_on_face(ref_case,
                                                                 f,
                                                                 s,
                                                                 f_or[f],
                                                                 f_fl[f],
                                                                 f_ro[f],
                                                                 GeometryInfo<dim>::face_refinement_case(ref_case,
                                                                     f,
                                                                     f_or[f],
                                                                     f_fl[f],
                                                                     f_ro[f]));
                        new_hexes[current_child]->set_face_orientation (f, f_or[f]);
                        new_hexes[current_child]->set_face_flip        (f, f_fl[f]);
                        new_hexes[current_child]->set_face_rotation    (f, f_ro[f]);
                      }

                  // now see if we have created cells that are
                  // distorted and if so add them to our list
                  if ((check_for_distorted_cells == true)
                      &&
                      has_distorted_children (hex,
                                              internal::int2type<dim>(),
                                              internal::int2type<spacedim>()))
                    cells_with_distorted_children.distorted_cells.push_back (hex);

                  // note that the refinement flag was already cleared
                  // at the beginning of this loop
                }
          }

        // clear user data on quads. we used some of this data to
        // indicate anisotropic refinemnt cases on faces. all data
        // should be cleared by now, but the information whether we
        // used indices or pointers is still present. reset it now to
        // enable the user to use whichever he likes later on.
        triangulation.faces->quads.clear_user_data();

        // return the list with distorted children
        return cells_with_distorted_children;
      }


      /**
       * At the boundary of the domain, the new point on the face may
       * be far inside the current cell, if the boundary has a strong
       * curvature. If we allow anisotropic refinement here, the
       * resulting cell may be strongly distorted. To prevent this,
       * this function flags such cells for isotropic refinement. It
       * is called automatically from
       * prepare_coarsening_and_refinement().
       *
       * This function does nothing in 1d (therefore the
       * specialization).
       */
      template <int spacedim>
      static
      void
      prevent_distorted_boundary_cells (const Triangulation<1,spacedim> &);


      template <int dim, int spacedim>
      static
      void
      prevent_distorted_boundary_cells (Triangulation<dim,spacedim> &triangulation)
      {
        // If the codimension is one, we cannot perform this check
        // yet.
        if (spacedim>dim) return;

        for (typename Triangulation<dim,spacedim>::cell_iterator
             cell=triangulation.begin(); cell!=triangulation.end(); ++cell)
          if (cell->at_boundary() &&
              cell->refine_flag_set() &&
              cell->refine_flag_set()!=RefinementCase<dim>::isotropic_refinement)
            {
              // The cell is at the boundary and it is flagged for
              // anisotropic refinement. Therefore, we have a closer
              // look
              const RefinementCase<dim> ref_case=cell->refine_flag_set();
              for (unsigned int face_no=0;
                   face_no<GeometryInfo<dim>::faces_per_cell;
                   ++face_no)
                if (cell->face(face_no)->at_boundary())
                  {
                    // this is the critical face at the boundary.
                    if (GeometryInfo<dim>::face_refinement_case(ref_case,face_no)
                        !=RefinementCase<dim-1>::isotropic_refinement)
                      {
                        // up to now, we do not want to refine this
                        // cell along the face under consideration
                        // here.
                        const typename Triangulation<dim,spacedim>::face_iterator
                        face = cell->face(face_no);
                        // the new point on the boundary would be this
                        // one.
                        const Point<spacedim> new_bound
                          = face->center(true);
                        // to check it, transform to the unit cell
                        // with Q1Mapping
                        const Point<dim> new_unit
                          = StaticMappingQ1<dim,spacedim>::mapping.
                            transform_real_to_unit_cell(cell,
                                                        new_bound);

                        // Now, we have to calculate the distance from
                        // the face in the unit cell.

                        // take the correct coordinate direction (0
                        // for faces 0 and 1, 1 for faces 2 and 3, 2
                        // for faces 4 and 5) and subtract the correct
                        // boundary value of the face (0 for faces 0,
                        // 2, and 4; 1 for faces 1, 3 and 5)
                        const double dist = std::fabs(new_unit[face_no/2] - face_no%2);

                        // compare this with the empirical value
                        // allowed. if it is too big, flag the face
                        // for isotropic refinement
                        const double allowed=0.25;

                        if (dist>allowed)
                          cell->flag_for_face_refinement(face_no);
                      }//if flagged for anistropic refinement
                  }//if (cell->face(face)->at_boundary())
            }//for all cells
      }


      /**
       * Some dimension dependent stuff for mesh smoothing.
       *
       * At present, this function does nothing in 1d and 2D, but
       * makes sure no two cells with a level difference greater than
       * one share one line in 3D. This is a requirement needed for
       * the interpolation of hanging nodes, since otherwise to steps
       * of interpolation would be necessary. This would make the
       * processes implemented in the @p ConstraintMatrix class much
       * more complex, since these two steps of interpolation do not
       * commute.
       */
      template <int dim, int spacedim>
      static
      void
      prepare_refinement_dim_dependent (const Triangulation<dim,spacedim> &)
      {
        Assert (dim < 3,
                ExcMessage ("Wrong function called -- there should "
                            "be a specialization."));
      }


      template <int spacedim>
      static
      void
      prepare_refinement_dim_dependent (Triangulation<3,spacedim> &triangulation)
      {
        const unsigned int dim = 3;

        // first clear flags on lines, since we need them to determine
        // which lines will be refined
        triangulation.clear_user_flags_line();

        // also clear flags on hexes, since we need them to mark those
        // cells which are to be coarsened
        triangulation.clear_user_flags_hex();

        // variable to store whether the mesh was changed in the
        // present loop and in the whole process
        bool mesh_changed      = false;

        do
          {
            mesh_changed = false;

            // for this following, we need to know which cells are
            // going to be coarsened, if we had to make a
            // decision. the following function sets these flags:
            triangulation.fix_coarsen_flags ();


            // flag those lines that are refined and will not be
            // coarsened and those that will be refined
            for (typename Triangulation<dim,spacedim>::cell_iterator
                 cell=triangulation.begin(); cell!=triangulation.end(); ++cell)
              if (cell->refine_flag_set())
                {
                  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                    if (GeometryInfo<dim>::line_refinement_case(cell->refine_flag_set(), line)
                        ==RefinementCase<1>::cut_x)
                      // flag a line, that will be
                      // refined
                      cell->line(line)->set_user_flag();
                }
              else if (cell->has_children() && !cell->child(0)->coarsen_flag_set())
                {
                  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                    if (GeometryInfo<dim>::line_refinement_case(cell->refinement_case(), line)
                        ==RefinementCase<1>::cut_x)
                      // flag a line, that is refined
                      // and will stay so
                      cell->line(line)->set_user_flag();
                }
              else if (cell->has_children() && cell->child(0)->coarsen_flag_set())
                cell->set_user_flag();


            // now check whether there are cells with lines that are
            // more than once refined or that will be more than once
            // refined. The first thing should never be the case, in
            // the second case we flag the cell for refinement
            for (typename Triangulation<dim,spacedim>::active_cell_iterator
                 cell=triangulation.last_active(); cell!=triangulation.end(); --cell)
              for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                {
                  if (cell->line(line)->has_children())
                    {
                      // if this line is refined, its children should
                      // not have further children
                      //
                      // however, if any of the children is flagged
                      // for further refinement, we need to refine
                      // this cell also (at least, if the cell is not
                      // already flagged)
                      bool offending_line_found = false;

                      for (unsigned int c=0; c<2; ++c)
                        {
                          Assert (cell->line(line)->child(c)->has_children() == false,
                                  ExcInternalError());

                          if (cell->line(line)->child(c)->user_flag_set () &&
                              (GeometryInfo<dim>::line_refinement_case(cell->refine_flag_set(),
                                                                       line)
                               ==RefinementCase<1>::no_refinement))
                            {
                              // tag this cell for refinement
                              cell->clear_coarsen_flag ();
                              // if anisotropic coarsening is allowed:
                              // extend the refine_flag in the needed
                              // direction, else set refine_flag
                              // (isotropic)
                              if (triangulation.smooth_grid &
                                  Triangulation<dim,spacedim>::allow_anisotropic_smoothing)
                                cell->flag_for_line_refinement(line);
                              else
                                cell->set_refine_flag();

                              for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
                                if (GeometryInfo<dim>::line_refinement_case(cell->refine_flag_set(), line)
                                    ==RefinementCase<1>::cut_x)
                                  // flag a line, that will be refined
                                  cell->line(l)->set_user_flag();

                              // note that we have changed the grid
                              offending_line_found = true;

                              // it may save us several loop
                              // iterations if we flag all lines of
                              // this cell now (and not at the outset
                              // of the next iteration) for refinement
                              for (unsigned int l=0;
                                   l<GeometryInfo<dim>::lines_per_cell; ++l)
                                if (!cell->line(l)->has_children() &&
                                    (GeometryInfo<dim>::line_refinement_case(cell->refine_flag_set(),
                                                                             l)
                                     !=RefinementCase<1>::no_refinement))
                                  cell->line(l)->set_user_flag();

                              break;
                            }
                        }

                      if (offending_line_found)
                        {
                          mesh_changed = true;
                          break;
                        }
                    }
                }


            // there is another thing here: if any of the lines will
            // be refined, then we may not coarsen the present cell
            // similarly, if any of the lines *is* already refined, we
            // may not coarsen the current cell. however, there's a
            // catch: if the line is refined, but the cell behind it
            // is going to be coarsened, then the situation
            // changes. if we forget this second condition, the
            // refine_and_coarsen_3d test will start to fail. note
            // that to know which cells are going to be coarsened, the
            // call for fix_coarsen_flags above is necessary
            for (typename Triangulation<dim,spacedim>::cell_iterator
                 cell=triangulation.last(); cell!=triangulation.end(); --cell)
              {
                if (cell->user_flag_set())
                  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
                    if (cell->line(line)->has_children() &&
                        (cell->line(line)->child(0)->user_flag_set() ||
                         cell->line(line)->child(1)->user_flag_set()))
                      {
                        for (unsigned int c=0; c<cell->n_children(); ++c)
                          cell->child(c)->clear_coarsen_flag ();
                        cell->clear_user_flag();
                        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
                          if (GeometryInfo<dim>::line_refinement_case(cell->refinement_case(), l)
                              ==RefinementCase<1>::cut_x)
                            // flag a line, that is refined
                            // and will stay so
                            cell->line(l)->set_user_flag();
                        mesh_changed = true;
                        break;
                      }
              }
          }
        while (mesh_changed == true);
      }



      /**
       * Helper function for @p fix_coarsen_flags. Return whether
       * coarsening of this cell is allowed.  Coarsening can be
       * forbidden if the neighboring cells are or will be refined
       * twice along the common face.
       */
      template <int dim, int spacedim>
      static
      bool
      coarsening_allowed (const typename Triangulation<dim,spacedim>::cell_iterator &cell)
      {
        // in 1d, coarsening is always allowed since we don't enforce
        // the 2:1 constraint there
        if (dim == 1)
          return true;

        const RefinementCase<dim> ref_case = cell->refinement_case();
        for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
          {

            // if the cell is not refined along that face, coarsening
            // will not change anything, so do nothing. the same
            // applies, if the face is at the boandary
            const RefinementCase<dim-1> face_ref_case =
              GeometryInfo<dim>::face_refinement_case(cell->refinement_case(), n);

            const unsigned int n_subfaces
              = GeometryInfo<dim-1>::n_children(face_ref_case);

            if (n_subfaces == 0 || cell->at_boundary(n))
              continue;
            for (unsigned int c=0; c<n_subfaces; ++c)
              {
                const typename Triangulation<dim,spacedim>::cell_iterator
                child = cell->child(GeometryInfo<dim>::
                                    child_cell_on_face(ref_case,
                                                       n,c));

                const typename Triangulation<dim,spacedim>::cell_iterator
                child_neighbor = child->neighbor(n);
                if (!child->neighbor_is_coarser(n))
                  // in 2d, if the child's neighbor is coarser, then
                  // it has no children. however, in 3d it might be
                  // otherwise. consider for example, that our face
                  // might be refined with cut_x, but the neighbor is
                  // refined with cut_xy at that face. then the
                  // neighbor pointers of the children of our cell
                  // will point to the common neighbor cell, not to
                  // its children. what we really want to know in the
                  // following is, whether the neighbor cell is
                  // refined twice with reference to our cell.  that
                  // only has to be asked, if the child's neighbor is
                  // not a coarser one.
                  if ((child_neighbor->has_children() &&
                       !child_neighbor->user_flag_set())||
                      // neighbor has children, which are further
                      // refined along the face, otherwise something
                      // went wrong in the construction of neighbor
                      // pointers.  then only allow coarsening if this
                      // neighbor will be coarsened as well
                      // (user_pointer is set).  the same applies, if
                      // the neighbors children are not refined but
                      // will be after refinement
                      child_neighbor->refine_flag_set())
                    return false;
              }
          }
        return true;
      }
    };
  }
}


template <int dim, int spacedim>
const StraightBoundary<dim,spacedim>
Triangulation<dim, spacedim>::straight_boundary = StraightBoundary<dim,spacedim>();



template <int dim, int spacedim>
const unsigned int
Triangulation<dim, spacedim>::dimension;



template <int dim, int spacedim>
Triangulation<dim, spacedim>::
Triangulation (const MeshSmoothing smooth_grid,
               const bool check_for_distorted_cells)
  :
  smooth_grid(smooth_grid),
  faces(NULL),
  anisotropic_refinement(false),
  check_for_distorted_cells(check_for_distorted_cells),
  vertex_to_boundary_id_map_1d (0),
  vertex_to_manifold_id_map_1d (0)
{
  if (dim == 1)
    {
      vertex_to_boundary_id_map_1d
        = new std::map<unsigned int, types::boundary_id>();
      vertex_to_manifold_id_map_1d
        = new std::map<unsigned int, types::manifold_id>();
    }

  // connect the any_change signal to the other signals
  signals.create.connect (signals.any_change);
  signals.post_refinement.connect (signals.any_change);
  signals.clear.connect (signals.any_change);
}


template <int dim, int spacedim>
Triangulation<dim, spacedim>::
Triangulation (const Triangulation<dim, spacedim> &other)
// do not set any subscriptors;
// anyway, calling this constructor
// is an error!
  :
  Subscriptor(),
  check_for_distorted_cells(other.check_for_distorted_cells),
  vertex_to_boundary_id_map_1d (0),
  vertex_to_manifold_id_map_1d (0)
{
  Assert (false, ExcMessage ("You are not allowed to call this constructor "
                             "because copying Triangulation objects is not "
                             "allowed. Use Triangulation::copy_from() instead."));
}



template <int dim, int spacedim>
Triangulation<dim, spacedim>::~Triangulation ()
{
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.clear ();
  delete faces;

  // the vertex_to_boundary_id_map_1d field
  // should be unused except in 1d
  Assert ((dim == 1)
          ||
          (vertex_to_boundary_id_map_1d == 0),
          ExcInternalError());
  delete vertex_to_boundary_id_map_1d;
  // the vertex_to_manifold_id_map_1d field
  // should be unused except in 1d
  Assert ((dim == 1)
          ||
          (vertex_to_manifold_id_map_1d == 0),
          ExcInternalError());
  delete vertex_to_manifold_id_map_1d;
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::clear ()
{
  clear_despite_subscriptions();
  signals.clear();
}



template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::set_mesh_smoothing(const MeshSmoothing mesh_smoothing)
{
  Assert (n_levels() == 0, ExcTriangulationNotEmpty ());
  smooth_grid=mesh_smoothing;
}



template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::set_boundary (const types::manifold_id m_number,
                                            const Boundary<dim, spacedim> &boundary_object)
{
  set_manifold(m_number, boundary_object);
}

template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::set_manifold (const types::manifold_id m_number,
                                            const Manifold<dim, spacedim> &manifold_object)
{
  Assert(m_number < numbers::invalid_manifold_id,
         ExcIndexRange(m_number,0,numbers::invalid_manifold_id));

  manifold[m_number] = &manifold_object;
}


template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::set_boundary (const types::manifold_id m_number)
{
  set_manifold(m_number);
}


template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::set_manifold (const types::manifold_id m_number)
{
  Assert(m_number < numbers::invalid_manifold_id,
         ExcIndexRange(m_number,0,numbers::invalid_manifold_id));

  //delete the entry located at number.
  manifold.erase(m_number);
}



template <int dim, int spacedim>
const Boundary<dim,spacedim> &
Triangulation<dim, spacedim>::get_boundary (const types::manifold_id m_number) const
{
  const Boundary<dim, spacedim> *man =
    dynamic_cast<const Boundary<dim, spacedim> *>(&get_manifold(m_number));
  Assert(man != NULL,
         ExcMessage("You tried to get a Boundary, but I only have a Manifold."));

  return *man;
}


template <int dim, int spacedim>
const Manifold<dim,spacedim> &
Triangulation<dim, spacedim>::get_manifold (const types::manifold_id m_number) const
{
  //look, if there is a manifold stored at
  //manifold_id number.
  typename std::map<types::manifold_id, SmartPointer<const Manifold<dim,spacedim>, Triangulation<dim, spacedim> > >::const_iterator it
    = manifold.find(m_number);

  if (it != manifold.end())
    {
      //if we have found an entry, return it
      return *(it->second);
    }
  else
    {
      //if we have not found an entry connected with number, we return
      //straight_boundary
      return straight_boundary;
    }
}




template <int dim, int spacedim>
std::vector<types::boundary_id>
Triangulation<dim, spacedim>::get_boundary_indicators () const
{
  // in 1d, we store a map of all used boundary indicators. use it for
  // our purposes
  if (dim == 1)
    {
      std::vector<types::boundary_id> boundary_indicators;
      for (std::map<unsigned int, types::boundary_id>::const_iterator
           p = vertex_to_boundary_id_map_1d->begin();
           p !=  vertex_to_boundary_id_map_1d->end();
           ++p)
        boundary_indicators.push_back (p->second);

      return boundary_indicators;
    }
  else
    {
      std::set<types::boundary_id> b_ids;
      active_cell_iterator cell=begin_active();
      for (; cell!=end(); ++cell)
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->at_boundary(face))
            b_ids.insert(cell->face(face)->boundary_indicator());
      std::vector<types::boundary_id> boundary_indicators(b_ids.begin(), b_ids.end());
      return boundary_indicators;
    }
}

template <int dim, int spacedim>
std::vector<types::manifold_id>
Triangulation<dim, spacedim>::get_manifold_ids () const
{
  std::set<types::manifold_id> m_ids;
  active_cell_iterator cell=begin_active();
  for (; cell!=end(); ++cell)
    {
      m_ids.insert(cell->manifold_id());
      if (dim>1)
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->at_boundary(face))
            m_ids.insert(cell->face(face)->manifold_id());
    }
  std::vector<types::manifold_id> manifold_indicators(m_ids.begin(), m_ids.end());
  return manifold_indicators;
}

/*-----------------------------------------------------------------*/


template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::
copy_triangulation (const Triangulation<dim, spacedim> &old_tria)
{
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size () == 0, ExcTriangulationNotEmpty());
  Assert (faces == NULL, ExcTriangulationNotEmpty());

  Assert (old_tria.levels.size() != 0, ExcInternalError());
  Assert (old_tria.vertices.size() != 0, ExcInternalError());
  Assert (dim == 1 || old_tria.faces != NULL, ExcInternalError());


  // copy normal elements
  vertices               = old_tria.vertices;
  vertices_used          = old_tria.vertices_used;
  anisotropic_refinement = old_tria.anisotropic_refinement;
  smooth_grid            = old_tria.smooth_grid;

  faces         = new internal::Triangulation::TriaFaces<dim>(*old_tria.faces);

  typename std::map<types::manifold_id,
           SmartPointer<const Manifold<dim,spacedim> , Triangulation<dim, spacedim> > >::const_iterator
           bdry_iterator = old_tria.manifold.begin();
  for (; bdry_iterator != old_tria.manifold.end() ; bdry_iterator++)
    manifold[bdry_iterator->first] = bdry_iterator->second;


  levels.reserve (old_tria.levels.size());
  for (unsigned int level=0; level<old_tria.levels.size(); ++level)
    levels.push_back (new
                      internal::Triangulation::
                      TriaLevel<dim>(*old_tria.levels[level]));

  number_cache = old_tria.number_cache;

  if (dim == 1)
    {
      delete vertex_to_boundary_id_map_1d;
      vertex_to_boundary_id_map_1d
        = (new std::map<unsigned int, types::boundary_id>
           (*old_tria.vertex_to_boundary_id_map_1d));

      delete vertex_to_manifold_id_map_1d;
      vertex_to_manifold_id_map_1d
        = (new std::map<unsigned int, types::manifold_id>
           (*old_tria.vertex_to_manifold_id_map_1d));
    }

  // inform those who are listening on old_tria of the copy operation
  old_tria.signals.copy (*this);
  // also inform all listeners of the current triangulation that the
  // triangulation has been created
  signals.create();

  // note that we need not copy the
  // subscriptor!
}



template <int dim, int spacedim>
void
Triangulation<dim,spacedim>::
create_triangulation_compatibility (const std::vector<Point<spacedim> > &v,
                                    const std::vector<CellData<dim> >   &cells,
                                    const SubCellData                   &subcelldata)
{
  std::vector<CellData<dim> > reordered_cells (cells);
  SubCellData                 reordered_subcelldata (subcelldata);

  // in-place reordering of data
  reorder_compatibility (reordered_cells, reordered_subcelldata);

  // now create triangulation from
  // reordered data
  create_triangulation(v, reordered_cells, reordered_subcelldata);
}


template <int dim, int spacedim>
void
Triangulation<dim,spacedim>::
create_triangulation (const std::vector<Point<spacedim> >    &v,
                      const std::vector<CellData<dim> > &cells,
                      const SubCellData &subcelldata)
{
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
  Assert (faces == NULL, ExcTriangulationNotEmpty());
  // check that no forbidden arrays
  // are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  // try to create a triangulation; if this fails, we still want to
  // throw an exception but if we just do so we'll get into trouble
  // because sometimes other objects are already attached to it:
  try
    {
      internal::Triangulation::Implementation::create_triangulation (v, cells, subcelldata, *this);
    }
  catch (...)
    {
      clear_despite_subscriptions();
      throw;
    }

  internal::Triangulation::Implementation
  ::compute_number_cache (*this, levels.size(), number_cache);

  // now verify that there are indeed no distorted cells. as per the
  // documentation of this class, we first collect all distorted cells
  // and then throw an exception if there are any
  if (check_for_distorted_cells == true)
    {
      DistortedCellList distorted_cells = collect_distorted_coarse_cells (*this);
      // throw the array (and fill the various location fields) if
      // there are distorted cells. otherwise, just fall off the end
      // of the function
      AssertThrow (distorted_cells.distorted_cells.size() == 0,
                   distorted_cells);
    }


  /*
      When the triangulation is a manifold (dim < spacedim), the normal field
      provided from the map class depends on the order of the vertices.
      It may happen that this normal field is discontinous.
      The following code takes care that this is not the case by setting the
      cell direction flag on those cell that produce the wrong orientation.

      To determine if 2 neighbours have the same or opposite orientation
      we use a table of truth.
      Its entries are indexes by the local indeces of the common face.
      For example if two elements share a face, and this face is
      face 0 for element 0 and face 1 for element 1, then
      table(0,1) will tell whether the orientation are the same (true) or
      opposite (false).

      Even though there may be a combinatorial/graph theory argument to get
      this table in any dimension, I tested by hand all the different possible
      cases in 1D and 2D to generate the table.

      Assuming that a surface respects the standard orientation for 2d meshes,
      the tables of truth are symmetric and their true values are the following
      1D curves:  (0,1)
      2D surface: (0,1),(0,2),(1,3),(2,3)

      We store this data using an n_faces x n_faces full matrix, which is actually
      much bigger than the minimal data required, but it makes the code more readable.

    */
  if (dim < spacedim)
    {
      Table<2,bool> correct(GeometryInfo< dim >::faces_per_cell,
                            GeometryInfo< dim >::faces_per_cell);
      switch (dim)
        {
        case 1:
        {
          bool values [][2] = {{false,true},
            {true,false}
          };
          for (unsigned int i=0; i< GeometryInfo< dim >::faces_per_cell; ++i)
            for (unsigned int j=0; j< GeometryInfo< dim >::faces_per_cell; ++j)
              correct(i,j) = ( values[i][j]);
          break;
        }
        case 2:
        {
          bool values [][4]= {{false,true ,true , false},
            {true ,false,false, true },
            {true ,false,false, true },
            {false,true ,true , false}
          };
          for (unsigned int i=0; i< GeometryInfo< dim >::faces_per_cell; ++i)
            for (unsigned int j=0; j< GeometryInfo< dim >::faces_per_cell; ++j)
              correct(i,j) = ( values[i][j]);
          break;
        }
        default:
          Assert (false, ExcNotImplemented());
        }


      std::list<active_cell_iterator> this_round, next_round;
      active_cell_iterator neighbor;

      this_round.push_back (begin_active());
      begin_active()->set_direction_flag (true);
      begin_active()->set_user_flag ();

      while (this_round.size() > 0)
        {
          for ( typename std::list<active_cell_iterator>::iterator cell = this_round.begin();
                cell != this_round.end(); ++cell)
            {
              for (unsigned int i = 0; i < GeometryInfo< dim >::faces_per_cell; ++i)
                {
                  if ( !((*cell)->face(i)->at_boundary()) )
                    {
                      neighbor = (*cell)->neighbor(i);

                      unsigned int cf = (*cell)->face_index(i);
                      unsigned int j = 0;
                      while (neighbor->face_index(j) != cf)
                        {
                          ++j;
                        }

                      if ( (correct(i,j) && !(*cell)->direction_flag())
                           ||
                           (!correct(i,j) && (*cell)->direction_flag()) )
                        {
                          if (neighbor->user_flag_set() == false)
                            {
                              neighbor->set_direction_flag (false);
                              neighbor->set_user_flag ();
                              next_round.push_back (neighbor);
                            }
                          else
                            Assert (neighbor->direction_flag() == false,
                                    ExcNonOrientableTriangulation());

                        }
                    }
                }
            }

          // Before we quit let's check
          // that if the triangulation
          // is disconnected that we
          // still get all cells
          if (next_round.size() == 0)
            for (active_cell_iterator cell = begin_active();
                 cell != end(); ++cell)
              if (cell->user_flag_set() == false)
                {
                  next_round.push_back (cell);
                  cell->set_direction_flag (true);
                  cell->set_user_flag ();
                  break;
                }

          this_round = next_round;
          next_round.clear();
        }
    }

  // inform all listeners that the triangulation has been created
  signals.create();
}




template <int dim, int spacedim>
void
Triangulation<dim,spacedim>::
flip_all_direction_flags()
{
  AssertThrow (dim+1 == spacedim, ExcMessage ("Only works for dim == spacedim-1"));
  for (active_cell_iterator cell = begin_active();
       cell != end(); ++cell)
    cell->set_direction_flag (!cell->direction_flag());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::distort_random (const double factor,
                                                   const bool   keep_boundary)
{
  internal::Triangulation::Implementation::distort_random (factor, keep_boundary, *this);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::set_all_refine_flags ()
{
  Assert(n_cells()>0, ExcMessage("Error: An empty Triangulation can not be refined."));
  active_cell_iterator cell = begin_active(),
                       endc = end();

  for (; cell != endc; ++cell)
    {
      cell->clear_coarsen_flag();
      cell->set_refine_flag ();
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::refine_global (const unsigned int times)
{
  for (unsigned int i=0; i<times; ++i)
    {
      set_all_refine_flags();
      execute_coarsening_and_refinement ();
    }
}



/*-------------------- refine/coarsen flags -------------------------*/



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_refine_flags (std::vector<bool> &v) const
{
  v.resize (dim*n_active_cells(), false);
  std::vector<bool>::iterator  i = v.begin();
  active_cell_iterator cell = begin_active(),
                       endc = end();
  for (; cell!=endc; ++cell)
    for (unsigned int j=0; j<dim; ++j,++i)
      if (cell->refine_flag_set() & (1<<j) )
        *i = true;

  Assert (i == v.end(), ExcInternalError());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_refine_flags (std::ostream &out) const
{
  std::vector<bool> v;
  save_refine_flags (v);
  write_bool_vector (mn_tria_refine_flags_begin, v, mn_tria_refine_flags_end,
                     out);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_refine_flags (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_refine_flags_begin, v, mn_tria_refine_flags_end,
                    in);
  load_refine_flags (v);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_refine_flags (const std::vector<bool> &v)
{
  AssertThrow (v.size() == dim*n_active_cells(), ExcGridReadError());

  active_cell_iterator cell = begin_active(),
                       endc = end();
  std::vector<bool>::const_iterator i = v.begin();
  for (; cell!=endc; ++cell)
    {
      unsigned int ref_case=0;

      for (unsigned int j=0; j<dim; ++j, ++i)
        if (*i == true)
          ref_case+=1<<j;
      Assert(ref_case<RefinementCase<dim>::isotropic_refinement+1,
             ExcGridReadError());
      if (ref_case>0)
        cell->set_refine_flag(RefinementCase<dim>(ref_case));
      else
        cell->clear_refine_flag();
    }

  Assert (i == v.end(), ExcInternalError());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_coarsen_flags (std::vector<bool> &v) const
{
  v.resize (n_active_cells(), false);
  std::vector<bool>::iterator  i = v.begin();
  active_cell_iterator cell = begin_active(),
                       endc = end();
  for (; cell!=endc; ++cell, ++i)
    *i = cell->coarsen_flag_set();

  Assert (i == v.end(), ExcInternalError());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_coarsen_flags (std::ostream &out) const
{
  std::vector<bool> v;
  save_coarsen_flags (v);
  write_bool_vector (mn_tria_coarsen_flags_begin, v, mn_tria_coarsen_flags_end,
                     out);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_coarsen_flags (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_coarsen_flags_begin, v, mn_tria_coarsen_flags_end,
                    in);
  load_coarsen_flags (v);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_coarsen_flags (const std::vector<bool> &v)
{
  Assert (v.size() == n_active_cells(), ExcGridReadError());

  active_cell_iterator cell = begin_active(),
                       endc = end();
  std::vector<bool>::const_iterator i = v.begin();
  for (; cell!=endc; ++cell, ++i)
    if (*i == true)
      cell->set_coarsen_flag();
    else
      cell->clear_coarsen_flag();

  Assert (i == v.end(), ExcInternalError());
}


template <int dim, int spacedim>
bool Triangulation<dim,spacedim>::get_anisotropic_refinement_flag() const
{
  return anisotropic_refinement;
}



/*-------------------- user data/flags -------------------------*/


namespace
{
  // clear user data of cells
  template <int dim>
  void clear_user_data (std::vector<internal::Triangulation::TriaLevel<dim>*> &levels)
  {
    for (unsigned int level=0; level<levels.size(); ++level)
      levels[level]->cells.clear_user_data();
  }


  // clear user data of faces
  void clear_user_data (internal::Triangulation::TriaFaces<1> *)
  {
    // nothing to do in 1d
  }


  void clear_user_data (internal::Triangulation::TriaFaces<2> *faces)
  {
    faces->lines.clear_user_data();
  }


  void clear_user_data (internal::Triangulation::TriaFaces<3> *faces)
  {
    faces->lines.clear_user_data();
    faces->quads.clear_user_data();
  }
}


template <int dim, int spacedim>
void Triangulation<dim,spacedim>::clear_user_data ()
{
  // let functions in anonymous namespace do their work
  dealii::clear_user_data (levels);
  dealii::clear_user_data (faces);
}



namespace
{
  void clear_user_flags_line (std::vector<internal::Triangulation::TriaLevel<1>*> &levels,
                              internal::Triangulation::TriaFaces<1> *)
  {
    for (unsigned int level=0; level<levels.size(); ++level)
      levels[level]->cells.clear_user_flags();
  }

  template <int dim>
  void clear_user_flags_line (std::vector<internal::Triangulation::TriaLevel<dim>*> &,
                              internal::Triangulation::TriaFaces<dim> *faces)
  {
    faces->lines.clear_user_flags();
  }
}


template <int dim, int spacedim>
void Triangulation<dim,spacedim>::clear_user_flags_line ()
{
  dealii::clear_user_flags_line (levels, faces);
}



namespace
{
  void clear_user_flags_quad (std::vector<internal::Triangulation::TriaLevel<1>*> &,
                              internal::Triangulation::TriaFaces<1> *)
  {
    // nothing to do in 1d
  }

  void clear_user_flags_quad (std::vector<internal::Triangulation::TriaLevel<2>*> &levels,
                              internal::Triangulation::TriaFaces<2> *)
  {
    for (unsigned int level=0; level<levels.size(); ++level)
      levels[level]->cells.clear_user_flags();
  }

  template <int dim>
  void clear_user_flags_quad (std::vector<internal::Triangulation::TriaLevel<dim>*> &,
                              internal::Triangulation::TriaFaces<dim> *faces)
  {
    faces->quads.clear_user_flags();
  }
}


template <int dim, int spacedim>
void Triangulation<dim,spacedim>::clear_user_flags_quad ()
{
  dealii::clear_user_flags_quad (levels, faces);
}



namespace
{
  void clear_user_flags_hex (std::vector<internal::Triangulation::TriaLevel<1>*> &,
                             internal::Triangulation::TriaFaces<1> *)
  {
    // nothing to do in 1d
  }


  void clear_user_flags_hex (std::vector<internal::Triangulation::TriaLevel<2>*> &,
                             internal::Triangulation::TriaFaces<2> *)
  {
    // nothing to do in 2d
  }

  void clear_user_flags_hex (std::vector<internal::Triangulation::TriaLevel<3>*> &levels,
                             internal::Triangulation::TriaFaces<3> *)
  {
    for (unsigned int level=0; level<levels.size(); ++level)
      levels[level]->cells.clear_user_flags();
  }
}


template <int dim, int spacedim>
void Triangulation<dim,spacedim>::clear_user_flags_hex ()
{
  dealii::clear_user_flags_hex (levels, faces);
}



template <int dim, int spacedim>
void Triangulation<dim,spacedim>::clear_user_flags ()
{
  clear_user_flags_line ();
  clear_user_flags_quad ();
  clear_user_flags_hex ();
}



template <int dim, int spacedim>
void Triangulation<dim,spacedim>::clear_user_pointers ()
{
  clear_user_data();
}




template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags (std::ostream &out) const
{
  save_user_flags_line (out);

  if (dim>=2)
    save_user_flags_quad (out);

  if (dim>=3)
    save_user_flags_hex (out);

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags (std::vector<bool> &v) const
{
  // clear vector and append
  // all the stuff later on
  v.clear ();

  std::vector<bool> tmp;

  save_user_flags_line (tmp);
  v.insert (v.end(), tmp.begin(), tmp.end());

  if (dim >= 2)
    {
      save_user_flags_quad (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }

  if (dim >= 3)
    {
      save_user_flags_hex (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags (std::istream &in)
{
  load_user_flags_line (in);

  if (dim>=2)
    load_user_flags_quad (in);

  if (dim>=3)
    load_user_flags_hex (in);

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags (const std::vector<bool> &v)
{
  Assert (v.size() == n_lines()+n_quads()+n_hexs(), ExcInternalError());
  std::vector<bool> tmp;

  // first extract the flags
  // belonging to lines
  tmp.insert (tmp.end(),
              v.begin(), v.begin()+n_lines());
  // and set the lines
  load_user_flags_line (tmp);

  if (dim >= 2)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
                  v.begin()+n_lines(), v.begin()+n_lines()+n_quads());
      load_user_flags_quad (tmp);
    }

  if (dim >= 3)
    {
      tmp.clear();
      tmp.insert (tmp.end(),
                  v.begin()+n_lines()+n_quads(), v.begin()+n_lines()+n_quads()+n_hexs());
      load_user_flags_hex (tmp);
    }

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags_line (std::vector<bool> &v) const
{
  v.resize (n_lines(), false);
  std::vector<bool>::iterator  i = v.begin();
  line_iterator line = begin_line(),
                endl = end_line();
  for (; line!=endl; ++line, ++i)
    *i = line->user_flag_set();

  Assert (i == v.end(), ExcInternalError());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags_line (std::ostream &out) const
{
  std::vector<bool> v;
  save_user_flags_line (v);
  write_bool_vector (mn_tria_line_user_flags_begin, v, mn_tria_line_user_flags_end,
                     out);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags_line (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_line_user_flags_begin, v, mn_tria_line_user_flags_end,
                    in);
  load_user_flags_line (v);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags_line (const std::vector<bool> &v)
{
  Assert (v.size() == n_lines(), ExcGridReadError());

  line_iterator line = begin_line(),
                endl = end_line();
  std::vector<bool>::const_iterator i = v.begin();
  for (; line!=endl; ++line, ++i)
    if (*i == true)
      line->set_user_flag();
    else
      line->clear_user_flag();

  Assert (i == v.end(), ExcInternalError());
}


namespace
{
  template <typename Iterator>
  bool get_user_flag (const Iterator &i)
  {
    return i->user_flag_set();
  }



  template <int structdim, int dim, int spacedim>
  bool get_user_flag (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &)
  {
    Assert (false, ExcInternalError());
    return false;
  }



  template <typename Iterator>
  void set_user_flag (const Iterator &i)
  {
    i->set_user_flag();
  }



  template <int structdim, int dim, int spacedim>
  void set_user_flag (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &)
  {
    Assert (false, ExcInternalError());
  }



  template <typename Iterator>
  void clear_user_flag (const Iterator &i)
  {
    i->clear_user_flag();
  }



  template <int structdim, int dim, int spacedim>
  void clear_user_flag (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &)
  {
    Assert (false, ExcInternalError());
  }
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags_quad (std::vector<bool> &v) const
{
  v.resize (n_quads(), false);

  if (dim >= 2)
    {
      std::vector<bool>::iterator  i = v.begin();
      quad_iterator quad = begin_quad(),
                    endq = end_quad();
      for (; quad!=endq; ++quad, ++i)
        *i = get_user_flag (quad);

      Assert (i == v.end(), ExcInternalError());
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags_quad (std::ostream &out) const
{
  std::vector<bool> v;
  save_user_flags_quad (v);
  write_bool_vector (mn_tria_quad_user_flags_begin, v, mn_tria_quad_user_flags_end,
                     out);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags_quad (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_quad_user_flags_begin, v, mn_tria_quad_user_flags_end,
                    in);
  load_user_flags_quad (v);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags_quad (const std::vector<bool> &v)
{
  Assert (v.size() == n_quads(), ExcGridReadError());

  if (dim >= 2)
    {
      quad_iterator quad = begin_quad(),
                    endq = end_quad();
      std::vector<bool>::const_iterator i = v.begin();
      for (; quad!=endq; ++quad, ++i)
        if (*i == true)
          set_user_flag(quad);
        else
          clear_user_flag(quad);

      Assert (i == v.end(), ExcInternalError());
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags_hex (std::vector<bool> &v) const
{
  v.resize (n_hexs(), false);

  if (dim >= 3)
    {
      std::vector<bool>::iterator  i = v.begin();
      hex_iterator hex = begin_hex(),
                   endh = end_hex();
      for (; hex!=endh; ++hex, ++i)
        *i = get_user_flag (hex);

      Assert (i == v.end(), ExcInternalError());
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_flags_hex (std::ostream &out) const
{
  std::vector<bool> v;
  save_user_flags_hex (v);
  write_bool_vector (mn_tria_hex_user_flags_begin, v, mn_tria_hex_user_flags_end,
                     out);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags_hex (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_hex_user_flags_begin, v, mn_tria_hex_user_flags_end,
                    in);
  load_user_flags_hex (v);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_flags_hex (const std::vector<bool> &v)
{
  Assert (v.size() == n_hexs(), ExcGridReadError());

  if (dim >= 3)
    {
      hex_iterator hex = begin_hex(),
                   endh = end_hex();
      std::vector<bool>::const_iterator i = v.begin();
      for (; hex!=endh; ++hex, ++i)
        if (*i == true)
          set_user_flag(hex);
        else
          clear_user_flag(hex);

      Assert (i == v.end(), ExcInternalError());
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_indices (std::vector<unsigned int> &v) const
{
  // clear vector and append all the
  // stuff later on
  v.clear ();

  std::vector<unsigned int> tmp;

  save_user_indices_line (tmp);
  v.insert (v.end(), tmp.begin(), tmp.end());

  if (dim >= 2)
    {
      save_user_indices_quad (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }

  if (dim >= 3)
    {
      save_user_indices_hex (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_indices (const std::vector<unsigned int> &v)
{
  Assert (v.size() == n_lines()+n_quads()+n_hexs(), ExcInternalError());
  std::vector<unsigned int> tmp;

  // first extract the indices
  // belonging to lines
  tmp.insert (tmp.end(),
              v.begin(), v.begin()+n_lines());
  // and set the lines
  load_user_indices_line (tmp);

  if (dim >= 2)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
                  v.begin()+n_lines(), v.begin()+n_lines()+n_quads());
      load_user_indices_quad (tmp);
    }

  if (dim >= 3)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
                  v.begin()+n_lines()+n_quads(), v.begin()+n_lines()+n_quads()+n_hexs());
      load_user_indices_hex (tmp);
    }

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



namespace
{
  template <typename Iterator>
  unsigned int get_user_index (const Iterator &i)
  {
    return i->user_index();
  }



  template <int structdim, int dim, int spacedim>
  unsigned int get_user_index (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &)
  {
    Assert (false, ExcInternalError());
    return numbers::invalid_unsigned_int;
  }



  template <typename Iterator>
  void set_user_index (const Iterator &i,
                       const unsigned int x)
  {
    i->set_user_index(x);
  }



  template <int structdim, int dim, int spacedim>
  void set_user_index (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &,
                       const unsigned int)
  {
    Assert (false, ExcInternalError());
  }
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_indices_line (std::vector<unsigned int> &v) const
{
  v.resize (n_lines(), 0);
  std::vector<unsigned int>::iterator  i = v.begin();
  line_iterator line = begin_line(),
                endl = end_line();
  for (; line!=endl; ++line, ++i)
    *i = line->user_index();
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_indices_line (const std::vector<unsigned int> &v)
{
  Assert (v.size() == n_lines(), ExcGridReadError());

  line_iterator line = begin_line(),
                endl = end_line();
  std::vector<unsigned int>::const_iterator i = v.begin();
  for (; line!=endl; ++line, ++i)
    line->set_user_index(*i);
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_indices_quad (std::vector<unsigned int> &v) const
{
  v.resize (n_quads(), 0);

  if (dim >= 2)
    {
      std::vector<unsigned int>::iterator  i = v.begin();
      quad_iterator quad = begin_quad(),
                    endq = end_quad();
      for (; quad!=endq; ++quad, ++i)
        *i = get_user_index(quad);
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_indices_quad (const std::vector<unsigned int> &v)
{
  Assert (v.size() == n_quads(), ExcGridReadError());

  if (dim >= 2)
    {
      quad_iterator quad = begin_quad(),
                    endq = end_quad();
      std::vector<unsigned int>::const_iterator i = v.begin();
      for (; quad!=endq; ++quad, ++i)
        set_user_index(quad, *i);
    }
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_indices_hex (std::vector<unsigned int> &v) const
{
  v.resize (n_hexs(), 0);

  if (dim >= 3)
    {
      std::vector<unsigned int>::iterator  i = v.begin();
      hex_iterator hex = begin_hex(),
                   endh = end_hex();
      for (; hex!=endh; ++hex, ++i)
        *i = get_user_index(hex);
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_indices_hex (const std::vector<unsigned int> &v)
{
  Assert (v.size() == n_hexs(), ExcGridReadError());

  if (dim >= 3)
    {
      hex_iterator hex = begin_hex(),
                   endh = end_hex();
      std::vector<unsigned int>::const_iterator i = v.begin();
      for (; hex!=endh; ++hex, ++i)
        set_user_index(hex, *i);
    }
}



//---------------- user pointers ----------------------------------------//


namespace
{
  template <typename Iterator>
  void *get_user_pointer (const Iterator &i)
  {
    return i->user_pointer();
  }



  template <int structdim, int dim, int spacedim>
  void *get_user_pointer (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &)
  {
    Assert (false, ExcInternalError());
    return 0;
  }



  template <typename Iterator>
  void set_user_pointer (const Iterator &i,
                         void *x)
  {
    i->set_user_pointer(x);
  }



  template <int structdim, int dim, int spacedim>
  void set_user_pointer (const TriaIterator<InvalidAccessor<structdim,dim,spacedim> > &,
                         void *)
  {
    Assert (false, ExcInternalError());
  }
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_pointers (std::vector<void *> &v) const
{
  // clear vector and append all the
  // stuff later on
  v.clear ();

  std::vector<void *> tmp;

  save_user_pointers_line (tmp);
  v.insert (v.end(), tmp.begin(), tmp.end());

  if (dim >= 2)
    {
      save_user_pointers_quad (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }

  if (dim >= 3)
    {
      save_user_pointers_hex (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_pointers (const std::vector<void *> &v)
{
  Assert (v.size() == n_lines()+n_quads()+n_hexs(), ExcInternalError());
  std::vector<void *> tmp;

  // first extract the pointers
  // belonging to lines
  tmp.insert (tmp.end(),
              v.begin(), v.begin()+n_lines());
  // and set the lines
  load_user_pointers_line (tmp);

  if (dim >= 2)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
                  v.begin()+n_lines(), v.begin()+n_lines()+n_quads());
      load_user_pointers_quad (tmp);
    }

  if (dim >= 3)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
                  v.begin()+n_lines()+n_quads(), v.begin()+n_lines()+n_quads()+n_hexs());
      load_user_pointers_hex (tmp);
    }

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_pointers_line (std::vector<void *> &v) const
{
  v.resize (n_lines(), 0);
  std::vector<void *>::iterator  i = v.begin();
  line_iterator line = begin_line(),
                endl = end_line();
  for (; line!=endl; ++line, ++i)
    *i = line->user_pointer();
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_pointers_line (const std::vector<void *> &v)
{
  Assert (v.size() == n_lines(), ExcGridReadError());

  line_iterator line = begin_line(),
                endl = end_line();
  std::vector<void *>::const_iterator i = v.begin();
  for (; line!=endl; ++line, ++i)
    line->set_user_pointer(*i);
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_pointers_quad (std::vector<void *> &v) const
{
  v.resize (n_quads(), 0);

  if (dim >= 2)
    {
      std::vector<void *>::iterator  i = v.begin();
      quad_iterator quad = begin_quad(),
                    endq = end_quad();
      for (; quad!=endq; ++quad, ++i)
        *i = get_user_pointer(quad);
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_pointers_quad (const std::vector<void *> &v)
{
  Assert (v.size() == n_quads(), ExcGridReadError());

  if (dim >= 2)
    {
      quad_iterator quad = begin_quad(),
                    endq = end_quad();
      std::vector<void *>::const_iterator i = v.begin();
      for (; quad!=endq; ++quad, ++i)
        set_user_pointer(quad, *i);
    }
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::save_user_pointers_hex (std::vector<void *> &v) const
{
  v.resize (n_hexs(), 0);

  if (dim >= 3)
    {
      std::vector<void *>::iterator  i = v.begin();
      hex_iterator hex = begin_hex(),
                   endh = end_hex();
      for (; hex!=endh; ++hex, ++i)
        *i = get_user_pointer(hex);
    }
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::load_user_pointers_hex (const std::vector<void *> &v)
{
  Assert (v.size() == n_hexs(), ExcGridReadError());

  if (dim >= 3)
    {
      hex_iterator hex = begin_hex(),
                   endh = end_hex();
      std::vector<void *>::const_iterator i = v.begin();
      for (; hex!=endh; ++hex, ++i)
        set_user_pointer(hex, *i);
    }
}



/*------------------------ Cell iterator functions ------------------------*/


template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::raw_cell_iterator
Triangulation<dim,spacedim>::begin_raw (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      return begin_raw_line (level);
    case 2:
      return begin_raw_quad (level);
    case 3:
      return begin_raw_hex (level);
    default:
      Assert (false, ExcNotImplemented());
      return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::cell_iterator
Triangulation<dim,spacedim>::begin (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      return begin_line (level);
    case 2:
      return begin_quad (level);
    case 3:
      return begin_hex (level);
    default:
      Assert (false, ExcImpossibleInDim(dim));
      return cell_iterator();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::active_cell_iterator
Triangulation<dim,spacedim>::begin_active (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      return begin_active_line (level);
    case 2:
      return begin_active_quad (level);
    case 3:
      return begin_active_hex (level);
    default:
      Assert (false, ExcNotImplemented());
      return active_cell_iterator();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::cell_iterator
Triangulation<dim,spacedim>::last () const
{
  const unsigned int level = levels.size()-1;

  Assert (level<n_global_levels() || level<levels.size(), ExcInvalidLevel(level));
  if (levels[level]->cells.cells.size() ==0)
    return end(level);

  // find the last raw iterator on
  // this level
  raw_cell_iterator ri (const_cast<Triangulation<dim,spacedim>*>(this),
                        level,
                        levels[level]->cells.cells.size()-1);

  // then move to the last used one
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::active_cell_iterator
Triangulation<dim,spacedim>::last_active () const
{
  // get the last used cell
  cell_iterator cell = last();

  if (cell != end())
    {
      // then move to the last active one
      if (cell->active()==true)
        return cell;
      while ((--cell).state() == IteratorState::valid)
        if (cell->active()==true)
          return cell;
    }
  return cell;
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::cell_iterator
Triangulation<dim,spacedim>::end () const
{
  return cell_iterator (const_cast<Triangulation<dim, spacedim>*>(this),
                        -1,
                        -1);
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::raw_cell_iterator
Triangulation<dim, spacedim>::end_raw (const unsigned int level) const
{
  Assert (level<n_global_levels(), ExcInvalidLevel(level));
  if (level < levels.size()-1)
    return begin_raw (level+1);
  else
    return end();
}


template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::cell_iterator
Triangulation<dim, spacedim>::end (const unsigned int level) const
{
  if (level < levels.size()-1)
    return begin (level+1);
  Assert (level<n_global_levels() || level<levels.size(), ExcInvalidLevel(level));
  return end();
}


template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::active_cell_iterator
Triangulation<dim, spacedim>::end_active (const unsigned int level) const
{
  Assert (level<n_global_levels() || level < levels.size(), ExcInvalidLevel(level));
  return (level >= levels.size()-1 ?
          active_cell_iterator(end()) :
          begin_active (level+1));
}



template <int dim, int spacedim>
IteratorRange<typename Triangulation<dim, spacedim>::cell_iterator>
Triangulation<dim, spacedim>::cell_iterators () const
{
  return
    IteratorRange<typename Triangulation<dim, spacedim>::cell_iterator>
    (begin(), end());
}


template <int dim, int spacedim>
IteratorRange<typename Triangulation<dim, spacedim>::active_cell_iterator>
Triangulation<dim, spacedim>::active_cell_iterators () const
{
  return
    IteratorRange<typename Triangulation<dim, spacedim>::active_cell_iterator>
    (begin_active(), end());
}



template <int dim, int spacedim>
IteratorRange<typename Triangulation<dim, spacedim>::cell_iterator>
Triangulation<dim, spacedim>::cell_iterators_on_level (const unsigned int level) const
{
  return
    IteratorRange<typename Triangulation<dim, spacedim>::cell_iterator>
    (begin(level), end(level));
}



template <int dim, int spacedim>
IteratorRange<typename Triangulation<dim, spacedim>::active_cell_iterator>
Triangulation<dim, spacedim>::active_cell_iterators_on_level (const unsigned int level) const
{
  return
    IteratorRange<typename Triangulation<dim, spacedim>::active_cell_iterator>
    (begin_active(level), end_active(level));
}


/*------------------------ Face iterator functions ------------------------*/


template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::face_iterator
Triangulation<dim,spacedim>::begin_face () const
{
  switch (dim)
    {
    case 1:
      Assert (false, ExcImpossibleInDim(1));
      return raw_face_iterator();
    case 2:
      return begin_line ();
    case 3:
      return begin_quad ();
    default:
      Assert (false, ExcNotImplemented());
      return face_iterator ();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::active_face_iterator
Triangulation<dim,spacedim>::begin_active_face () const
{
  switch (dim)
    {
    case 1:
      Assert (false, ExcImpossibleInDim(1));
      return raw_face_iterator();
    case 2:
      return begin_active_line ();
    case 3:
      return begin_active_quad ();
    default:
      Assert (false, ExcNotImplemented());
      return active_face_iterator ();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::face_iterator
Triangulation<dim,spacedim>::end_face () const
{
  switch (dim)
    {
    case 1:
      Assert (false, ExcImpossibleInDim(1));
      return raw_face_iterator();
    case 2:
      return end_line ();
    case 3:
      return end_quad ();
    default:
      Assert (false, ExcNotImplemented());
      return raw_face_iterator ();
    }
}




/*------------------------ Line iterator functions ------------------------*/



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::raw_line_iterator
Triangulation<dim, spacedim>::begin_raw_line (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      Assert (level<n_global_levels() || level<levels.size(), ExcInvalidLevel(level));

      if (level >= levels.size() || levels[level]->cells.cells.size() == 0)
        return end_line();

      return raw_line_iterator (const_cast<Triangulation<dim,spacedim>*>(this),
                                level,
                                0);

    default:
      Assert (level == 0, ExcFacesHaveNoLevel());
      return raw_line_iterator (const_cast<Triangulation<dim, spacedim>*>(this),
                                0,
                                0);
    }
}


template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::line_iterator
Triangulation<dim, spacedim>::begin_line (const unsigned int level) const
{
  // level is checked in begin_raw
  raw_line_iterator ri = begin_raw_line (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::active_line_iterator
Triangulation<dim, spacedim>::begin_active_line (const unsigned int level) const
{
  // level is checked in begin_raw
  line_iterator i = begin_line (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::line_iterator
Triangulation<dim, spacedim>::end_line () const
{
  return raw_line_iterator (const_cast<Triangulation<dim, spacedim>*>(this),
                            -1,
                            -1);
}



/*------------------------ Quad iterator functions ------------------------*/


template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::raw_quad_iterator
Triangulation<dim,spacedim>::begin_raw_quad (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      Assert (false, ExcImpossibleInDim(1));
      return raw_hex_iterator();
    case 2:
    {
      Assert (level<n_global_levels() || level<levels.size(), ExcInvalidLevel(level));

      if (level >= levels.size() || levels[level]->cells.cells.size() == 0)
        return end_quad();

      return raw_quad_iterator (const_cast<Triangulation<dim,spacedim>*>(this),
                                level,
                                0);
    }

    case 3:
    {
      Assert (level == 0, ExcFacesHaveNoLevel());

      return raw_quad_iterator (const_cast<Triangulation<dim,spacedim>*>(this),
                                0,
                                0);
    }


    default:
      Assert (false, ExcNotImplemented());
      return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::quad_iterator
Triangulation<dim,spacedim>::begin_quad (const unsigned int level) const
{
  // level is checked in begin_raw
  raw_quad_iterator ri = begin_raw_quad (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::active_quad_iterator
Triangulation<dim,spacedim>::begin_active_quad (const unsigned int level) const
{
  // level is checked in begin_raw
  quad_iterator i = begin_quad (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::quad_iterator
Triangulation<dim,spacedim>::end_quad () const
{
  return raw_quad_iterator (const_cast<Triangulation<dim, spacedim>*>(this),
                            -1,
                            -1);
}


/*------------------------ Hex iterator functions ------------------------*/


template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::raw_hex_iterator
Triangulation<dim,spacedim>::begin_raw_hex (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
    case 2:
      Assert (false, ExcImpossibleInDim(1));
      return raw_hex_iterator();
    case 3:
    {
      Assert (level<n_global_levels() || level<levels.size(), ExcInvalidLevel(level));

      if (level >= levels.size() || levels[level]->cells.cells.size() == 0)
        return end_hex();

      return raw_hex_iterator (const_cast<Triangulation<dim,spacedim>*>(this),
                               level,
                               0);
    }

    default:
      Assert (false, ExcNotImplemented());
      return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::hex_iterator
Triangulation<dim,spacedim>::begin_hex (const unsigned int level) const
{
  // level is checked in begin_raw
  raw_hex_iterator ri = begin_raw_hex (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::active_hex_iterator
Triangulation<dim, spacedim>::begin_active_hex (const unsigned int level) const
{
  // level is checked in begin_raw
  hex_iterator i = begin_hex (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::hex_iterator
Triangulation<dim, spacedim>::end_hex () const
{
  return raw_hex_iterator (const_cast<Triangulation<dim,spacedim>*>(this),
                           -1,
                           -1);
}




// -------------------------------- number of cells etc ---------------


namespace internal
{
  namespace Triangulation
  {
    inline
    unsigned int
    n_cells (const internal::Triangulation::NumberCache<1> &c)
    {
      return c.n_lines;
    }


    inline
    unsigned int
    n_active_cells (const internal::Triangulation::NumberCache<1> &c)
    {
      return c.n_active_lines;
    }


    inline
    unsigned int
    n_cells (const internal::Triangulation::NumberCache<2> &c)
    {
      return c.n_quads;
    }


    inline
    unsigned int
    n_active_cells (const internal::Triangulation::NumberCache<2> &c)
    {
      return c.n_active_quads;
    }


    inline
    unsigned int
    n_cells (const internal::Triangulation::NumberCache<3> &c)
    {
      return c.n_hexes;
    }


    inline
    unsigned int
    n_active_cells (const internal::Triangulation::NumberCache<3> &c)
    {
      return c.n_active_hexes;
    }
  }
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_cells () const
{
  return internal::Triangulation::n_cells (number_cache);
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_cells () const
{
  return internal::Triangulation::n_active_cells (number_cache);
}

template <int dim, int spacedim>
types::global_dof_index Triangulation<dim, spacedim>::n_global_active_cells () const
{
  return n_active_cells();
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_faces () const
{
  switch (dim)
    {
    case 1:
      return 0;
    case 2:
      return n_lines();
    case 3:
      return n_quads();
    default:
      Assert (false, ExcNotImplemented());
    }
  return 0;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_raw_faces () const
{
  switch (dim)
    {
    case 2:
      return n_raw_lines();
    case 3:
      return n_raw_quads();
    default:
      Assert (false, ExcNotImplemented());
    }
  return 0;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_faces () const
{
  switch (dim)
    {
    case 1:
      return 0;
    case 2:
      return n_active_lines();
    case 3:
      return n_active_quads();
    default:
      Assert (false, ExcNotImplemented());
    }
  return 0;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_raw_cells (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      return n_raw_lines(level);
    case 2:
      return n_raw_quads(level);
    case 3:
      return n_raw_hexs(level);
    default:
      Assert (false, ExcNotImplemented());
    }
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_cells (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      return n_lines(level);
    case 2:
      return n_quads(level);
    case 3:
      return n_hexs(level);
    default:
      Assert (false, ExcNotImplemented());
    }
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_cells (const unsigned int level) const
{
  switch (dim)
    {
    case 1:
      return n_active_lines(level);
    case 2:
      return n_active_quads(level);
    case 3:
      return n_active_hexs(level);
    default:
      Assert (false, ExcNotImplemented());
    }
  return 0;
}


template <int dim, int spacedim>
bool Triangulation<dim, spacedim>::has_hanging_nodes () const
{
  for (unsigned int lvl = 0; lvl<n_global_levels()-1; lvl++)
    if (n_active_cells(lvl) != 0)
      return true;

  return false;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_lines () const
{
  return number_cache.n_lines;
}


//TODO: Merge the following 6 functions somehow
template <>
unsigned int Triangulation<1,1>::n_raw_lines (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->cells.cells.size();
}


template <>
unsigned int Triangulation<1,1>::n_raw_lines () const
{
  Assert(false, ExcNotImplemented());
  return 0;
}



template <>
unsigned int Triangulation<1,2>::n_raw_lines (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->cells.cells.size();
}


template <>
unsigned int Triangulation<1,2>::n_raw_lines () const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_raw_lines (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->cells.cells.size();
}

template <>
unsigned int Triangulation<1,3>::n_raw_lines () const
{
  Assert(false, ExcNotImplemented());
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_raw_lines (const unsigned int) const
{
  Assert(false, ExcFacesHaveNoLevel());
  return 0;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_raw_lines () const
{
  return faces->lines.cells.size();
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_lines (const unsigned int level) const
{
  Assert (level < number_cache.n_lines_level.size(),
          ExcIndexRange (level, 0, number_cache.n_lines_level.size()));
  Assert (dim == 1, ExcFacesHaveNoLevel());
  return number_cache.n_lines_level[level];
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_lines () const
{
  return number_cache.n_active_lines;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_lines (const unsigned int level) const
{
  Assert (level < number_cache.n_lines_level.size(),
          ExcIndexRange (level, 0, number_cache.n_lines_level.size()));
  Assert (dim == 1, ExcFacesHaveNoLevel());

  return number_cache.n_active_lines_level[level];
}


template <>
unsigned int Triangulation<1,1>::n_quads () const
{
  return 0;
}


template <>
unsigned int Triangulation<1,1>::n_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,1>::n_raw_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,1>::n_raw_hexs (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,1>::n_active_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,1>::n_active_quads () const
{
  return 0;
}




template <>
unsigned int Triangulation<1,2>::n_quads () const
{
  return 0;
}


template <>
unsigned int Triangulation<1,2>::n_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,2>::n_raw_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,2>::n_raw_hexs (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,2>::n_active_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,2>::n_active_quads () const
{
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_quads () const
{
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_raw_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_raw_hexs (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_active_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1,3>::n_active_quads () const
{
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_quads () const
{
  return number_cache.n_quads;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_quads (const unsigned int level) const
{
  Assert (dim == 2, ExcFacesHaveNoLevel());
  Assert (level < number_cache.n_quads_level.size(),
          ExcIndexRange (level, 0, number_cache.n_quads_level.size()));
  return number_cache.n_quads_level[level];
}



template <>
unsigned int Triangulation<2,2>::n_raw_quads (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->cells.cells.size();
}



template <>
unsigned int Triangulation<2,3>::n_raw_quads (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->cells.cells.size();
}


template <>
unsigned int Triangulation<3,3>::n_raw_quads (const unsigned int) const
{
  Assert(false, ExcFacesHaveNoLevel());
  return 0;
}





template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_raw_quads () const
{
  Assert (false, ExcNotImplemented());
  return 0;
}



template <>
unsigned int Triangulation<3,3>::n_raw_quads () const
{
  return faces->quads.cells.size();
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_quads () const
{
  return number_cache.n_active_quads;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_quads (const unsigned int level) const
{
  Assert (level < number_cache.n_quads_level.size(),
          ExcIndexRange (level, 0, number_cache.n_quads_level.size()));
  Assert (dim == 2, ExcFacesHaveNoLevel());

  return number_cache.n_active_quads_level[level];
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_hexs () const
{
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_hexs (const unsigned int) const
{
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_raw_hexs (const unsigned int) const
{
  return 0;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_hexs () const
{
  return 0;
}



template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::n_active_hexs (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<3,3>::n_hexs () const
{
  return number_cache.n_hexes;
}



template <>
unsigned int Triangulation<3,3>::n_hexs (const unsigned int level) const
{
  Assert (level < number_cache.n_hexes_level.size(),
          ExcIndexRange (level, 0, number_cache.n_hexes_level.size()));

  return number_cache.n_hexes_level[level];
}



template <>
unsigned int Triangulation<3,3>::n_raw_hexs (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->cells.cells.size();
}


template <>
unsigned int Triangulation<3,3>::n_active_hexs () const
{
  return number_cache.n_active_hexes;
}



template <>
unsigned int Triangulation<3,3>::n_active_hexs (const unsigned int level) const
{
  Assert (level < number_cache.n_hexes_level.size(),
          ExcIndexRange (level, 0, number_cache.n_hexes_level.size()));

  return number_cache.n_active_hexes_level[level];
}



template <int dim, int spacedim>
unsigned int
Triangulation<dim, spacedim>::n_used_vertices () const
{
  return std::count_if (vertices_used.begin(), vertices_used.end(),
                        std::bind2nd (std::equal_to<bool>(), true));
}



template <int dim, int spacedim>
const std::vector<bool> &
Triangulation<dim, spacedim>::get_used_vertices () const
{
  return vertices_used;
}




template <>
unsigned int Triangulation<1,1>::max_adjacent_cells () const
{
  return 2;
}



template <>
unsigned int Triangulation<1,2>::max_adjacent_cells () const
{
  return 2;
}


template <>
unsigned int Triangulation<1,3>::max_adjacent_cells () const
{
  return 2;
}


template <int dim, int spacedim>
unsigned int Triangulation<dim, spacedim>::max_adjacent_cells () const
{
  cell_iterator cell = begin(0),
                endc = (n_levels() > 1 ? begin(1) : cell_iterator(end()));
  // store the largest index of the
  // vertices used on level 0
  unsigned int max_vertex_index = 0;
  for (; cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      if (cell->vertex_index(vertex) > max_vertex_index)
        max_vertex_index = cell->vertex_index(vertex);

  // store the number of times a cell
  // touches a vertex. An unsigned
  // int should suffice, even for
  // larger dimensions
  std::vector<unsigned short int> usage_count (max_vertex_index+1, 0);
  // touch a vertex's usage count
  // every time we find an adjacent
  // element
  for (cell=begin(); cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      ++usage_count[cell->vertex_index(vertex)];

  return std::max (GeometryInfo<dim>::vertices_per_cell,
                   static_cast<unsigned int>(*std::max_element (usage_count.begin(),
                                             usage_count.end())));
}



template <int dim, int spacedim>
types::subdomain_id
Triangulation<dim,spacedim>::locally_owned_subdomain () const
{
  return numbers::invalid_subdomain_id;
}



template <int dim, int spacedim>
void
Triangulation<dim, spacedim>::execute_coarsening_and_refinement ()
{
  prepare_coarsening_and_refinement ();

  // verify a case with which we have had
  // some difficulty in the past (see the
  // deal.II/coarsening_* tests)
  if (smooth_grid & limit_level_difference_at_vertices)
    Assert (satisfies_level1_at_vertex_rule (*this) == true,
            ExcInternalError());

  // Inform all listeners about beginning of refinement.
  signals.pre_refinement();

  execute_coarsening();

  const DistortedCellList
  cells_with_distorted_children = execute_refinement();

  // verify a case with which we have had
  // some difficulty in the past (see the
  // deal.II/coarsening_* tests)
  if (smooth_grid & limit_level_difference_at_vertices)
    Assert (satisfies_level1_at_vertex_rule (*this) == true,
            ExcInternalError());

  // finally build up neighbor connectivity
  // information
  update_neighbors(*this);

  // Inform all listeners about end of refinement.
  signals.post_refinement();

  AssertThrow (cells_with_distorted_children.distorted_cells.size() == 0,
               cells_with_distorted_children);
}


template<int dim, int spacedim>
void
Triangulation<dim, spacedim>::clear_despite_subscriptions()
{
  // This is the former function
  // clear without the assertion in
  // the beginning.
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.clear ();

  delete faces;
  faces = NULL;

  vertices.clear ();
  vertices_used.clear ();

  manifold.clear();

  number_cache = internal::Triangulation::NumberCache<dim>();
}


template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::DistortedCellList
Triangulation<dim,spacedim>::execute_refinement ()
{
  const DistortedCellList
  cells_with_distorted_children
    =
      internal::Triangulation::Implementation::
      execute_refinement (*this,check_for_distorted_cells);



  // re-compute number of lines
  internal::Triangulation::Implementation
  ::compute_number_cache (*this, levels.size(), number_cache);

#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level)
    levels[level]->cells.monitor_memory (dim);

  // check whether really all
  // refinement flags are reset (also
  // of previously non-active cells
  // which we may not have
  // touched. If the refinement flag
  // of a non-active cell is set,
  // something went wrong since the
  // cell-accessors should have
  // caught this)
  cell_iterator cell = begin(),
                endc = end();
  while (cell != endc)
    Assert (!(cell++)->refine_flag_set(), ExcInternalError ());
#endif

  return cells_with_distorted_children;
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::execute_coarsening ()
{
  // create a vector counting for each line how
  // many cells contain this line. in 3D, this
  // is used later on to decide which lines can
  // be deleted after coarsening a cell. in
  // other dimensions it will be ignored
  std::vector<unsigned int> line_cell_count = count_cells_bounded_by_line (*this);
  std::vector<unsigned int> quad_cell_count = count_cells_bounded_by_quad (*this);

  // loop over all cells. Flag all
  // cells of which all children are
  // flagged for
  // coarsening and delete the childrens'
  // flags. In effect, only those
  // cells are flagged of which originally
  // all children were flagged and for which
  // all children are on the same refinement
  // level. For flagging, the user flags are
  // used, to avoid confusion and because
  // non-active cells can't be flagged for
  // coarsening. Note that because of the
  // effects of @p{fix_coarsen_flags}, of a
  // cell either all or no children must
  // be flagged for coarsening, so it is
  // ok to only check the first child
  clear_user_flags ();

  cell_iterator cell = begin(),
                endc = end();
  for (; cell!=endc; ++cell)
    if (!cell->active())
      if (cell->child(0)->coarsen_flag_set())
        {
          cell->set_user_flag();
          for (unsigned int child=0; child<cell->n_children(); ++child)
            {
              Assert (cell->child(child)->coarsen_flag_set(),
                      ExcInternalError());
              cell->child(child)->clear_coarsen_flag();
            }
        }


  // now do the actual coarsening
  // step. Since the loop goes over
  // used cells we only need not
  // worry about deleting some cells
  // since the ++operator will then
  // just hop over them if we should
  // hit one. Do the loop in the
  // reverse way since we may only
  // delete some cells if their
  // neighbors have already been
  // deleted (if the latter are on a
  // higher level for example)
  //
  // since we delete the *children* of cells, we can ignore cells
  // on the highest level, i.e., level must be less than or equal
  // to n_levels()-2.
  if (levels.size() >= 2)
    for (cell = last(); cell!=endc; --cell)
      if (cell->level()<=static_cast<int>(levels.size()-2) && cell->user_flag_set())
        // use a separate function,
        // since this is dimension
        // specific
        internal::Triangulation::Implementation
        ::delete_children (*this, cell, line_cell_count, quad_cell_count);

  // re-compute number of lines and
  // quads
  internal::Triangulation::Implementation
  ::compute_number_cache (*this, levels.size(), number_cache);

  // in principle no user flags
  // should be
  // set any more at this point
#if DEBUG
  for (cell=begin(); cell!=endc; ++cell)
    Assert (cell->user_flag_set() == false, ExcInternalError());
#endif
}



template <int dim, int spacedim>
void Triangulation<dim, spacedim>::fix_coarsen_flags ()
{
  // copy a piece of code from prepare_coarsening_and_refinement that
  // ensures that the level difference at vertices is limited if so
  // desired. we need this code here since at least in 1d we don't
  // call the dimension-independent version of
  // prepare_coarsening_and_refinement function. in 2d and 3d, having
  // this hunk here makes our lives a bit easier as well as it takes
  // care of these cases earlier than it would otherwise happen.
  //
  // the main difference to the code in p_c_and_r is that here we
  // absolutely have to make sure that we get things right, i.e. that
  // in particular we set flags right if
  // limit_level_difference_at_vertices is set. to do so we iterate
  // until the flags don't change any more
  std::vector<bool> previous_coarsen_flags (n_active_cells());
  save_coarsen_flags (previous_coarsen_flags);

  std::vector<int> vertex_level (vertices.size(), 0);

  bool continue_iterating = true;

  do
    {
      if (smooth_grid & limit_level_difference_at_vertices)
        {
          Assert(!anisotropic_refinement,
                 ExcMessage("In case of anisotropic refinement the "
                            "limit_level_difference_at_vertices flag for "
                            "mesh smoothing must not be set!"));

          // store highest level one of the cells adjacent to a vertex
          // belongs to
          std::fill (vertex_level.begin(), vertex_level.end(), 0);
          active_cell_iterator cell = begin_active(),
                               endc = end();
          for (; cell!=endc; ++cell)
            {
              if (cell->refine_flag_set())
                for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                     ++vertex)
                  vertex_level[cell->vertex_index(vertex)]
                    = std::max (vertex_level[cell->vertex_index(vertex)],
                                cell->level()+1);
              else if (!cell->coarsen_flag_set())
                for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                     ++vertex)
                  vertex_level[cell->vertex_index(vertex)]
                    = std::max (vertex_level[cell->vertex_index(vertex)],
                                cell->level());
              else
                {
                  // if coarsen flag is set then tentatively assume
                  // that the cell will be coarsened. this isn't
                  // always true (the coarsen flag could be removed
                  // again) and so we may make an error here. we try
                  // to correct this by iterating over the entire
                  // process until we are converged
                  Assert (cell->coarsen_flag_set(), ExcInternalError());
                  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                       ++vertex)
                    vertex_level[cell->vertex_index(vertex)]
                      = std::max (vertex_level[cell->vertex_index(vertex)],
                                  cell->level()-1);
                }
            }


          // loop over all cells in reverse order. do so because we
          // can then update the vertex levels on the adjacent
          // vertices and maybe already flag additional cells in this
          // loop
          //
          // note that not only may we have to add additional
          // refinement flags, but we will also have to remove
          // coarsening flags on cells adjacent to vertices that will
          // see refinement
          for (cell=last_active(); cell != endc; --cell)
            if (cell->refine_flag_set() == false)
              {
                for (unsigned int vertex=0;
                     vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
                  if (vertex_level[cell->vertex_index(vertex)] >=
                      cell->level()+1)
                    {
                      // remove coarsen flag...
                      cell->clear_coarsen_flag();

                      // ...and if necessary also refine the current
                      // cell, at the same time updating the level
                      // information about vertices
                      if (vertex_level[cell->vertex_index(vertex)] >
                          cell->level()+1)
                        {
                          cell->set_refine_flag();

                          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell;
                               ++v)
                            vertex_level[cell->vertex_index(v)]
                              = std::max (vertex_level[cell->vertex_index(v)],
                                          cell->level()+1);
                        }

                      // continue and see whether we may, for example,
                      // go into the inner 'if' above based on a
                      // different vertex
                    }
              }
        }

      // loop over all cells. Flag all cells of which all children are
      // flagged for coarsening and delete the childrens' flags. Also
      // delete all flags of cells for which not all children of a
      // cell are flagged. In effect, only those cells are flagged of
      // which originally all children were flagged and for which all
      // children are on the same refinement level. For flagging, the
      // user flags are used, to avoid confusion and because
      // non-active cells can't be flagged for coarsening
      //
      // In effect, all coarsen flags are turned into user flags of
      // the mother cell if coarsening is possible or deleted
      // otherwise.
      clear_user_flags ();
      // Coarsen flags of cells with no mother cell, i.e. on the
      // coarsest level are deleted explicitly.
      active_cell_iterator acell  = begin_active(0),
                           end_ac = end_active(0);
      for (; acell!=end_ac; ++acell)
        acell->clear_coarsen_flag();

      cell_iterator cell = begin(),
                    endc = end();
      for (; cell!=endc; ++cell)
        {
          // nothing to do if we are already on the finest level
          if (cell->active())
            continue;

          const unsigned int n_children=cell->n_children();
          unsigned int flagged_children=0;
          for (unsigned int child=0; child<n_children; ++child)
            if (cell->child(child)->active() &&
                cell->child(child)->coarsen_flag_set())
              {
                ++flagged_children;
                // clear flag since we don't need it anymore
                cell->child(child)->clear_coarsen_flag();
              }

          // flag this cell for coarsening if all children were
          // flagged
          if (flagged_children == n_children)
            cell->set_user_flag();
        }

      // in principle no coarsen flags should be set any more at this
      // point
#if DEBUG
      for (cell=begin(); cell!=endc; ++cell)
        Assert (cell->coarsen_flag_set() == false, ExcInternalError());
#endif

      // now loop over all cells which have the user flag set. their
      // children were flagged for coarsening. set the coarsen flag
      // again if we are sure that none of the neighbors of these
      // children are refined, or will be refined, since then we would
      // get a two-level jump in refinement. on the other hand, if one
      // of the children's neighbors has their user flag set, then we
      // know that its children will go away by coarsening, and we
      // will be ok.
      //
      // note on the other hand that we do allow level-2 jumps in
      // refinement between neighbors in 1d, so this whole procedure
      // is only necessary if we are not in 1d
      //
      // since we remove some coarsening/user flags in the process, we
      // have to work from the finest level to the coarsest one, since
      // we occasionally inspect user flags of cells on finer levels
      // and need to be sure that these flags are final
      for (cell=last(); cell!=endc; --cell)
        if (cell->user_flag_set())
          // if allowed: flag the
          // children for coarsening
          if (internal::Triangulation::Implementation::template coarsening_allowed<dim,spacedim>(cell))
            for (unsigned int c=0; c<cell->n_children(); ++c)
              {
                Assert (cell->child(c)->refine_flag_set()==false,
                        ExcInternalError());

                cell->child(c)->set_coarsen_flag();
              }

      // clear all user flags again, now that we don't need them any
      // more
      clear_user_flags ();


      // now see if anything has changed in the last iteration of this
      // function
      std::vector<bool> current_coarsen_flags (n_active_cells());
      save_coarsen_flags (current_coarsen_flags);

      continue_iterating = (current_coarsen_flags != previous_coarsen_flags);
      previous_coarsen_flags = current_coarsen_flags;
    }
  while (continue_iterating == true);
}


//TODO: merge the following 3 functions since they are the same
template <>
bool Triangulation<1,1>::prepare_coarsening_and_refinement ()
{
  // save the flags to determine whether something was changed in the
  // course of this function
  std::vector<bool> flags_before;
  save_coarsen_flags (flags_before);

  // do nothing in 1d, except setting the coarsening flags correctly
  fix_coarsen_flags ();

  std::vector<bool> flags_after;
  save_coarsen_flags (flags_after);

  return (flags_before != flags_after);
}


template <>
bool Triangulation<1,2>::prepare_coarsening_and_refinement ()
{
  // save the flags to determine whether something was changed in the
  // course of this function
  std::vector<bool> flags_before;
  save_coarsen_flags (flags_before);

  // do nothing in 1d, except setting the coarsening flags correctly
  fix_coarsen_flags ();

  std::vector<bool> flags_after;
  save_coarsen_flags (flags_after);

  return (flags_before != flags_after);
}


template <>
bool Triangulation<1,3>::prepare_coarsening_and_refinement ()
{
  // save the flags to determine whether something was changed in the
  // course of this function
  std::vector<bool> flags_before;
  save_coarsen_flags (flags_before);

  // do nothing in 1d, except setting the coarsening flags correctly
  fix_coarsen_flags ();

  std::vector<bool> flags_after;
  save_coarsen_flags (flags_after);

  return (flags_before != flags_after);
}




namespace
{

  // check if the given @param cell marked for coarsening would
  // produce an unrefined island. To break up long chains of these
  // cells we recursively check our neighbors in case we change this
  // cell. This reduces the number of outer iterations dramatically.
  template <int dim, int spacedim>
  void
  possibly_do_not_produce_unrefined_islands(
    const typename Triangulation<dim,spacedim>::cell_iterator &cell)
  {
    Assert (cell->has_children(), ExcInternalError());

    unsigned int n_neighbors=0;
    // count all neighbors that will be refined along the face of our
    // cell after the next step
    unsigned int count=0;
    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
      {
        const typename Triangulation<dim,spacedim>::cell_iterator neighbor = cell->neighbor(n);
        if (neighbor.state() == IteratorState::valid)
          {
            ++n_neighbors;
            if (face_will_be_refined_by_neighbor(cell,n))
              ++count;
          }
      }
    // clear coarsen flags if either all existing neighbors will be
    // refined or all but one will be and the cell is in the interior
    // of the domain
    if (count==n_neighbors ||
        (count>=n_neighbors-1 &&
         n_neighbors == GeometryInfo<dim>::faces_per_cell) )
      {
        for (unsigned int c=0; c<cell->n_children(); ++c)
          cell->child(c)->clear_coarsen_flag();

        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (!cell->at_boundary(face)
              &&
              ( !cell->neighbor(face)->active() )
              && (cell_will_be_coarsened(cell->neighbor(face))) )
            possibly_do_not_produce_unrefined_islands<dim,spacedim>( cell->neighbor(face) );
      }
  }


  // see if the current cell needs to be refined to avoid unrefined
  // islands.
  //
  // there are sometimes chains of cells that induce refinement of
  // each other. to avoid running the loop in
  // prepare_coarsening_and_refinement over and over again for each
  // one of them, at least for the isotropic refinement case we seek
  // to flag neighboring elements as well as necessary. this takes
  // care of (slightly pathological) cases like
  // deal.II/mesh_smoothing_03
  template <int dim, int spacedim>
  void
  possibly_refine_unrefined_island
  (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
   const bool allow_anisotropic_smoothing)
  {
    Assert (cell->has_children() == false, ExcInternalError());
    Assert (cell->refine_flag_set() == false, ExcInternalError());


    // now we provide two algorithms. the first one is the standard
    // one, coming from the time, where only isotropic refinement was
    // possible. it simply counts the neighbors that are or will be
    // refined and compares to the number of other ones. the second
    // one does this check independently for each direction: if all
    // neighbors in one direction (normally two, at the boundary only
    // one) are refined, the current cell is flagged to be refined in
    // an according direction.

    if (allow_anisotropic_smoothing == false)
      {
        // use first algorithm
        unsigned int refined_neighbors = 0,
                     unrefined_neighbors = 0;
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (!cell->at_boundary(face))
            {
              if (face_will_be_refined_by_neighbor(cell,face))
                ++refined_neighbors;
              else
                ++unrefined_neighbors;
            }

        if (unrefined_neighbors < refined_neighbors)
          {
            cell->clear_coarsen_flag();
            cell->set_refine_flag ();

            // ok, so now we have flagged this cell. if we know that
            // there were any unrefined neighbors at all, see if any
            // of those will have to be refined as well
            if (unrefined_neighbors > 0)
              for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                if (!cell->at_boundary(face)
                    &&
                    (face_will_be_refined_by_neighbor(cell,face) == false)
                    &&
                    (cell->neighbor(face)->has_children() == false)
                    &&
                    (cell->neighbor(face)->refine_flag_set() == false))
                  possibly_refine_unrefined_island<dim,spacedim>
                  (cell->neighbor(face),
                   allow_anisotropic_smoothing);
          }
      }
    else
      {
        // variable to store the cell refine case needed to fulfill
        // all smoothing requirements
        RefinementCase<dim> smoothing_cell_refinement_case
          = RefinementCase<dim>::no_refinement;

        // use second algorithm, do the check individually for each
        // direction
        for (unsigned int face_pair=0;
             face_pair<GeometryInfo<dim>::faces_per_cell/2; ++face_pair)
          {
            // variable to store the cell refine case needed to refine
            // at the current face pair in the same way as the
            // neighbors do...
            RefinementCase<dim> directional_cell_refinement_case
              = RefinementCase<dim>::isotropic_refinement;

            for (unsigned int face_index=0; face_index<2; ++face_index)
              {
                unsigned int face=2*face_pair+face_index;
                // variable to store the refine case (to come) of the
                // face under consideration
                RefinementCase<dim-1> expected_face_ref_case
                  = RefinementCase<dim-1>::no_refinement;

                if (cell->neighbor(face).state() == IteratorState::valid)
                  face_will_be_refined_by_neighbor<dim,spacedim>(cell,face,expected_face_ref_case);
                // now extract which refine case would be necessary to
                // achieve the same face refinement. set the
                // intersection with other requirements for the same
                // direction.

                // note: using the intersection is not an obvious
                // decision, we could also argue that it is more
                // natural to use the union. however, intersection is
                // the less aggressive tactic and favours a smaller
                // number of refined cells over an intensive
                // smoothing. this way we try not to loose too much of
                // the effort we put in anisotropic refinement
                // indicators due to overly aggressive smoothing...
                directional_cell_refinement_case
                  = (directional_cell_refinement_case &
                     GeometryInfo<dim>::min_cell_refinement_case_for_face_refinement(
                       expected_face_ref_case,
                       face,
                       cell->face_orientation(face),
                       cell->face_flip(face),
                       cell->face_rotation(face)));
              }//for both face indices
            // if both requirements sum up to something useful, add
            // this to the refine case for smoothing. note: if
            // directional_cell_refinement_case is isotropic still,
            // then something went wrong...
            Assert(directional_cell_refinement_case <
                   RefinementCase<dim>::isotropic_refinement,
                   ExcInternalError());
            smoothing_cell_refinement_case = smoothing_cell_refinement_case |
                                             directional_cell_refinement_case;
          }//for all face_pairs
        // no we collected contributions from all directions. combine
        // the new flags with the existing refine case, but only if
        // smoothing is required
        if (smoothing_cell_refinement_case)
          {
            cell->clear_coarsen_flag();
            cell->set_refine_flag(cell->refine_flag_set() |
                                  smoothing_cell_refinement_case);
          }
      }
  }
}


template <int dim, int spacedim>
bool Triangulation<dim,spacedim>::prepare_coarsening_and_refinement ()
{
  // save the flags to determine whether something was changed in the
  // course of this function
  std::vector<bool> flags_before[2];
  save_coarsen_flags (flags_before[0]);
  save_refine_flags (flags_before[1]);

  // save the flags at the outset of each loop. we do so in order to
  // find out whether something was changed in the present loop, in
  // which case we would have to re-run the loop. the other
  // possibility to find this out would be to set a flag
  // @p{something_changed} to true each time we change something.
  // however, sometimes one change in one of the parts of the loop is
  // undone by another one, so we might end up in an endless loop. we
  // could be tempted to break this loop at an arbitrary number of
  // runs, but that would not be a clean solution, since we would
  // either have to 1/ break the loop too early, in which case the
  // promise that a second call to this function immediately after the
  // first one does not change anything, would be broken, or 2/ we do
  // as many loops as there are levels. we know that information is
  // transported over one level in each run of the loop, so this is
  // enough. Unfortunately, each loop is rather expensive, so we chose
  // the way presented here
  std::vector<bool> flags_before_loop[2] = {flags_before[0],
                                            flags_before[1]
                                           };

  // now for what is done in each loop: we have to fulfill several
  // tasks at the same time, namely several mesh smoothing algorithms
  // and mesh regularisation, by which we mean that the next mesh
  // fulfills several requirements such as no double refinement at
  // each face or line, etc.
  //
  // since doing these things at once seems almost impossible (in the
  // first year of this library, they were done in two functions, one
  // for refinement and one for coarsening, and most things within
  // these were done at once, so the code was rather impossible to
  // join into this, only, function), we do them one after each
  // other. the order in which we do them is such that the important
  // tasks, namely regularisation, are done last and the least
  // important things are done the first. the following order is
  // chosen:
  //
  // 0/ Only if coarsest_level_1 or patch_level_1 is set: clear all
  //    coarsen flags on level 1 to avoid level 0 cells being created
  //    by coarsening.  As coarsen flags will never be added, this can
  //    be done once and for all before the actual loop starts.
  //
  // 1/ do not coarsen a cell if 'most of the neighbors' will be
  //    refined after the step. This is to prevent occurrence of
  //    unrefined islands.
  //
  // 2/ eliminate refined islands in the interior and at the
  //    boundary. since they don't do much harm besides increasing the
  //    number of degrees of freedom, doing this has a rather low
  //    priority.
  //
  // 3/ limit the level difference of neighboring cells at each
  //    vertex.
  //
  // 4/ eliminate unrefined islands. this has higher priority since
  //    this diminishes the approximation properties not only of the
  //    unrefined island, but also of the surrounding patch.
  //
  // 5/ ensure patch level 1. Then the triangulation consists of
  //    patches, i.e. of cells that are refined once. It follows that
  //    if at least one of the children of a cell is or will be
  //    refined than all children need to be refined. This step only
  //    sets refinement flags and does not set coarsening flags.  If
  //    the patch_level_1 flag is set, then
  //    eliminate_unrefined_islands, eliminate_refined_inner_islands
  //    and eliminate_refined_boundary_islands will be fulfilled
  //    automatically and do not need to be enforced separately.
  //
  // 6/ take care of the requirement that no double refinement is done
  //    at each face
  //
  // 7/ take care that no double refinement is done at each line in 3d
  //    or higher dimensions.
  //
  // 8/ make sure that all children of each cell are either flagged
  //    for coarsening or none of the children is
  //
  // For some of these steps, it is known that they interact. Namely,
  // it is not possible to guarantee that after step 6 another step 5
  // would have no effect; the same holds for the opposite order and
  // also when taking into account step 7. however, it is important to
  // guarantee that step five or six do not undo something that step 5
  // did, and step 7 not something of step 6, otherwise the
  // requirements will not be satisfied even if the loop
  // terminates. this is accomplished by the fact that steps 5 and 6
  // only *add* refinement flags and delete coarsening flags
  // (therefore, step 6 can't undo something that step 4 already did),
  // and step 7 only deletes coarsening flags, never adds some. step 7
  // needs also take care that it won't tag cells for refinement for
  // which some neighbors are more refined or will be refined.

  //////////////////////////////////////
  // STEP 0:
  //    Only if coarsest_level_1 or patch_level_1 is set: clear all
  //    coarsen flags on level 1 to avoid level 0 cells being created
  //    by coarsening.
  if (((smooth_grid & coarsest_level_1) ||
       (smooth_grid & patch_level_1)) && n_levels()>=2)
    {
      active_cell_iterator
      cell=begin_active(1),
      endc=end_active(1);

      for (; cell!=endc; ++cell)
        cell->clear_coarsen_flag();
    }

  bool mesh_changed_in_this_loop = false;
  do
    {
      //////////////////////////////////////
      // STEP 1:
      //    do not coarsen a cell if 'most of the neighbors' will be
      //    refined after the step. This is to prevent the occurrence
      //    of unrefined islands.  If patch_level_1 is set, this will
      //    be automatically fulfilled.
      if (smooth_grid & do_not_produce_unrefined_islands &&
          !(smooth_grid & patch_level_1))
        {
          cell_iterator       cell;
          const cell_iterator endc = end();

          for (cell=begin(); cell!=endc; ++cell)
            {
              // only do something if this
              // cell will be coarsened
              if (!cell->active() && cell_will_be_coarsened(cell))
                possibly_do_not_produce_unrefined_islands<dim,spacedim>(cell);
            }
        }


      //////////////////////////////////////
      // STEP 2:
      //    eliminate refined islands in the interior and at the
      //    boundary. since they don't do much harm besides increasing
      //    the number of degrees of freedom, doing this has a rather
      //    low priority.  If patch_level_1 is set, this will be
      //    automatically fulfilled.
      //
      //    there is one corner case to consider: if this is a
      //    distributed triangulation, there may be refined islands on
      //    the boundary of which we own only part (e.g. a single cell
      //    in the corner of a domain). the rest of the island is
      //    ghost cells and it *looks* like the area around it
      //    (artificial cells) are coarser but this is only because
      //    they may actually be equally fine on other
      //    processors. it's hard to detect this case but we can do
      //    the following: only set coarsen flags to remove this
      //    refined island if all cells we want to set flags on are
      //    locally owned
      if (smooth_grid & (eliminate_refined_inner_islands |
                         eliminate_refined_boundary_islands) &&
          !(smooth_grid & patch_level_1))
        {
          cell_iterator       cell;
          const cell_iterator endc = end();

          for (cell=begin(); cell!=endc; ++cell)
            if (!cell->active() ||
                (cell->active() &&
                 cell->refine_flag_set() &&
                 cell->is_locally_owned()))
              {
                // check whether all children are active, i.e. not
                // refined themselves. This is a precondition that the
                // children may be coarsened away. If the cell is only
                // flagged for refinement, then all future children
                // will be active
                bool all_children_active = true;
                if (!cell->active())
                  for (unsigned int c=0; c<cell->n_children(); ++c)
                    if (!cell->child(c)->active() ||
                        cell->child(c)->is_ghost() ||
                        cell->child(c)->is_artificial())
                      {
                        all_children_active = false;
                        break;
                      }

                if (all_children_active)
                  {
                    // count number of refined and unrefined neighbors
                    // of cell.  neighbors on lower levels are counted
                    // as unrefined since they can only get to the
                    // same level as this cell by the next refinement
                    // cycle
                    unsigned int unrefined_neighbors = 0,
                                 total_neighbors = 0;

                    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
                      {
                        const cell_iterator neighbor = cell->neighbor(n);
                        if (neighbor.state() == IteratorState::valid)
                          {
                            ++total_neighbors;

                            if (!face_will_be_refined_by_neighbor(cell,n))
                              ++unrefined_neighbors;
                          }

                      }

                    // if all neighbors unrefined: mark this cell for
                    // coarsening or don't refine if marked for that
                    //
                    // also do the distinction between the two
                    // versions of the eliminate_refined_*_islands
                    // flag
                    //
                    // the last check is whether there are any
                    // neighbors at all. if not so, then we are (e.g.)
                    // on the coarsest grid with one cell, for which,
                    // of course, we do not remove the refine flag.
                    if ((unrefined_neighbors == total_neighbors)
                        &&
                        (((unrefined_neighbors==GeometryInfo<dim>::faces_per_cell) &&
                          (smooth_grid & eliminate_refined_inner_islands)) ||
                         ((unrefined_neighbors<GeometryInfo<dim>::faces_per_cell) &&
                          (smooth_grid & eliminate_refined_boundary_islands)) )
                        &&
                        (total_neighbors != 0))
                      {
                        if (!cell->active())
                          for (unsigned int c=0; c<cell->n_children(); ++c)
                            {
                              cell->child(c)->clear_refine_flag ();
                              cell->child(c)->set_coarsen_flag ();
                            }
                        else
                          cell->clear_refine_flag();
                      }
                  }
              }
        }

      //////////////////////////////////////
      // STEP 3:
      //    limit the level difference of neighboring cells at each
      //    vertex.
      //
      //    in case of anisotropic refinement this does not make
      //    sense. as soon as one cell is anisotropically refined, an
      //    Assertion is thrown. therefore we can ignore this problem
      //    later on
      if (smooth_grid & limit_level_difference_at_vertices)
        {
          Assert(!anisotropic_refinement,
                 ExcMessage("In case of anisotropic refinement the "
                            "limit_level_difference_at_vertices flag for "
                            "mesh smoothing must not be set!"));

          // store highest level one of the cells adjacent to a vertex
          // belongs to
          std::vector<int> vertex_level (vertices.size(), 0);
          active_cell_iterator cell = begin_active(),
                               endc = end();
          for (; cell!=endc; ++cell)
            {
              if (cell->refine_flag_set())
                for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                     ++vertex)
                  vertex_level[cell->vertex_index(vertex)]
                    = std::max (vertex_level[cell->vertex_index(vertex)],
                                cell->level()+1);
              else if (!cell->coarsen_flag_set())
                for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                     ++vertex)
                  vertex_level[cell->vertex_index(vertex)]
                    = std::max (vertex_level[cell->vertex_index(vertex)],
                                cell->level());
              else
                {
                  // if coarsen flag is set then tentatively assume
                  // that the cell will be coarsened. this isn't
                  // always true (the coarsen flag could be removed
                  // again) and so we may make an error here
                  Assert (cell->coarsen_flag_set(), ExcInternalError());
                  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
                       ++vertex)
                    vertex_level[cell->vertex_index(vertex)]
                      = std::max (vertex_level[cell->vertex_index(vertex)],
                                  cell->level()-1);
                }
            }


          // loop over all cells in reverse order. do so because we
          // can then update the vertex levels on the adjacent
          // vertices and maybe already flag additional cells in this
          // loop
          //
          // note that not only may we have to add additional
          // refinement flags, but we will also have to remove
          // coarsening flags on cells adjacent to vertices that will
          // see refinement
          for (cell=last_active(); cell != endc; --cell)
            if (cell->refine_flag_set() == false)
              {
                for (unsigned int vertex=0;
                     vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
                  if (vertex_level[cell->vertex_index(vertex)] >=
                      cell->level()+1)
                    {
                      // remove coarsen flag...
                      cell->clear_coarsen_flag();

                      // ...and if necessary also refine the current
                      // cell, at the same time updating the level
                      // information about vertices
                      if (vertex_level[cell->vertex_index(vertex)] >
                          cell->level()+1)
                        {
                          cell->set_refine_flag();

                          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell;
                               ++v)
                            vertex_level[cell->vertex_index(v)]
                              = std::max (vertex_level[cell->vertex_index(v)],
                                          cell->level()+1);
                        }

                      // continue and see whether we may, for example,
                      // go into the inner'if'
                      // above based on a
                      // different vertex
                    }
              }
        }

      /////////////////////////////////////
      // STEP 4:
      //    eliminate unrefined islands. this has higher priority
      //    since this diminishes the approximation properties not
      //    only of the unrefined island, but also of the surrounding
      //    patch.
      //
      //    do the loop from finest to coarsest cells since we may
      //    trigger a cascade by marking cells for refinement which
      //    may trigger more cells further down below
      if (smooth_grid & eliminate_unrefined_islands)
        {
          active_cell_iterator cell=last_active(),
                               endc=end();

          for (; cell != endc; --cell)
            // only do something if cell is not already flagged for
            // (isotropic) refinement
            if (cell->refine_flag_set() != RefinementCase<dim>::isotropic_refinement)
              possibly_refine_unrefined_island<dim,spacedim>
              (cell,
               (smooth_grid & allow_anisotropic_smoothing) != 0);
        }

      /////////////////////////////////
      // STEP 5:
      //    ensure patch level 1.
      //
      //    Introduce some terminology:
      //    - a cell that is refined
      //      once is a patch of
      //      level 1 simply called patch.
      //    - a cell that is globally
      //      refined twice is called
      //      a patch of level 2.
      //    - patch level n says that
      //      the triangulation consists
      //      of patches of level n.
      //      This makes sense only
      //      if the grid is already at
      //      least n times globally
      //      refined.
      //
      //    E.g. from patch level 1 follows: if at least one of the
      //    children of a cell is or will be refined than enforce all
      //    children to be refined.

      //    This step 4 only sets refinement flags and does not set
      //    coarsening flags.
      if (smooth_grid & patch_level_1)
        {

          // An important assumption (A) is that before calling this
          // function the grid was already of patch level 1.

          // loop over all cells whose children are all active.  (By
          // assumption (A) either all or none of the children are
          // active).  If the refine flag of at least one of the
          // children is set then set_refine_flag and
          // clear_coarsen_flag of all children.
          for (cell_iterator cell = begin(); cell != end(); ++cell)
            if (!cell->active())
              {
                // ensure the invariant. we can then check whether all
                // of its children are further refined or not by
                // simply looking at the first child
                Assert (cell_is_patch_level_1(cell),
                        ExcInternalError());
                if (cell->child(0)->has_children() == true)
                  continue;

                // cell is found to be a patch.  combine the refine
                // cases of all children
                RefinementCase<dim> combined_ref_case = RefinementCase<dim>::no_refinement;
                for (unsigned int i=0; i<cell->n_children(); ++i)
                  combined_ref_case = combined_ref_case |
                                      cell->child(i)->refine_flag_set();
                if (combined_ref_case != RefinementCase<dim>::no_refinement)
                  for (unsigned int i=0; i<cell->n_children(); ++i)
                    {
                      cell_iterator child = cell->child(i);

                      child->clear_coarsen_flag();
                      child->set_refine_flag(combined_ref_case);
                    }
              }

          // The code above dealt with the case where we may get a
          // non-patch_level_1 mesh from refinement. Now also deal
          // with the case where we could get such a mesh by
          // coarsening.  Coarsen the children (and remove the
          // grandchildren) only if all cell->grandchild(i)
          // ->coarsen_flag_set() are set.
          //
          // for a case where this is a bit tricky, take a look at the
          // mesh_smoothing_0[12] testcases
          for (cell_iterator cell = begin(); cell != end(); ++cell)
            {
              // check if this cell has active grandchildren. note
              // that we know that it is patch_level_1, i.e. if one of
              // its children is active then so are all, and it isn't
              // going to have any grandchildren at all:
              if (cell->active()
                  ||
                  cell->child(0)->active())
                continue;

              // cell is not active, and so are none of its
              // children. check the grandchildren. note that the
              // children are also patch_level_1, and so we only ever
              // need to check their first child
              const unsigned int n_children=cell->n_children();
              bool has_active_grandchildren = false;

              for (unsigned int i=0; i<n_children; ++i)
                if (cell->child(i)->child(0)->active())
                  {
                    has_active_grandchildren = true;
                    break;
                  }

              if (has_active_grandchildren == false)
                continue;


              // ok, there are active grandchildren. see if either all
              // or none of them are flagged for coarsening
              unsigned int n_grandchildren=0;

              // count all coarsen flags of the grandchildren.
              unsigned int n_coarsen_flags=0;

              // cell is not a patch (of level 1) as it has a
              // grandchild.  Is cell a patch of level 2??  Therefore:
              // find out whether all cell->child(i) are patches
              for (unsigned int c=0; c<n_children; ++c)
                {
                  // get at the child. by assumption (A), and the
                  // check by which we got here, the child is not
                  // active
                  cell_iterator child=cell->child(c);

                  const unsigned int nn_children=child->n_children();
                  n_grandchildren += nn_children;

                  // if child is found to be a patch of active cells
                  // itself, then add up how many of its children are
                  // supposed to be coarsened
                  if (child->child(0)->active())
                    for (unsigned int cc=0; cc<nn_children; ++cc)
                      if (child->child(cc)->coarsen_flag_set())
                        ++n_coarsen_flags;
                }

              // if not all grandchildren are supposed to be coarsened
              // (e.g. because some simply don't have the flag set, or
              // because they are not active and therefore cannot
              // carry the flag), then remove the coarsen flag from
              // all of the active grandchildren. note that there may
              // be coarsen flags on the grandgrandchildren -- we
              // don't clear them here, but we'll get to them in later
              // iterations if necessary
              //
              // there is nothing we have to do if no coarsen flags
              // have been set at all
              if ((n_coarsen_flags != n_grandchildren)
                  &&
                  (n_coarsen_flags > 0))
                for (unsigned int c=0; c<n_children; ++c)
                  {
                    const cell_iterator child = cell->child(c);
                    if (child->child(0)->active())
                      for (unsigned int cc=0; cc<child->n_children(); ++cc)
                        child->child(cc)->clear_coarsen_flag();
                  }
            }
        }

      //////////////////////////////////
      //
      //  at the boundary we could end up with cells with negative
      //  volume or at least with a part, that is negative, if the
      //  cell is refined anisotropically. we have to check, whether
      //  that can happen
      internal::Triangulation::Implementation::prevent_distorted_boundary_cells(*this);

      /////////////////////////////////
      // STEP 6:
      //    take care of the requirement that no
      //    double refinement is done at each face
      //
      //    in case of anisotropic refinement it is only likely, but
      //    not sure, that the cells, which are more refined along a
      //    certain face common to two cells are on a higher
      //    level. therefore we cannot be sure, that the requirement
      //    of no double refinement is fulfilled after a single pass
      //    of the following actions. We could just wait for the next
      //    global loop. when this function terminates, the
      //    requirement will be fulfilled. However, it might be faster
      //    to insert an inner loop here.
      bool changed = true;
      while (changed)
        {
          changed=false;
          active_cell_iterator cell=last_active(),
                               endc=end();

          for (; cell != endc; --cell)
            if (cell->refine_flag_set())
              {
                // loop over neighbors of cell
                for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
                  {
                    // only do something if the face is not at the
                    // boundary and if the face will be refined with
                    // the RefineCase currently flagged for
                    if (cell->neighbor(i).state() == IteratorState::valid &&
                        (GeometryInfo<dim>::face_refinement_case(cell->refine_flag_set(),
                                                                 i)
                         != RefinementCase<dim-1>::no_refinement))
                      {
                        // 1) if the neighbor has children: nothing to
                        // worry about.  2) if the neighbor is active
                        // and a coarser one, ensure, that its
                        // refine_flag is set 3) if the neighbor is
                        // active and as refined along the face as our
                        // current cell, make sure, that no
                        // coarsen_flag is set. if we remove the
                        // coarsen flag of our neighbor,
                        // fix_coarsen_flags() makes sure, that the
                        // mother cell will not be coarsened
                        if (cell->neighbor(i)->active())
                          {
                            if (cell->neighbor_is_coarser(i))
                              {
                                if (cell->neighbor(i)->coarsen_flag_set())
                                  cell->neighbor(i)->clear_coarsen_flag();
                                // we'll set the refine flag for this
                                // neighbor below. we note, that we
                                // have changed something by setting
                                // the changed flag to true. We do not
                                // need to do so, if we just removed
                                // the coarsen flag, as the changed
                                // flag only indicates the need to
                                // re-run the inner loop. however, we
                                // only loop over cells flagged for
                                // refinement here, so nothing to
                                // worry about if we remove coarsen
                                // flags

                                if (dim==2)
                                  {
                                    if (smooth_grid & allow_anisotropic_smoothing)
                                      changed=cell->neighbor(i)->flag_for_face_refinement(cell->neighbor_of_coarser_neighbor(i).first,
                                                                                          RefinementCase<dim-1>::cut_x);
                                    else
                                      {
                                        if (!cell->neighbor(i)->refine_flag_set())
                                          changed=true;
                                        cell->neighbor(i)->set_refine_flag();
                                      }
                                  }
                                else //i.e. if (dim==3)
                                  {
// ugly situations might arise here, consider the following situation, which
// shows neighboring cells at the common face, where the upper right element is
// coarser at the given face. Now the upper child element of the lower left
// wants to refine according to cut_z, such that there is a 'horizontal'
// refinement of the face marked with #####
//
//                            /               /
//                           /               /
//                          *---------------*
//                          |               |
//                          |               |
//                          |               |
//                          |               |
//                          |               |
//                          |               | /
//                          |               |/
//                          *---------------*
//
//
//     *---------------*
//    /|              /|
//   / |     #####   / |
//     |               |
//     *---------------*
//    /|              /|
//   / |             / |
//     |               |
//     *---------------*
//    /               /
//   /               /
//
// this introduces too many hanging nodes and the neighboring (coarser) cell
// (upper right) has to be refined. If it is only refined according to cut_z,
// then everything is ok:
//
//                            /               /
//                           /               /
//                          *---------------*
//                          |               |
//                          |               | /
//                          |               |/
//                          *---------------*
//                          |               |
//                          |               | /
//                          |               |/
//                          *---------------*
//
//
//     *---------------*
//    /|              /|
//   / *---------------*
//    /|              /|
//     *---------------*
//    /|              /|
//   / |             / |
//     |               |
//     *---------------*
//    /               /
//   /               /
//
// if however the cell wants to refine itself in an other way, or if we disallow
// anisotropic smoothing, then simply refining the neighbor isotropically is not
// going to work, since this introduces a refinement of face ##### with both
// cut_x and cut_y, which is not possible:
//
//                            /       /       /
//                           /       /       /
//                          *-------*-------*
//                          |       |       |
//                          |       |       | /
//                          |       |       |/
//                          *-------*-------*
//                          |       |       |
//                          |       |       | /
//                          |       |       |/
//                          *-------*-------*
//
//
//     *---------------*
//    /|              /|
//   / *---------------*
//    /|              /|
//     *---------------*
//    /|              /|
//   / |             / |
//     |               |
//     *---------------*
//    /               /
//   /               /
//
// thus, in this case we also need to refine our current cell in the new
// direction:
//
//                            /       /       /
//                           /       /       /
//                          *-------*-------*
//                          |       |       |
//                          |       |       | /
//                          |       |       |/
//                          *-------*-------*
//                          |       |       |
//                          |       |       | /
//                          |       |       |/
//                          *-------*-------*
//
//
//     *-------*-------*
//    /|      /|      /|
//   / *-------*-------*
//    /|      /|      /|
//     *-------*-------*
//    /|      /       /|
//   / |             / |
//     |               |
//     *---------------*
//    /               /
//   /               /

                                    std::pair<unsigned int, unsigned int> nb_indices
                                      =cell->neighbor_of_coarser_neighbor(i);
                                    unsigned int refined_along_x=0,
                                                 refined_along_y=0,
                                                 to_be_refined_along_x=0,
                                                 to_be_refined_along_y=0;

                                    const int this_face_index=cell->face_index(i);

// step 1: detect, along which axis the face is currently refined
                                    if ((this_face_index
                                         == cell->neighbor(i)->face(nb_indices.first)->child_index(0)) ||
                                        (this_face_index
                                         == cell->neighbor(i)->face(nb_indices.first)->child_index(1)))
                                      {
                                        // this might be an
                                        // anisotropic child. get the
                                        // face refine case of the
                                        // neighbors face and count
                                        // refinements in x and y
                                        // direction.
                                        RefinementCase<dim-1> frc=cell->neighbor(i)->face(nb_indices.first)->refinement_case();
                                        if (frc & RefinementCase<dim>::cut_x)
                                          ++refined_along_x;
                                        if (frc & RefinementCase<dim>::cut_y)
                                          ++refined_along_y;
                                      }
                                    else
                                      // this has to be an isotropic
                                      // child
                                      {
                                        ++refined_along_x;
                                        ++refined_along_y;
                                      }
// step 2: detect, along which axis the face has to be refined given the current
// refine flag
                                    RefinementCase<dim-1> flagged_frc=
                                      GeometryInfo<dim>::face_refinement_case(cell->refine_flag_set(),
                                                                              i,
                                                                              cell->face_orientation(i),
                                                                              cell->face_flip(i),
                                                                              cell->face_rotation(i));
                                    if (flagged_frc & RefinementCase<dim>::cut_x)
                                      ++to_be_refined_along_x;
                                    if (flagged_frc & RefinementCase<dim>::cut_y)
                                      ++to_be_refined_along_y;

// step 3: set the refine flag of the (coarser and active) neighbor.
                                    if ((smooth_grid & allow_anisotropic_smoothing) ||
                                        cell->neighbor(i)->refine_flag_set())
                                      {
                                        if (refined_along_x + to_be_refined_along_x > 1)
                                          changed |= cell->neighbor(i)->flag_for_face_refinement(nb_indices.first,
                                                                                                 RefinementCase<dim-1>::cut_axis(0));
                                        if (refined_along_y + to_be_refined_along_y > 1)
                                          changed |= cell->neighbor(i)->flag_for_face_refinement(nb_indices.first,
                                                                                                 RefinementCase<dim-1>::cut_axis(1));
                                      }
                                    else
                                      {
                                        if (cell->neighbor(i)->refine_flag_set()!=RefinementCase<dim>::isotropic_refinement)
                                          changed=true;
                                        cell->neighbor(i)->set_refine_flag();
                                      }

// step 4: if necessary (see above) add to the refine flag of the current cell
                                    cell_iterator nb=cell->neighbor(i);
                                    RefinementCase<dim-1> nb_frc
                                      = GeometryInfo<dim>::face_refinement_case(nb->refine_flag_set(),
                                                                                nb_indices.first,
                                                                                nb->face_orientation(nb_indices.first),
                                                                                nb->face_flip(nb_indices.first),
                                                                                nb->face_rotation(nb_indices.first));
                                    if ((nb_frc & RefinementCase<dim>::cut_x) &&
                                        !(refined_along_x || to_be_refined_along_x))
                                      changed |= cell->flag_for_face_refinement(i,RefinementCase<dim-1>::cut_axis(0));
                                    if ((nb_frc & RefinementCase<dim>::cut_y) &&
                                        !(refined_along_y || to_be_refined_along_y))
                                      changed |= cell->flag_for_face_refinement(i,RefinementCase<dim-1>::cut_axis(1));
                                  }
                              }// if neighbor is coarser
                            else // -> now the neighbor is not coarser
                              {
                                cell->neighbor(i)->clear_coarsen_flag();
                                const unsigned int nb_nb=cell->neighbor_of_neighbor(i);
                                const cell_iterator neighbor=cell->neighbor(i);
                                RefinementCase<dim-1> face_ref_case=
                                  GeometryInfo<dim>::face_refinement_case(neighbor->refine_flag_set(),
                                                                          nb_nb,
                                                                          neighbor->face_orientation(nb_nb),
                                                                          neighbor->face_flip(nb_nb),
                                                                          neighbor->face_rotation(nb_nb));
                                RefinementCase<dim-1> needed_face_ref_case
                                  =GeometryInfo<dim>::face_refinement_case(cell->refine_flag_set(),
                                                                           i,
                                                                           cell->face_orientation(i),
                                                                           cell->face_flip(i),
                                                                           cell->face_rotation(i));
                                // if the neighbor wants to refine the
                                // face with cut_x and we want cut_y
                                // or vice versa, we have to refine
                                // isotropically at the given face
                                if ((face_ref_case==RefinementCase<dim>::cut_x && needed_face_ref_case==RefinementCase<dim>::cut_y) ||
                                    (face_ref_case==RefinementCase<dim>::cut_y && needed_face_ref_case==RefinementCase<dim>::cut_x))
                                  {
                                    changed=cell->flag_for_face_refinement(i, face_ref_case);
                                    neighbor->flag_for_face_refinement(nb_nb, needed_face_ref_case);
                                  }
                              }
                          }
                        else //-> the neighbor is not active
                          {
                            RefinementCase<dim-1> face_ref_case = cell->face(i)->refinement_case(),
                                                  needed_face_ref_case = GeometryInfo<dim>::face_refinement_case(cell->refine_flag_set(),
                                                                         i,
                                                                         cell->face_orientation(i),
                                                                         cell->face_flip(i),
                                                                         cell->face_rotation(i));
                            // if the face is refined with cut_x and
                            // we want cut_y or vice versa, we have to
                            // refine isotropically at the given face
                            if ((face_ref_case==RefinementCase<dim>::cut_x && needed_face_ref_case==RefinementCase<dim>::cut_y) ||
                                (face_ref_case==RefinementCase<dim>::cut_y && needed_face_ref_case==RefinementCase<dim>::cut_x))
                              changed=cell->flag_for_face_refinement(i, face_ref_case);
                          }
                      }
                  }
              }
        }

      //////////////////////////////////////
      // STEP 7:
      //    take care that no double refinement
      //    is done at each line in 3d or higher
      //    dimensions.
      internal::Triangulation::Implementation::prepare_refinement_dim_dependent (*this);

      //////////////////////////////////////
      // STEP 8:
      //    make sure that all children of each
      //    cell are either flagged for coarsening
      //    or none of the children is
      fix_coarsen_flags ();
      // get the refinement and coarsening
      // flags
      std::vector<bool> flags_after_loop[2];
      save_coarsen_flags (flags_after_loop[0]);
      save_refine_flags (flags_after_loop[1]);

      // find out whether something was
      // changed in this loop
      mesh_changed_in_this_loop
        = ((flags_before_loop[0] != flags_after_loop[0]) ||
           (flags_before_loop[1] != flags_after_loop[1]));

      // set the flags for the next loop
      // already
      flags_before_loop[0].swap(flags_after_loop[0]);
      flags_before_loop[1].swap(flags_after_loop[1]);
    }
  while (mesh_changed_in_this_loop);


  // find out whether something was really changed in this
  // function. Note that @p{flags_before_loop} represents the state
  // after the last loop, i.e.  the present state
  return ((flags_before[0] != flags_before_loop[0]) ||
          (flags_before[1] != flags_before_loop[1]));
}




template <int dim, int spacedim>
void Triangulation<dim, spacedim>::write_bool_vector (const unsigned int  magic_number1,
                                                      const std::vector<bool> &v,
                                                      const unsigned int  magic_number2,
                                                      std::ostream            &out)
{
  const unsigned int N = v.size();
  unsigned char *flags = new unsigned char[N/8+1];
  for (unsigned int i=0; i<N/8+1; ++i) flags[i]=0;

  for (unsigned int position=0; position<N; ++position)
    flags[position/8] |= (v[position] ? (1<<(position%8)) : 0);

  AssertThrow (out, ExcIO());

  // format:
  // 0. magic number
  // 1. number of flags
  // 2. the flags
  // 3. magic number
  out << magic_number1 << ' ' << N << std::endl;
  for (unsigned int i=0; i<N/8+1; ++i)
    out << static_cast<unsigned int>(flags[i]) << ' ';

  out << std::endl << magic_number2 << std::endl;

  delete[] flags;

  AssertThrow (out, ExcIO());
}


template <int dim, int spacedim>
void Triangulation<dim, spacedim>::read_bool_vector (const unsigned int  magic_number1,
                                                     std::vector<bool>       &v,
                                                     const unsigned int  magic_number2,
                                                     std::istream            &in)
{
  AssertThrow (in, ExcIO());

  unsigned int magic_number;
  in >> magic_number;
  AssertThrow (magic_number==magic_number1, ExcGridReadError());

  unsigned int N;
  in >> N;
  v.resize (N);

  unsigned char *flags = new unsigned char[N/8+1];
  unsigned short int tmp;
  for (unsigned int i=0; i<N/8+1; ++i)
    {
      in >> tmp;
      flags[i] = tmp;
    }

  for (unsigned int position=0; position!=N; ++position)
    v[position] = (flags[position/8] & (1<<(position%8)));

  in >> magic_number;
  AssertThrow (magic_number==magic_number2, ExcGridReadError());

  delete[] flags;

  AssertThrow (in, ExcIO());
}



template <int dim, int spacedim>
std::size_t
Triangulation<dim, spacedim>::memory_consumption () const
{
  std::size_t mem = 0;
  mem += MemoryConsumption::memory_consumption(levels);
  for (unsigned int i=0; i<levels.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*levels[i]);
  mem += MemoryConsumption::memory_consumption (vertices);
  mem += MemoryConsumption::memory_consumption (vertices_used);
  mem += sizeof(manifold);
  mem += sizeof(smooth_grid);
  mem += MemoryConsumption::memory_consumption (number_cache);
  mem += sizeof (faces);
  mem += MemoryConsumption::memory_consumption (*faces);

  return mem;
}




template<int dim, int spacedim>
Triangulation<dim, spacedim>::RefinementListener::~RefinementListener ()
{}



template<int dim, int spacedim>
Triangulation<dim, spacedim>::DistortedCellList::~DistortedCellList () throw ()
{
  // don't do anything here. the compiler will automatically convert
  // any exceptions created by the destructors of the member variables
  // into abort() in order to satisfy the throw() specification
}




template<int dim, int spacedim>
void Triangulation<dim, spacedim>::
RefinementListener::pre_refinement_notification (const Triangulation<dim, spacedim> &)
{}



template<int dim, int spacedim>
void Triangulation<dim, spacedim>::
RefinementListener::post_refinement_notification (const Triangulation<dim, spacedim> &)
{}



template<int dim, int spacedim>
void Triangulation<dim, spacedim>::
RefinementListener::copy_notification (const Triangulation<dim, spacedim> &,
                                       const Triangulation<dim, spacedim> &)
{}



template<int dim, int spacedim>
void Triangulation<dim, spacedim>::
RefinementListener::create_notification (const Triangulation<dim, spacedim> &)
{}



template<int dim, int spacedim>
void
Triangulation<dim, spacedim>::add_refinement_listener (RefinementListener &listener) const
{
  // in this compatibility mode with the old-style refinement
  // listeners, an external class presents itself as one that may or
  // may not have overloaded all of the functions that the
  // RefinementListener class has. consequently, we need to connect
  // each of its functions to the relevant signals. for those
  // functions that haven't been overloaded, that means that
  // triggering the signal yields a call to the function in the
  // RefinementListener base class which simply does nothing
  std::vector<boost::signals2::connection> connections;

  connections.push_back
  (signals.create.connect (std_cxx11::bind (&RefinementListener::create_notification,
                                            std_cxx11::ref(listener),
                                            std_cxx11::cref(*this))));
  connections.push_back
  (signals.copy.connect (std_cxx11::bind (&RefinementListener::copy_notification,
                                          std_cxx11::ref(listener),
                                          std_cxx11::cref(*this),
                                          std_cxx11::_1)));
  connections.push_back
  (signals.pre_refinement.connect (std_cxx11::bind (&RefinementListener::pre_refinement_notification,
                                                    std_cxx11::ref(listener),
                                                    std_cxx11::cref(*this))));
  connections.push_back
  (signals.post_refinement.connect (std_cxx11::bind (&RefinementListener::post_refinement_notification,
                                                     std_cxx11::ref(listener),
                                                     std_cxx11::cref(*this))));

  // now push the set of connections into the map
  refinement_listener_map.insert (std::make_pair(&listener, connections));
}



template<int dim, int spacedim>
void
Triangulation<dim, spacedim>::remove_refinement_listener (RefinementListener &listener) const
{
  Assert (refinement_listener_map.find (&listener) != refinement_listener_map.end(),
          ExcMessage("You try to remove a refinement listener that does "
                     "not appear to have been added previously."));

  // get the element of the map, and terminate these connections. then
  // erase the element from the list
  std::vector<boost::signals2::connection> connections
    = refinement_listener_map.find(&listener)->second;
  for (unsigned int i=0; i<connections.size(); ++i)
    connections[i].disconnect ();

  refinement_listener_map.erase (refinement_listener_map.find (&listener));
}

template <>
const Manifold<2,1> &Triangulation<2, 1>::get_manifold(const types::manifold_id) const
{
  Assert(false, ExcImpossibleInDim(1));
  // We cannot simply create a temporary Manifold<2,1> because it is not
  // instantiated and would lead to unresolved symbols. Given the fact that
  // this function should be unreachable anyaway, just dereference a
  // nullptr:
  return *static_cast<FlatManifold<2,1>*>(0);
}

template <>
const Manifold<3,1> &Triangulation<3, 1>::get_manifold(const types::manifold_id) const
{
  Assert(false, ExcImpossibleInDim(1));
  // We cannot simply create a temporary Manifold<2,1> because it is not
  // instantiated and would lead to unresolved symbols. Given the fact that
  // this function should be unreachable anyaway, just dereference a
  // nullptr:
  return *static_cast<FlatManifold<3,1>*>(0);
}

template <>
const Manifold<3,2> &Triangulation<3, 2>::get_manifold(const types::manifold_id) const
{
  Assert(false, ExcImpossibleInDim(2));
  // We cannot simply create a temporary Manifold<2,1> because it is not
  // instantiated and would lead to unresolved symbols. Given the fact that
  // this function should be unreachable anyaway, just dereference a
  // nullptr:
  return *static_cast<FlatManifold<3,2>*>(0);
}

// explicit instantiations
#include "tria.inst"

DEAL_II_NAMESPACE_CLOSE

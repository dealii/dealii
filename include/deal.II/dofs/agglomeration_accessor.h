// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef dealii_agglomeration_accessor_h
#define dealii_agglomeration_accessor_h

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/iterator_range.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/filtered_iterator.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declarations
#ifndef DOXYGEN
template <int, int>
class AgglomerationHandler;
template <int, int>
class AgglomerationIterator;
#endif


/**
 * Accessor class used by AgglomerationIterator to access agglomeration data.
 */
template <int dim, int spacedim = dim>
class AgglomerationAccessor
{
public:
  /**
   * Type for storing the polygons in an agglomerate.
   */
  using AgglomerationContainer =
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>;


  /**
   * Get the DoFs indices associated to the agglomerate.
   */
  void
  get_dof_indices(std::vector<types::global_dof_index> &) const;

  /**
   * Return, for a cell, the number of faces. In case the cell is a standard
   * cell, then the number of faces is the classical one. If it's a master cell,
   * then it returns the number of faces of the agglomeration identified by the
   * master cell itself.
   */
  unsigned int
  n_faces() const;

  /**
   * Return the number of deal.II faces that is building a polygon face.
   */
  unsigned int
  n_agglomerated_faces() const;

  /**
   * Return the agglomerate which shares face f.
   */
  const AgglomerationIterator<dim, spacedim>
  neighbor(const unsigned int f) const;

  /**
   * Return the present index (seen from the neighboring agglomerate) of the
   * present face f.
   */
  unsigned int
  neighbor_of_agglomerated_neighbor(const unsigned int f) const;

  /**
   *
   * This function generalizes the behaviour of cell->face(f)->at_boundary()
   * in the case where f is an index out of the range [0,..., n_faces).
   * In practice, if you call this function with a standard deal.II cell, you
   * have precisely the same result as calling cell->face(f)->at_boundary().
   * Otherwise, if the cell is a master one, you have a boolean returning true
   * is that face for the agglomeration is on the boundary or not.
   */
  bool
  at_boundary(const unsigned int f) const;

  /**
   * Return a vector of face iterators describing the boundary of agglomerate.
   */
  const std::vector<typename Triangulation<dim>::active_face_iterator> &
  polytope_boundary() const;

  /**
   *
   * Return the volume of a polytope.
   */
  double
  volume() const;

  /**
   * Return the diameter of the present polytopal element.
   */
  double
  diameter() const;

  /**
   * Returns the deal.II cells that build the agglomerate.
   */
  AgglomerationContainer
  get_agglomerate() const;

  /**
   * Return the BoundingBox which bounds the present polytope. In case the
   * present polytope is not locally owned, it returns the BoundingBox of that
   * ghosted polytope.
   */
  const BoundingBox<dim> &
  get_bounding_box() const;

  /**
   * Return the index of the present polytope.
   */
  types::global_cell_index
  index() const;

  /**
   * Returns an active cell iterator for the dof_handler, matching the polytope
   * referenced by the input iterator. The type of the returned object is a
   * DoFHandler::active_cell_iterator which can be used to initialize
   * FiniteElement data.
   */
  typename DoFHandler<dim, spacedim>::active_cell_iterator
  as_dof_handler_iterator(const DoFHandler<dim, spacedim> &dof_handler) const;

  /**
   * Returns the number of classical deal.II cells that are building the present
   * polygon.
   */
  unsigned int
  n_background_cells() const;

  /* Returns true if this polygon is owned by the current processor. On a serial
   * Triangulation this returs always true, but may yield false for a
   * parallel::distributed::Triangulation.
   */
  bool
  is_locally_owned() const;

  /**
   * The polytopal analogue of CellAccessor::id(). It provides a way to uniquely
   * identify cells in a parallel Triangulation such as a
   * parallel::distributed::Triangulation.
   */
  CellId
  id() const;

  /**
   * The polytopal analogue of CellAccessor::subdomain_id(). In case of a serial
   * Triangulation, it returns the numbers::invalid_subdomain_id.
   */
  types::subdomain_id
  subdomain_id() const;

  /**
   * Returns a vector of indices identifying the children polytopes.
   */
  inline const std::vector<types::global_cell_index> &
  children() const;

  /**
   * Returns the FiniteElement object used by the current polytope.
   * This function should only be called after the corresponding agglomeration
   * handler has invoked distribute_agglomerated_dofs().
   */
  const FiniteElement<dim, spacedim> &
  get_fe() const;

  /**
   * Sets the active finite element index.
   * This function should be called when using hp::FECollection to specify
   * which finite element in the collection is assigned to the current polytope.
   */
  void
  set_active_fe_index(const types::fe_index index) const;

  /**
   * Returns the index of the active finite element.
   * When using hp::FECollection, this function retrieves the index
   * of the finite element assigned to the current polytope.
   */
  types::fe_index
  active_fe_index() const;

private:
  /**
   * Private default constructor. This is not supposed to be used and hence will
   * throw.
   */
  AgglomerationAccessor();

  /**
   * Private constructor for an agglomerate. This is meant to be invoked by
   * the AgglomerationIterator class. It takes as input the master cell of the
   * agglomerate and a pointer to the handler.
   */
  AgglomerationAccessor(
    const typename Triangulation<dim, spacedim>::active_cell_iterator
                                              &master_cell,
    const AgglomerationHandler<dim, spacedim> *ah);

  /**
   * Same as above, but needed when the argument @p cells is a ghost cell.
   */
  AgglomerationAccessor(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
    const CellId                                                      &cell_id,
    const AgglomerationHandler<dim, spacedim>                         *ah);

  /**
   * Default destructor.
   */
  ~AgglomerationAccessor() = default;


  /**
   * The unique deal.II cell associated to the present polytope.
   */
  typename Triangulation<dim, spacedim>::active_cell_iterator master_cell;

  /**
   * The index of the present polytope.
   */
  types::global_cell_index present_index;

  /**
   * The index of the present polytope.
   */
  CellId present_id;

  /**
   * The rank owning of the present polytope.
   */
  types::subdomain_id present_subdomain_id;

  /**
   * A pointer to the Handler.
   */
  AgglomerationHandler<dim, spacedim> *handler;

  /**
   * Comparison operator for Accessor. Two accessors are equal if they refer to
   * the same polytopal element.
   */
  bool
  operator==(const AgglomerationAccessor<dim, spacedim> &other) const;

  /**
   * Compare for inequality.
   */
  bool
  operator!=(const AgglomerationAccessor<dim, spacedim> &other) const;

  /**
   * Move to the next cell in the polytopal mesh.
   */
  void
  next();

  /**
   * Move to the previous cell in the polytopal mesh.
   */
  void
  prev();

  /**
   * Returns the slaves of the present agglomeration.
   */
  const AgglomerationContainer &
  get_slaves() const;

  unsigned int
  n_agglomerated_faces_per_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
    const;

  template <int, int>
  friend class AgglomerationIterator;
};



template <int dim, int spacedim>
unsigned int
AgglomerationAccessor<dim, spacedim>::n_agglomerated_faces_per_cell(
  const typename Triangulation<dim, spacedim>::active_cell_iterator &cell) const
{
  unsigned int n_neighbors = 0;
  for (const auto &f : cell->face_indices())
    {
      const auto &neighboring_cell = cell->neighbor(f);
      if ((cell->face(f)->at_boundary()) ||
          (neighboring_cell->is_active() &&
           !handler->are_cells_agglomerated(cell, neighboring_cell)))
        {
          ++n_neighbors;
        }
    }
  return n_neighbors;
}



template <int dim, int spacedim>
unsigned int
AgglomerationAccessor<dim, spacedim>::n_faces() const
{
  Assert(!handler->is_slave_cell(master_cell),
         ExcMessage("You cannot pass a slave cell."));
  return handler->number_of_agglomerated_faces[present_index];
}



template <int dim, int spacedim>
const AgglomerationIterator<dim, spacedim>
AgglomerationAccessor<dim, spacedim>::neighbor(const unsigned int f) const
{
  if (!at_boundary(f))
    {
      if (master_cell->is_ghost())
        {
          // The following path is needed when the present function is called
          // from neighbor_of_neighbor()

          const unsigned int sender_rank = master_cell->subdomain_id();

          const CellId &master_id_ghosted_neighbor =
            handler->recv_ghosted_master_id.at(sender_rank)
              .at(present_id)
              .at(f);

          // Use the id of the master cell to uniquely identify the neighboring
          // agglomerate

          return {master_cell, master_id_ghosted_neighbor, handler};
        }

      const types::global_cell_index polytope_index =
        handler->master2polygon.at(master_cell->active_cell_index());


      const auto &neigh =
        handler->polytope_cache.cell_face_at_boundary.at({polytope_index, f})
          .second;


      if (neigh->is_locally_owned())
        {
          typename DoFHandler<dim, spacedim>::active_cell_iterator cell_dh(
            *neigh, &(handler->agglo_dh));
          return {cell_dh, handler};
        }
      else
        {
          // Get master_id from the neighboring ghost polytope. This uniquely
          // identifies the neighboring polytope among all processors.
          const CellId &master_id_neighbor =
            handler->polytope_cache.ghosted_master_id.at({present_id, f});

          // Use the id of the master cell to uniquely identify the neighboring
          // agglomerate
          return {neigh, master_id_neighbor, handler};
        }
    }
  else
    {
      return {};
    }
}



template <int dim, int spacedim>
unsigned int
AgglomerationAccessor<dim, spacedim>::neighbor_of_agglomerated_neighbor(
  const unsigned int f) const
{
  // First, make sure it's not a boundary face.
  if (!at_boundary(f))
    {
      const auto &neigh_polytope =
        neighbor(f); // returns the neighboring master and id

      AssertThrow(neigh_polytope.state() == IteratorState::valid,
                  ExcInternalError());

      unsigned int n_faces_agglomerated_neighbor;

      // if it is locally owned, retrieve the number of faces
      if (neigh_polytope->is_locally_owned())
        {
          n_faces_agglomerated_neighbor = neigh_polytope->n_faces();
        }
      else
        {
          // The neighboring polytope is not locally owned. We need to get the
          // number of its faces from the neighboring rank.

          // First, retrieve the CellId of the neighboring polytope.
          const CellId &master_id_neighbor = neigh_polytope->id();

          // Then, get the neighboring rank
          const unsigned int sender_rank = neigh_polytope->subdomain_id();

          // From the neighboring rank, use the CellId of the neighboring
          // polytope to get the number of its faces.
          n_faces_agglomerated_neighbor =
            handler->recv_n_faces.at(sender_rank).at(master_id_neighbor);
        }


      // Loop over all faces of neighboring agglomerate
      for (unsigned int f_out = 0; f_out < n_faces_agglomerated_neighbor;
           ++f_out)
        {
          // Check if same CellId
          if (neigh_polytope->neighbor(f_out).state() == IteratorState::valid)
            if (neigh_polytope->neighbor(f_out)->id() == present_id)
              return f_out;
        }
      return numbers::invalid_unsigned_int;
    }
  else
    {
      // Face is at boundary
      return numbers::invalid_unsigned_int;
    }
}

// ------------------------------ inline functions -------------------------

template <int dim, int spacedim>
inline AgglomerationAccessor<dim, spacedim>::AgglomerationAccessor()
{}



template <int dim, int spacedim>
inline AgglomerationAccessor<dim, spacedim>::AgglomerationAccessor(
  const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
  const AgglomerationHandler<dim, spacedim>                         *ah)
{
  handler = const_cast<AgglomerationHandler<dim, spacedim> *>(ah);
  if (&(*handler->master_cells_container.end()) == std::addressof(cell))
    {
      present_index        = handler->master_cells_container.size();
      master_cell          = *handler->master_cells_container.end();
      present_id           = CellId(); // invalid id (TODO)
      present_subdomain_id = numbers::invalid_subdomain_id;
    }
  else
    {
      present_index = handler->master2polygon.at(cell->active_cell_index());
      master_cell   = cell;
      present_id    = master_cell->id();
      present_subdomain_id = master_cell->subdomain_id();
    }
}



template <int dim, int spacedim>
inline AgglomerationAccessor<dim, spacedim>::AgglomerationAccessor(
  const typename Triangulation<dim, spacedim>::active_cell_iterator &neigh_cell,
  const CellId                              &master_cell_id,
  const AgglomerationHandler<dim, spacedim> *ah)
{
  Assert(neigh_cell->is_ghost(), ExcInternalError());
  // neigh_cell is ghosted

  handler       = const_cast<AgglomerationHandler<dim, spacedim> *>(ah);
  master_cell   = neigh_cell;
  present_index = numbers::invalid_unsigned_int;
  // neigh_cell is ghosted, use the CellId of that agglomerate
  present_id           = master_cell_id;
  present_subdomain_id = master_cell->subdomain_id();
}



template <int dim, int spacedim>
inline void
AgglomerationAccessor<dim, spacedim>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices) const
{
  Assert(dof_indices.size() > 0,
         ExcMessage(
           "The vector of DoFs indices must be already properly resized."));
  if (is_locally_owned())
    {
      // Forward the call to the master cell
      typename DoFHandler<dim, spacedim>::cell_iterator master_cell_dh(
        *master_cell, &(handler->agglo_dh));
      master_cell_dh->get_dof_indices(dof_indices);
    }
  else
    {
      const std::vector<types::global_dof_index> &recv_dof_indices =
        handler->recv_ghost_dofs.at(present_subdomain_id).at(present_id);

      std::copy(recv_dof_indices.cbegin(),
                recv_dof_indices.cend(),
                dof_indices.begin());
    }
}



template <int dim, int spacedim>
inline typename AgglomerationAccessor<dim, spacedim>::AgglomerationContainer
AgglomerationAccessor<dim, spacedim>::get_agglomerate() const
{
  auto agglomeration = get_slaves();
  agglomeration.push_back(master_cell);
  return agglomeration;
}



template <int dim, int spacedim>
inline const std::vector<typename Triangulation<dim>::active_face_iterator> &
AgglomerationAccessor<dim, spacedim>::polytope_boundary() const
{
  return handler->polygon_boundary[master_cell];
}



template <int dim, int spacedim>
inline double
AgglomerationAccessor<dim, spacedim>::diameter() const
{
  Assert(!handler->is_slave_cell(master_cell),
         ExcMessage("The present function cannot be called for slave cells."));

  if (handler->is_master_cell(master_cell))
    {
      // Get the bounding box associated with the master cell
      const auto &bdary_pts =
        handler->bboxes[present_index].get_boundary_points();
      return (bdary_pts.second - bdary_pts.first).norm();
    }
  else
    {
      // Standard deal.II way to get the measure of a cell.
      return master_cell->diameter();
    }
}



template <int dim, int spacedim>
inline const BoundingBox<dim> &
AgglomerationAccessor<dim, spacedim>::get_bounding_box() const
{
  if (is_locally_owned())
    return handler->bboxes[present_index];
  else
    return handler->recv_ghosted_bbox.at(present_subdomain_id).at(present_id);
}



template <int dim, int spacedim>
inline double
AgglomerationAccessor<dim, spacedim>::volume() const
{
  Assert(!handler->is_slave_cell(master_cell),
         ExcMessage("The present function cannot be called for slave cells."));

  if (handler->is_master_cell(master_cell))
    {
      return handler->bboxes[present_index].volume();
    }
  else
    {
      return master_cell->measure();
    }
}



template <int dim, int spacedim>
inline void
AgglomerationAccessor<dim, spacedim>::next()
{
  // Increment the present index and update the polytope
  ++present_index;

  // Make sure not to query the CellId if it's past the last
  if (present_index < handler->master_cells_container.size())
    {
      master_cell          = handler->master_cells_container[present_index];
      present_id           = master_cell->id();
      present_subdomain_id = master_cell->subdomain_id();
    }
}



template <int dim, int spacedim>
inline void
AgglomerationAccessor<dim, spacedim>::prev()
{
  // Decrement the present index and update the polytope
  --present_index;
  master_cell = handler->master_cells_container[present_index];
  present_id  = master_cell->id();
}


template <int dim, int spacedim>
inline bool
AgglomerationAccessor<dim, spacedim>::operator==(
  const AgglomerationAccessor<dim, spacedim> &other) const
{
  return present_index == other.present_index;
}

template <int dim, int spacedim>
inline bool
AgglomerationAccessor<dim, spacedim>::operator!=(
  const AgglomerationAccessor<dim, spacedim> &other) const
{
  return !(*this == other);
}



template <int dim, int spacedim>
inline types::global_cell_index
AgglomerationAccessor<dim, spacedim>::index() const
{
  return present_index;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_cell_iterator
AgglomerationAccessor<dim, spacedim>::as_dof_handler_iterator(
  const DoFHandler<dim, spacedim> &dof_handler) const
{
  // Forward the call to the master cell using the right DoFHandler.
  return master_cell->as_dof_handler_iterator(dof_handler);
}



template <int dim, int spacedim>
inline const typename AgglomerationAccessor<dim,
                                            spacedim>::AgglomerationContainer &
AgglomerationAccessor<dim, spacedim>::get_slaves() const
{
  return handler->master2slaves.at(master_cell->active_cell_index());
}



template <int dim, int spacedim>
inline unsigned int
AgglomerationAccessor<dim, spacedim>::n_background_cells() const
{
  AssertThrow(get_agglomerate().size() > 0, ExcMessage("Empty agglomeration."));
  return get_agglomerate().size();
}



template <int dim, int spacedim>
unsigned int
AgglomerationAccessor<dim, spacedim>::n_agglomerated_faces() const
{
  const auto  &agglomeration = get_agglomerate();
  unsigned int n_neighbors   = 0;
  for (const auto &cell : agglomeration)
    n_neighbors += n_agglomerated_faces_per_cell(cell);
  return n_neighbors;
}



template <int dim, int spacedim>
inline bool
AgglomerationAccessor<dim, spacedim>::at_boundary(const unsigned int f) const
{
  if (master_cell->is_ghost())
    {
      const unsigned int sender_rank = master_cell->subdomain_id();
      return handler->recv_bdary_info.at(sender_rank).at(present_id).at(f);
    }
  else
    {
      Assert(!handler->is_slave_cell(master_cell),
             ExcMessage(
               "This function should not be called for a slave cell."));


      typename DoFHandler<dim, spacedim>::active_cell_iterator cell_dh(
        *master_cell, &(handler->agglo_dh));
      return handler->at_boundary(cell_dh, f);
    }
}



template <int dim, int spacedim>
inline bool
AgglomerationAccessor<dim, spacedim>::is_locally_owned() const
{
  return master_cell->is_locally_owned();
}



template <int dim, int spacedim>
inline CellId
AgglomerationAccessor<dim, spacedim>::id() const
{
  return present_id;
}



template <int dim, int spacedim>
inline types::subdomain_id
AgglomerationAccessor<dim, spacedim>::subdomain_id() const
{
  return present_subdomain_id;
}

template <int dim, int spacedim>
inline const std::vector<types::global_cell_index> &
AgglomerationAccessor<dim, spacedim>::children() const
{
  Assert(!handler->parent_child_info.empty(), ExcInternalError());
  return handler->parent_child_info.at(
    {present_index, handler->present_extraction_level});
}

template <int dim, int spacedim>
inline const FiniteElement<dim, spacedim> &
AgglomerationAccessor<dim, spacedim>::get_fe() const
{
  typename DoFHandler<dim, spacedim>::active_cell_iterator
    master_cell_as_dof_handler_iterator =
      master_cell->as_dof_handler_iterator(handler->agglo_dh);
  return master_cell_as_dof_handler_iterator->get_fe();
}

template <int dim, int spacedim>
inline void
AgglomerationAccessor<dim, spacedim>::set_active_fe_index(
  const types::fe_index index) const
{
  Assert(!handler->is_slave_cell(master_cell),
         ExcMessage("The present function cannot be called for slave cells."));
  typename DoFHandler<dim, spacedim>::active_cell_iterator
    master_cell_as_dof_handler_iterator =
      master_cell->as_dof_handler_iterator(handler->agglo_dh);
  master_cell_as_dof_handler_iterator->set_active_fe_index(index);
}

template <int dim, int spacedim>
inline types::fe_index
AgglomerationAccessor<dim, spacedim>::active_fe_index() const
{
  typename DoFHandler<dim, spacedim>::active_cell_iterator
    master_cell_as_dof_handler_iterator =
      master_cell->as_dof_handler_iterator(handler->agglo_dh);
  return master_cell_as_dof_handler_iterator->active_fe_index();
}

DEAL_II_NAMESPACE_CLOSE

#endif

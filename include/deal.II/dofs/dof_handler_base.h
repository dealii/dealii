// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_dof_handler_base_h
#define dealii_dof_handler_base_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/distributed/cell_data_transfer.templates.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/deprecated_function_map.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/number_cache.h>

#include <deal.II/hp/dof_faces.h>
#include <deal.II/hp/dof_level.h>
#include <deal.II/hp/fe_collection.h>

#include <boost/serialization/split_member.hpp>

#include <map>
#include <memory>
#include <set>
#include <unordered_set>
#include <vector>

DEAL_II_NAMESPACE_OPEN


#ifndef DOXYGEN

/* ----------------------- Inline functions ----------------------------------
 */


template <int dim, int spacedim>
inline bool
DoFHandler<dim, spacedim>::has_level_dofs() const
{
  return mg_number_cache.size() > 0;
}



template <int dim, int spacedim>
inline bool
DoFHandler<dim, spacedim>::has_active_dofs() const
{
  return number_cache.n_global_dofs > 0;
}



template <int dim, int spacedim>
inline types::global_dof_index
DoFHandler<dim, spacedim>::n_dofs() const
{
  return number_cache.n_global_dofs;
}



template <int dim, int spacedim>
inline types::global_dof_index
DoFHandler<dim, spacedim>::n_dofs(const unsigned int level) const
{
  Assert(has_level_dofs(),
         ExcMessage(
           "n_dofs(level) can only be called after distribute_mg_dofs()"));
  Assert(level < mg_number_cache.size(), ExcInvalidLevel(level));
  return mg_number_cache[level].n_global_dofs;
}



template <int dim, int spacedim>
types::global_dof_index
DoFHandler<dim, spacedim>::n_locally_owned_dofs() const
{
  return number_cache.n_locally_owned_dofs;
}



template <int dim, int spacedim>
const IndexSet &
DoFHandler<dim, spacedim>::locally_owned_dofs() const
{
  return number_cache.locally_owned_dofs;
}



template <int dim, int spacedim>
const IndexSet &
DoFHandler<dim, spacedim>::locally_owned_mg_dofs(const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(),
         ExcMessage("The given level index exceeds the number of levels "
                    "present in the triangulation"));
  Assert(
    mg_number_cache.size() == this->get_triangulation().n_global_levels(),
    ExcMessage(
      "The level dofs are not set up properly! Did you call distribute_mg_dofs()?"));
  return mg_number_cache[level].locally_owned_dofs;
}



template <int dim, int spacedim>
const std::vector<types::global_dof_index> &
DoFHandler<dim, spacedim>::n_locally_owned_dofs_per_processor() const
{
  if (number_cache.n_locally_owned_dofs_per_processor.empty() &&
      number_cache.n_global_dofs > 0)
    {
      const_cast<dealii::internal::DoFHandlerImplementation::NumberCache &>(
        number_cache)
        .n_locally_owned_dofs_per_processor =
        compute_n_locally_owned_dofs_per_processor();
    }
  return number_cache.n_locally_owned_dofs_per_processor;
}



template <int dim, int spacedim>
const std::vector<IndexSet> &
DoFHandler<dim, spacedim>::locally_owned_dofs_per_processor() const
{
  if (number_cache.locally_owned_dofs_per_processor.empty() &&
      number_cache.n_global_dofs > 0)
    {
      const_cast<dealii::internal::DoFHandlerImplementation::NumberCache &>(
        number_cache)
        .locally_owned_dofs_per_processor =
        compute_locally_owned_dofs_per_processor();
    }
  return number_cache.locally_owned_dofs_per_processor;
}



template <int dim, int spacedim>
const std::vector<IndexSet> &
DoFHandler<dim, spacedim>::locally_owned_mg_dofs_per_processor(
  const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(),
         ExcMessage("The given level index exceeds the number of levels "
                    "present in the triangulation"));
  Assert(
    mg_number_cache.size() == this->get_triangulation().n_global_levels(),
    ExcMessage(
      "The level dofs are not set up properly! Did you call distribute_mg_dofs()?"));
  if (mg_number_cache[level].locally_owned_dofs_per_processor.empty() &&
      mg_number_cache[level].n_global_dofs > 0)
    {
      const_cast<dealii::internal::DoFHandlerImplementation::NumberCache &>(
        mg_number_cache[level])
        .locally_owned_dofs_per_processor =
        compute_locally_owned_mg_dofs_per_processor(level);
    }
  return mg_number_cache[level].locally_owned_dofs_per_processor;
}



template <int dim, int spacedim>
std::vector<types::global_dof_index>
DoFHandler<dim, spacedim>::compute_n_locally_owned_dofs_per_processor() const
{
  const parallel::TriangulationBase<dim, spacedim> *tr =
    (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
      &this->get_triangulation()));
  if (tr != nullptr)
    return number_cache.get_n_locally_owned_dofs_per_processor(
      tr->get_communicator());
  else
    return number_cache.get_n_locally_owned_dofs_per_processor(MPI_COMM_SELF);
}



template <int dim, int spacedim>
std::vector<IndexSet>
DoFHandler<dim, spacedim>::compute_locally_owned_dofs_per_processor() const
{
  const parallel::TriangulationBase<dim, spacedim> *tr =
    (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
      &this->get_triangulation()));
  if (tr != nullptr)
    return number_cache.get_locally_owned_dofs_per_processor(
      tr->get_communicator());
  else
    return number_cache.get_locally_owned_dofs_per_processor(MPI_COMM_SELF);
}



template <int dim, int spacedim>
std::vector<IndexSet>
DoFHandler<dim, spacedim>::compute_locally_owned_mg_dofs_per_processor(
  const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(),
         ExcMessage("The given level index exceeds the number of levels "
                    "present in the triangulation"));
  Assert(
    mg_number_cache.size() == this->get_triangulation().n_global_levels(),
    ExcMessage(
      "The level dofs are not set up properly! Did you call distribute_mg_dofs()?"));
  const parallel::TriangulationBase<dim, spacedim> *tr =
    (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
      &this->get_triangulation()));
  if (tr != nullptr)
    return mg_number_cache[level].get_locally_owned_dofs_per_processor(
      tr->get_communicator());
  else
    return mg_number_cache[level].get_locally_owned_dofs_per_processor(
      MPI_COMM_SELF);
}



template <int dim, int spacedim>
inline const FiniteElement<dim, spacedim> &
DoFHandler<dim, spacedim>::get_fe(const unsigned int index) const
{
  (void)index;
  Assert(index == 0,
         ExcMessage(
           "There is only one FiniteElement stored. The index must be zero!"));
  return get_fe_collection()[0];
}



template <int dim, int spacedim>
inline const hp::FECollection<dim, spacedim> &
DoFHandler<dim, spacedim>::get_fe_collection() const
{
  Assert(
    fe_collection.size() > 0,
    ExcMessage(
      "You are trying to access the DoFHandler's FECollection object before it has been initialized."));
  return fe_collection;
}



template <int dim, int spacedim>
inline const Triangulation<dim, spacedim> &
DoFHandler<dim, spacedim>::get_triangulation() const
{
  Assert(tria != nullptr,
         ExcMessage("This DoFHandler object has not been associated "
                    "with a triangulation."));
  return *tria;
}



template <int dim, int spacedim>
inline const BlockInfo &
DoFHandler<dim, spacedim>::block_info() const
{
  return block_info_object;
}



template <int dim, int spacedim>
template <typename number>
types::global_dof_index
DoFHandler<dim, spacedim>::n_boundary_dofs(
  const std::map<types::boundary_id, const Function<spacedim, number> *>
    &boundary_ids) const
{
  // extract the set of boundary ids and forget about the function object
  // pointers
  std::set<types::boundary_id> boundary_ids_only;
  for (typename std::map<types::boundary_id,
                         const Function<spacedim, number> *>::const_iterator p =
         boundary_ids.begin();
       p != boundary_ids.end();
       ++p)
    boundary_ids_only.insert(p->first);

  // then just hand everything over to the other function that does the work
  return n_boundary_dofs(boundary_ids_only);
}



namespace internal
{
  /**
   * Return a string representing the dynamic type of the given argument.
   * This is basically the same what typeid(...).name() does, but it turns out
   * this is broken on Intel 13+.
   *
   * Defined in dof_handler.cc.
   */
  template <int dim, int spacedim>
  std::string
  policy_to_string(const dealii::internal::DoFHandlerImplementation::Policy::
                     PolicyBase<dim, spacedim> &policy);
} // namespace internal



template <int dim, int spacedim>
template <class Archive>
void
DoFHandler<dim, spacedim>::save(Archive &ar, const unsigned int) const
{
  ar &block_info_object;
  ar &vertex_dofs;
  ar &number_cache;

  // some versions of gcc have trouble with loading vectors of
  // std::unique_ptr objects because std::unique_ptr does not
  // have a copy constructor. do it one level at a time
  unsigned int n_levels = levels.size();
  ar &         n_levels;
  for (unsigned int i = 0; i < levels.size(); ++i)
    ar &levels[i];

  // boost dereferences a nullptr when serializing a nullptr
  // at least up to 1.65.1. This causes problems with clang-5.
  // Therefore, work around it.
  bool faces_is_nullptr = (faces.get() == nullptr);
  ar & faces_is_nullptr;
  if (!faces_is_nullptr)
    ar &faces;

  // write out the number of triangulation cells and later check during
  // loading that this number is indeed correct; same with something that
  // identifies the FE and the policy
  unsigned int n_cells     = tria->n_cells();
  std::string  fe_name     = this->get_fe(0).get_name();
  std::string  policy_name = internal::policy_to_string(*policy);

  ar &n_cells &fe_name &policy_name;
}



template <int dim, int spacedim>
template <class Archive>
void
DoFHandler<dim, spacedim>::load(Archive &ar, const unsigned int)
{
  ar &block_info_object;
  ar &vertex_dofs;
  ar &number_cache;

  // boost::serialization can restore pointers just fine, but if the
  // pointer object still points to something useful, that object is not
  // destroyed and we end up with a memory leak. consequently, first delete
  // previous content before re-loading stuff
  levels.clear();
  faces.reset();

  // some versions of gcc have trouble with loading vectors of
  // std::unique_ptr objects because std::unique_ptr does not
  // have a copy constructor. do it one level at a time
  unsigned int size;
  ar &         size;
  levels.resize(size);
  for (unsigned int i = 0; i < levels.size(); ++i)
    {
      std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<dim>> level;
      ar &                                                               level;
      levels[i] = std::move(level);
    }

  // Workaround for nullptr, see in save().
  bool faces_is_nullptr = true;
  ar & faces_is_nullptr;
  if (!faces_is_nullptr)
    ar &faces;

  // these are the checks that correspond to the last block in the save()
  // function
  unsigned int n_cells;
  std::string  fe_name;
  std::string  policy_name;

  ar &n_cells &fe_name &policy_name;

  AssertThrow(n_cells == tria->n_cells(),
              ExcMessage(
                "The object being loaded into does not match the triangulation "
                "that has been stored previously."));
  AssertThrow(
    fe_name == this->get_fe(0).get_name(),
    ExcMessage(
      "The finite element associated with this DoFHandler does not match "
      "the one that was associated with the DoFHandler previously stored."));
  AssertThrow(policy_name == internal::policy_to_string(*policy),
              ExcMessage(
                "The policy currently associated with this DoFHandler (" +
                internal::policy_to_string(*policy) +
                ") does not match the one that was associated with the "
                "DoFHandler previously stored (" +
                policy_name + ")."));
}



template <int dim, int spacedim>
inline types::global_dof_index
DoFHandler<dim, spacedim>::MGVertexDoFs::get_index(
  const unsigned int level,
  const unsigned int dof_number,
  const unsigned int dofs_per_vertex) const
{
  Assert((level >= coarsest_level) && (level <= finest_level),
         ExcInvalidLevel(level));
  return indices[dofs_per_vertex * (level - coarsest_level) + dof_number];
}



template <int dim, int spacedim>
inline void
DoFHandler<dim, spacedim>::MGVertexDoFs::set_index(
  const unsigned int            level,
  const unsigned int            dof_number,
  const unsigned int            dofs_per_vertex,
  const types::global_dof_index index)
{
  Assert((level >= coarsest_level) && (level <= finest_level),
         ExcInvalidLevel(level));
  indices[dofs_per_vertex * (level - coarsest_level) + dof_number] = index;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

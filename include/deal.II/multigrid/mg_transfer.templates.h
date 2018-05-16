// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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


#ifndef dealii_mg_transfer_templates_h
#define dealii_mg_transfer_templates_h

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector_view.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/distributed/tria.h>

#include <algorithm>

// Here you can turn on some cout statements and MPI Barriers for debugging:
//#define DEBUG_OUTPUT

DEAL_II_NAMESPACE_OPEN


namespace
{
  /**
   * Adjust vectors on all levels to correct size.  Here, we just count the
   * numbers of degrees of freedom on each level and @p reinit each level
   * vector to this length. For compatibility reasons with the next function
   * the target_component is added here but is not used.
   */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector (const dealii::DoFHandler<dim,spacedim> &mg_dof,
                 const std::vector<unsigned int> &,
                 MGLevelObject<dealii::Vector<number> > &v)
  {
    for (unsigned int level=v.min_level();
         level<=v.max_level(); ++level)
      {
        unsigned int n = mg_dof.n_dofs (level);
        v[level].reinit(n);
      }

  }

  /**
   * Adjust vectors on all levels to correct size.  Here, we just count the
   * numbers of degrees of freedom on each level and @p reinit each level
   * vector to this length. The target_component is handed to
   * MGTools::count_dofs_per_block. See for documentation there.
   */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector (const dealii::DoFHandler<dim,spacedim> &mg_dof,
                 std::vector<unsigned int> target_component,
                 MGLevelObject<BlockVector<number> > &v)
  {
    const unsigned int n_blocks = mg_dof.get_fe().n_blocks();
    if (target_component.size()==0)
      {
        target_component.resize(n_blocks);
        for (unsigned int i=0; i<n_blocks; ++i)
          target_component[i] = i;
      }
    Assert(target_component.size()==n_blocks,
           ExcDimensionMismatch(target_component.size(),n_blocks));
    const unsigned int max_block
      = *std::max_element (target_component.begin(),
                           target_component.end());
    const unsigned int n_target_blocks = max_block + 1;

    std::vector<std::vector<types::global_dof_index> >
    ndofs(mg_dof.get_triangulation().n_levels(),
          std::vector<types::global_dof_index>(n_target_blocks));
    MGTools::count_dofs_per_block (mg_dof, ndofs, target_component);

    for (unsigned int level=v.min_level();
         level<=v.max_level(); ++level)
      {
        v[level].reinit(n_target_blocks);
        for (unsigned int b=0; b<n_target_blocks; ++b)
          v[level].block(b).reinit(ndofs[level][b]);
        v[level].collect_sizes();
      }
  }

  /**
   * Adjust vectors on all levels to correct size.  Here, we just count the
   * numbers of degrees of freedom on each level and @p reinit each level
   * vector to this length.
   */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector (const dealii::DoFHandler<dim,spacedim> &mg_dof,
                 const std::vector<unsigned int> &,
                 MGLevelObject<LinearAlgebra::distributed::Vector<number> > &v)
  {
    const parallel::Triangulation<dim,spacedim> *tria =
      (dynamic_cast<const parallel::Triangulation<dim,spacedim>*>
       (&mg_dof.get_triangulation()));

    for (unsigned int level=v.min_level(); level<=v.max_level(); ++level)
      {
        if (v[level].size() != mg_dof.locally_owned_mg_dofs(level).size() ||
            v[level].local_size() != mg_dof.locally_owned_mg_dofs(level).n_elements())
          v[level].reinit(mg_dof.locally_owned_mg_dofs(level), tria != nullptr ?
                          tria->get_communicator() : MPI_COMM_SELF);
        else
          v[level] = 0.;
      }
  }


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Adjust vectors on all levels to correct size.  Here, we just count the
   * numbers of degrees of freedom on each level and @p reinit each level
   * vector to this length.
   */
  template <int dim, int spacedim>
  void
  reinit_vector (const dealii::DoFHandler<dim,spacedim> &mg_dof,
                 const std::vector<unsigned int> &,
                 MGLevelObject<TrilinosWrappers::MPI::Vector> &v)
  {
    const dealii::parallel::Triangulation<dim,spacedim> *tria =
      (dynamic_cast<const parallel::Triangulation<dim,spacedim>*>
       (&mg_dof.get_triangulation()));
    AssertThrow(tria!=nullptr, ExcMessage("multigrid with Trilinos vectors only works with a parallel Triangulation!"));

    for (unsigned int level=v.min_level();
         level<=v.max_level(); ++level)
      {
        v[level].reinit(mg_dof.locally_owned_mg_dofs(level), tria->get_communicator());
      }
  }
#endif
}



/* ------------------ MGLevelGlobalTransfer<VectorType> ----------------- */


namespace internal
{
  // generic copy function of two different vectors -> need to access each
  // individual entry
  template <typename T, typename V>
  void
  copy_vector (const std::vector<std::pair<types::global_dof_index,types::global_dof_index> > &copy_indices,
               const T &src,
               V &dst)
  {
    // we should have i->second == i->first, therefore we can use the same
    // function for both copying to mg as well as copying from mg
    for (std::vector<std::pair<types::global_dof_index, types::global_dof_index> >::
         const_iterator i = copy_indices.begin(); i != copy_indices.end(); ++i)
      dst(i->first) = src(i->first);
    dst.compress(VectorOperation::insert);
  }

  // specialized copy function for the same vector
  template <typename T>
  void
  copy_vector (const std::vector<std::pair<types::global_dof_index,types::global_dof_index> > &,
               const T &src,
               T &dst)
  {
    dst = src;
  }
}


template <typename VectorType>
template <int dim, class InVector, int spacedim>
void
MGLevelGlobalTransfer<VectorType>::copy_to_mg
(const DoFHandler<dim,spacedim> &mg_dof_handler,
 MGLevelObject<VectorType>      &dst,
 const InVector                 &src) const
{
  AssertIndexRange(dst.max_level(), mg_dof_handler.get_triangulation().n_global_levels());
  AssertIndexRange(dst.min_level(), dst.max_level()+1);
  reinit_vector(mg_dof_handler, component_to_block_map, dst);
#ifdef DEBUG_OUTPUT
  std::cout << "copy_to_mg src " << src.l2_norm() << std::endl;
  int ierr = MPI_Barrier(MPI_COMM_WORLD);
  AssertThrowMPI(ierr);
#endif

  if (perform_plain_copy)
    {
      // if the finest multigrid level covers the whole domain (i.e., no
      // adaptive refinement) and the numbering of the finest level DoFs and
      // the global DoFs are the same, we can do a plain copy
      AssertDimension(dst[dst.max_level()].size(), src.size());
      internal::copy_vector(copy_indices[dst.max_level()], src, dst[dst.max_level()]);
      return;
    }

  for (unsigned int level=dst.max_level()+1; level != dst.min_level(); )
    {
      --level;
#ifdef DEBUG_OUTPUT
      ierr = MPI_Barrier(MPI_COMM_WORLD);
      AssertThrowMPI(ierr);
#endif

      typedef std::vector<std::pair<types::global_dof_index, types::global_dof_index> >::const_iterator dof_pair_iterator;
      VectorType &dst_level = dst[level];

      // first copy local unknowns
      for (dof_pair_iterator i = copy_indices[level].begin();
           i != copy_indices[level].end(); ++i)
        dst_level(i->second) = src(i->first);

      // Do the same for the indices where the global index is local, but the
      // local index is not
      for (dof_pair_iterator i = copy_indices_global_mine[level].begin();
           i != copy_indices_global_mine[level].end(); ++i)
        dst_level(i->second) = src(i->first);

      dst_level.compress(VectorOperation::insert);

#ifdef DEBUG_OUTPUT
      ierr = MPI_Barrier(MPI_COMM_WORLD);
      AssertThrowMPI(ierr);
      std::cout << "copy_to_mg dst " << level << " " << dst_level.l2_norm() << std::endl;
#endif
    }
}



template <typename VectorType>
template <int dim, class OutVector, int spacedim>
void
MGLevelGlobalTransfer<VectorType>::copy_from_mg
(const DoFHandler<dim,spacedim>  &mg_dof_handler,
 OutVector                       &dst,
 const MGLevelObject<VectorType> &src) const
{
  (void)mg_dof_handler;
  AssertIndexRange(src.max_level(), mg_dof_handler.get_triangulation().n_global_levels());
  AssertIndexRange(src.min_level(), src.max_level()+1);
  if (perform_plain_copy)
    {
      AssertDimension(dst.size(), src[src.max_level()].size());
      internal::copy_vector(copy_indices[src.max_level()], src[src.max_level()], dst);
      return;
    }

  // For non-DG: degrees of freedom in the refinement face may need special
  // attention, since they belong to the coarse level, but have fine level
  // basis functions
  dst = 0;
  for (unsigned int level=src.min_level(); level<=src.max_level(); ++level)
    {
#ifdef DEBUG_OUTPUT
      int ierr = MPI_Barrier(MPI_COMM_WORLD);
      AssertThrowMPI(ierr);
      std::cout << "copy_from_mg src " << level << " " << src[level].l2_norm() << std::endl;
      ierr = MPI_Barrier(MPI_COMM_WORLD);
      AssertThrowMPI(ierr);
#endif

      typedef std::vector<std::pair<types::global_dof_index, types::global_dof_index> >::const_iterator dof_pair_iterator;
      const VectorType &src_level = src[level];

      // First copy all indices local to this process
      for (dof_pair_iterator i = copy_indices[level].begin();
           i != copy_indices[level].end(); ++i)
        dst(i->first) = src_level(i->second);

      // Do the same for the indices where the level index is local, but the
      // global index is not
      for (dof_pair_iterator i = copy_indices_level_mine[level].begin();
           i != copy_indices_level_mine[level].end(); ++i)
        dst(i->first) = src_level(i->second);

#ifdef DEBUG_OUTPUT
      {
        dst.compress(VectorOperation::insert);
        ierr = MPI_Barrier(MPI_COMM_WORLD);
        AssertThrowMPI(ierr);
        std::cout << "copy_from_mg level=" << level << " " << dst.l2_norm() << std::endl;
      }
#endif
    }
  dst.compress(VectorOperation::insert);
#ifdef DEBUG_OUTPUT
  const int ierr = MPI_Barrier(MPI_COMM_WORLD);
  AssertThrowMPI(ierr);
  std::cout << "copy_from_mg " << dst.l2_norm() << std::endl;
#endif
}



template <typename VectorType>
template <int dim, class OutVector, int spacedim>
void
MGLevelGlobalTransfer<VectorType>::copy_from_mg_add
(const DoFHandler<dim,spacedim>  &/*mg_dof_handler*/,
 OutVector                       &dst,
 const MGLevelObject<VectorType> &src) const
{
  // For non-DG: degrees of freedom in the refinement face may need special
  // attention, since they belong to the coarse level, but have fine level
  // basis functions
  for (unsigned int level=src.min_level(); level<=src.max_level(); ++level)
    {
      typedef std::vector<std::pair<types::global_dof_index, types::global_dof_index> >::const_iterator dof_pair_iterator;
      const VectorType &src_level = src[level];

      // First add all indices local to this process
      for (dof_pair_iterator i = copy_indices[level].begin();
           i != copy_indices[level].end(); ++i)
        dst(i->first) += src_level(i->second);

      // Do the same for the indices where the level index is local, but the
      // global index is not
      for (dof_pair_iterator i = copy_indices_level_mine[level].begin();
           i != copy_indices_level_mine[level].end(); ++i)
        dst(i->first) += src_level(i->second);
    }
  dst.compress(VectorOperation::add);
}



template <typename VectorType>
void
MGLevelGlobalTransfer<VectorType>::
set_component_to_block_map (const std::vector<unsigned int> &map)
{
  component_to_block_map = map;
}



/* --------- MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector> ------- */

template <typename Number>
template <int dim, typename Number2, int spacedim>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::copy_to_mg
(const DoFHandler<dim,spacedim>                        &mg_dof_handler,
 MGLevelObject<LinearAlgebra::distributed::Vector<Number> > &dst,
 const LinearAlgebra::distributed::Vector<Number2>          &src) const
{
  copy_to_mg (mg_dof_handler,
              dst,
              src,
              false);
}


template <typename Number>
template <int dim, typename Number2, int spacedim>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::copy_to_mg
(const DoFHandler<dim,spacedim>                             &mg_dof_handler,
 MGLevelObject<LinearAlgebra::distributed::Vector<Number> > &dst,
 const LinearAlgebra::distributed::Vector<Number2>          &src,
 const bool solution_transfer) const
{
  LinearAlgebra::distributed::Vector<Number> &this_ghosted_global_vector =
    solution_transfer ? solution_ghosted_global_vector : ghosted_global_vector;
  const std::vector<std::vector<std::pair<unsigned int, unsigned int> > > &this_copy_indices =
    solution_transfer ? solution_copy_indices : copy_indices;
  const std::vector<std::vector<std::pair<unsigned int, unsigned int> > > &this_copy_indices_level_mine =
    solution_transfer ? solution_copy_indices_level_mine : copy_indices_level_mine;

  AssertIndexRange(dst.max_level(), mg_dof_handler.get_triangulation().n_global_levels());
  AssertIndexRange(dst.min_level(), dst.max_level()+1);
  reinit_vector(mg_dof_handler, component_to_block_map, dst);

  if (perform_plain_copy)
    {
      // In this case, we can simply copy the local range.
      AssertDimension(dst[dst.max_level()].local_size(), src.local_size());
      dst[dst.max_level()].copy_locally_owned_data_from(src);
      return;
    }
  else if (perform_renumbered_plain_copy)
    {
      Assert(copy_indices_level_mine.back().empty(), ExcInternalError());
      LinearAlgebra::distributed::Vector<Number> &dst_level = dst[dst.max_level()];
      for (std::pair<unsigned int,unsigned int> i : this_copy_indices.back())
        dst_level.local_element(i.second) = src.local_element(i.first);
      return;
    }

  // the ghosted vector should already have the correct local size (but
  // different parallel layout)
  AssertDimension(ghosted_global_vector.local_size(), src.local_size());

  // copy the source vector to the temporary vector that we hold for the
  // purpose of data exchange
  this_ghosted_global_vector = src;
  this_ghosted_global_vector.update_ghost_values();

  for (unsigned int level=dst.max_level()+1; level != dst.min_level();)
    {
      --level;

      typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator dof_pair_iterator;
      LinearAlgebra::distributed::Vector<Number> &dst_level = dst[level];

      // first copy local unknowns
      for (dof_pair_iterator i = this_copy_indices[level].begin();
           i != this_copy_indices[level].end(); ++i)
        dst_level.local_element(i->second) = this_ghosted_global_vector.local_element(i->first);

      // Do the same for the indices where the level index is local, but the
      // global index is not
      for (dof_pair_iterator i = this_copy_indices_level_mine[level].begin();
           i != this_copy_indices_level_mine[level].end(); ++i)
        dst_level.local_element(i->second) = this_ghosted_global_vector.local_element(i->first);

      dst_level.compress(VectorOperation::insert);
    }
}



template <typename Number>
template <int dim, typename Number2, int spacedim>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::copy_from_mg
(const DoFHandler<dim,spacedim>                                   &mg_dof_handler,
 LinearAlgebra::distributed::Vector<Number2>                      &dst,
 const MGLevelObject<LinearAlgebra::distributed::Vector<Number> > &src) const
{
  (void)mg_dof_handler;
  AssertIndexRange(src.max_level(), mg_dof_handler.get_triangulation().n_global_levels());
  AssertIndexRange(src.min_level(), src.max_level()+1);
  if (perform_plain_copy)
    {
      // In this case, we can simply copy the local range. To avoid having
      // stray data in ghost entries of the destination, make sure to clear
      // them here.
      dst.zero_out_ghosts();
      AssertDimension(src[src.max_level()].local_size(), dst.local_size());
      dst.copy_locally_owned_data_from(src[src.max_level()]);
      return;
    }
  else if (perform_renumbered_plain_copy)
    {
      Assert(copy_indices_global_mine.back().empty(), ExcInternalError());
      const LinearAlgebra::distributed::Vector<Number> &src_level = src[src.max_level()];
      dst.zero_out_ghosts();
      for (std::pair<unsigned int,unsigned int> i : copy_indices.back())
        dst.local_element(i.first) = src_level.local_element(i.second);
      return;
    }

  // For non-DG: degrees of freedom in the refinement face may need special
  // attention, since they belong to the coarse level, but have fine level
  // basis functions

  dst = 0;
  for (unsigned int level=src.min_level(); level<=src.max_level(); ++level)
    {
      typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator dof_pair_iterator;

      // the ghosted vector should already have the correct local size (but
      // different parallel layout)
      AssertDimension(ghosted_level_vector[level].local_size(),
                      src[level].local_size());

      // the first time around, we copy the source vector to the temporary
      // vector that we hold for the purpose of data exchange
      LinearAlgebra::distributed::Vector<Number> &ghosted_vector =
        ghosted_level_vector[level];
      ghosted_vector = src[level];
      ghosted_vector.update_ghost_values();

      // first copy local unknowns
      for (dof_pair_iterator i = copy_indices[level].begin();
           i != copy_indices[level].end(); ++i)
        dst.local_element(i->first) = ghosted_vector.local_element(i->second);

      // Do the same for the indices where the level index is local, but the
      // global index is not
      for (dof_pair_iterator i = copy_indices_global_mine[level].begin();
           i != copy_indices_global_mine[level].end(); ++i)
        dst.local_element(i->first) = ghosted_vector.local_element(i->second);
    }
  dst.compress(VectorOperation::insert);
}



template <typename Number>
template <int dim, typename Number2, int spacedim>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::copy_from_mg_add
(const DoFHandler<dim,spacedim>                                   &/*mg_dof_handler*/,
 LinearAlgebra::distributed::Vector<Number2>                      &dst,
 const MGLevelObject<LinearAlgebra::distributed::Vector<Number> > &src) const
{
  // For non-DG: degrees of freedom in the refinement face may need special
  // attention, since they belong to the coarse level, but have fine level
  // basis functions

  dst.zero_out_ghosts();
  for (unsigned int level=src.min_level(); level<=src.max_level(); ++level)
    {
      typedef std::vector<std::pair<unsigned int, unsigned int> >::const_iterator dof_pair_iterator;

      // the ghosted vector should already have the correct local size (but
      // different parallel layout)
      AssertDimension(ghosted_level_vector[level].local_size(),
                      src[level].local_size());

      // the first time around, we copy the source vector to the temporary
      // vector that we hold for the purpose of data exchange
      LinearAlgebra::distributed::Vector<Number> &ghosted_vector =
        ghosted_level_vector[level];
      ghosted_vector = src[level];
      ghosted_vector.update_ghost_values();

      // first add local unknowns
      for (dof_pair_iterator i= copy_indices[level].begin();
           i != copy_indices[level].end(); ++i)
        dst.local_element(i->first) += ghosted_vector.local_element(i->second);

      // Do the same for the indices where the level index is local, but the
      // global index is not
      for (dof_pair_iterator i= copy_indices_global_mine[level].begin();
           i != copy_indices_global_mine[level].end(); ++i)
        dst.local_element(i->first) += ghosted_vector.local_element(i->second);
    }
  dst.compress(VectorOperation::add);
}



template <typename Number>
void
MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >::
set_component_to_block_map (const std::vector<unsigned int> &map)
{
  component_to_block_map = map;
}


DEAL_II_NAMESPACE_CLOSE

#endif

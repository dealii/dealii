// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#ifndef __deal2__block_info_h
#define __deal2__block_info_h

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/block_indices.h>

#include <iomanip>

DEAL_II_NAMESPACE_OPEN

// Forward declarations

template <int dim, int spacedim> class DoFHandler;
namespace hp
{
  template <int dim, int spacedim> class DoFHandler;
}


/**
 * @brief A small class collecting the different BlockIndices involved in
 * global, multilevel and local computations.
 *
 * Once a DoFHandler has been initialized with an FESystem, a data
 * object of type BlockInfo (accessed by DoFHandler::block_info() ) is
 * filled, which reflects the block structure of the degrees of
 * freedom.
 *
 * BlockInfo consists of several BlockIndices objects. The member
 * global() reflects the block structure of the system on the active
 * cell level, usually referred to as the global system. As soon as
 * DoFHandler::distribute_dofs() has been called, the function
 * BlockIndices::block_size() in global() will return the correct
 * sizes of each block. After DoFRenumbering::block_wise(),
 * BlockIndices::block_start() will return the start index for each of
 * the blocks.
 *
 * When a DoFHandler with levels is used, the same structure is automatically
 * generated for each level. The level blocks can be accessed through
 * level().
 *
 * Finally, there are local() BlockIndices, which describe the block
 * structure on a single cell. This is used for instance by
 * MeshWorker::Assembler::MatrixLocalBlocksToGlobalBlocks. The local
 * indices are not filled automatically, since they change the
 * behavior of the MeshWorker::Assembler classes relying on
 * BlockInfo. They must be initialized by hand through initialize_local().
 *
 * <h3>Usage</h3>
 *
 * The most common usage for this object is initializing vectors as in
 * the following code:
 *
 * <code>
 * MGDofHandler<dim> dof_handler(triangulation);
 * dof_handler.distribute_dofs(fesystem);
 * DoFRenumbering::block_wise(dof_handler);
 *
 * BlockVector<double> solution(dof_handler.block_info().global());
 *
 * MGLevelObject<BlockVector<double> > mg_vector(0, triangulation.n_levels()-1);
 * for (unsigned int i=0;i<triangulation.n_levels();++i)
 *   mg_vector[i].reinit(dof_handler.block_info().level(i));
 *</code>
 * In this example, <tt>solution</tt> obtains the block structure needed to
 * represent a finite element function on the DoFHandler. Similarly,
 * all levels of <tt>mg_vector</tt> will have the block structure
 * needed on that level.
 *
 * @todo Extend the functions local() and renumber() to the concept to
 * hpDoFHandler.
 *
 * @ingroup dofs
 * @author Guido Kanschat, 2009
 */
class BlockInfo : public Subscriptor
{
public:
  /**
   * @brief Fill the object with values
   * describing block structure
   * of the DoFHandler.
   *
  * By default, this function will
  * attempt to initialize whatever
  * is possible. If active dofs
  * have been assigned int the
  * DoFHandler argument, they
  * BlockIndices for those will be
  * generated. The same for level
  * dofs.
  *
  * This default behavior can be
  * overridden by the two
  * parameters, which can switch
  * off active dofs or level dofs.
  *
   * This function will also clear
   * the local() indices.
   */
  template <int dim, int spacedim>
  void initialize(const DoFHandler<dim, spacedim> &, bool levels_only = false, bool active_only = false);

  /**
   * @brief Initialize block structure
   * on cells and compute
   * renumbering between cell
   * dofs and block cell dofs.
   */
  template <int dim, int spacedim>
  void initialize_local(const DoFHandler<dim, spacedim> &);

  /**
   * Access the BlockIndices
   * structure of the global
   * system.
   */
  const BlockIndices &global() const;

  /**
   * Access BlockIndices for the
   * local system on a cell.
   */
  const BlockIndices &local() const;

  /**
   * Access the BlockIndices
   * structure of a level in the
   * multilevel hierarchy.
   */
  const BlockIndices &level(unsigned int level) const;

  /**
   * Return the index after local
   * renumbering.
   *
   * The input of this function is
   * an index between zero and the
   * number of dofs per cell,
   * numbered in local block
   * ordering, that is first all
   * indices of the first system
   * block, then all of the second
   * block and so forth. The
   * function then outputs the index
   * in the standard local
   * numbering of DoFAccessor.
   */
  types::global_dof_index renumber (const unsigned int i) const;

  /**
   * The number of base elements.
   */
  unsigned int n_base_elements() const;

  /**
   * Return the base element of
   * this index.
   */
  unsigned int base_element (const unsigned int i) const;

  /**
   * Write a summary of the block
   * structure to the stream.
   */
  template <class OS>
  void
  print(OS &stream) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Read or write the data of this object to or
   * from a stream for the purpose of serialization
   */
  template <class Archive>
  void serialize (Archive &ar,
                  const unsigned int version);

private:
  /**
   * @brief The block structure
   * of the global system.
   */
  BlockIndices bi_global;
  /**
   * @brief The multilevel block structure.
   */
  std::vector<BlockIndices> levels;

  /**
   * @brief The block structure
   * of the cell systems.
   */
  BlockIndices bi_local;

  /**
   * The base element associated
   * with each block.
   */
  std::vector<unsigned int> base_elements;

  /**
   * A vector containing the
   * renumbering from the
   * standard order of degrees of
   * freedom on a cell to a
   * component wise
   * ordering. Filled by
   * initialize().
   */
  std::vector<types::global_dof_index> local_renumbering;
};



//----------------------------------------------------------------------//

inline
const BlockIndices &
BlockInfo::global() const
{
  return bi_global;
}


inline
const BlockIndices &
BlockInfo::local() const
{
  return bi_local;
}


inline
const BlockIndices &
BlockInfo::level (const unsigned int l) const
{
  AssertIndexRange(l, levels.size());
  return levels[l];
}


inline
types::global_dof_index BlockInfo::renumber (const unsigned int i) const
{
  AssertIndexRange(i, static_cast<unsigned int>(local_renumbering.size()));
  return local_renumbering[i];
}


inline
unsigned int
BlockInfo::base_element (const unsigned int i) const
{
  AssertIndexRange(i, base_elements.size());

  return base_elements[i];
}


inline
unsigned int
BlockInfo::n_base_elements() const
{
  return base_elements.size();
}



template <class OS>
inline
void
BlockInfo::print (OS &os) const
{
  os << "global   dofs " << std::setw(5) << global().total_size() << " blocks";
  for (unsigned int i=0; i<global().size(); ++i)
    os << ' ' << std::setw(5) << global().block_size(i);
  os << std::endl;

  if (local().size() == 0)
    {
      os << "local dofs not initialized" << std::endl;
    }
  else
    {
      os << "local    dofs " << std::setw(5) << local().total_size() << " blocks";
      for (unsigned int i=0; i<local().size(); ++i)
        os << ' '  << std::setw(5) << local().block_size(i);
      os << std::endl;
    }

  for (unsigned int l=0; l<levels.size(); ++l)
    {
      os << "level " << std::setw(2) << l << " dofs " << std::setw(5) << level(l).total_size() << " blocks";
      for (unsigned int i=0; i<level(l).size(); ++i)
        os << ' '  << std::setw(5) << level(l).block_size(i);
      os << std::endl;
    }
}


inline
std::size_t
BlockInfo::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (bi_global) +
          MemoryConsumption::memory_consumption (levels) +
          MemoryConsumption::memory_consumption (bi_local) +
          MemoryConsumption::memory_consumption (base_elements));
}


template <class Archive>
void BlockInfo::serialize (Archive &ar,
                           const unsigned int version)
{
  ar &bi_global;
  ar &levels;
  ar &bi_local;
  ar &base_elements;
  ar &local_renumbering;
}


DEAL_II_NAMESPACE_CLOSE

#endif

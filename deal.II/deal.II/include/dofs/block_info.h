//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_info_h
#define __deal2__block_info_h

#include <base/subscriptor.h>
#include <lac/block_indices.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations

template <int dim, int spacedim> class DoFHandler;
template <int dim, int spacedim> class MGDoFHandler;
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
 * BlockInfo consists of deveral BlockIndices objects. The member
 * global() reflects the block structure of the system on the active
 * cell level, usually referred to as the global system. As soon as
 * DoFHandler::distribute_dofs() has been called, the function
 * BlockIndices::block_size() in global() will return the correct
 * sizes of each block. After DoFRenumbering::block_wise(),
 * BlockIndices::block_start() will return the start index for each of
 * the blocks.
 *
 * When an MGDoFHandler is used, the same structure is automatically
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
 * @todo Extend the functions local() and renumber() to the concept to
 * hpDoFHandler.
 *
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
				      * This function will also clear
				      * the local() indices.
				      */
    template <int dim, int spacedim>
    void initialize(const DoFHandler<dim, spacedim>&);
    
				     /**
				      * @brief Fill the object with values
				      * describing level block
				      * structure of the
				      * MGDoFHandler. If
				      * <tt>levels_only</tt> is false,
				      * the other initialize() is
				      * called as well.
				      *
				      * This function will also clear
				      * the local() indices.
				      */
    template <int dim, int spacedim>
    void initialize(const MGDoFHandler<dim, spacedim>&, bool levels_only = false);
    
				     /**
				      * @brief Initialize block structure
				      * on cells and compute
				      * renumbering between cell
				      * dofs and block cell dofs.
				      */
    template <int dim, int spacedim>
    void initialize_local(const DoFHandler<dim, spacedim>&);
    
				     /**
				      * Access the BlockIndices
				      * structure of the global
				      * system.
				      */
    const BlockIndices& global() const;
    
				     /**
				      * Access BlockIndices for the
				      * local system on a cell.
				      */
    const BlockIndices& local() const;
    
				     /**
				      * Access the BlockIndices
				      * structure of a level in the
				      * multilevel hierarchy.
				      */
    const BlockIndices& level(unsigned int level) const;

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
				      * function then outpus the index
				      * in the standard local
				      * numbering of DoFAccessor.
				      */
    unsigned int renumber (unsigned int i) const;

				     /**
				      * The number of base elements.
				      */
    unsigned int n_base_elements() const;
    
				     /**
				      * Return the base element of
				      * this index.
				      */
    unsigned int base_element(unsigned int i) const;
    
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
    std::vector<unsigned int> local_renumbering;
};



//----------------------------------------------------------------------//

inline
const BlockIndices&
BlockInfo::global() const
{
  return bi_global;
}


inline
const BlockIndices&
BlockInfo::local() const
{
  return bi_local;
}


inline
const BlockIndices&
BlockInfo::level(unsigned int l) const
{
  AssertIndexRange(l, levels.size());
  return levels[l];
}


inline
unsigned int BlockInfo::renumber(unsigned int i) const
{
  AssertIndexRange(i, local_renumbering.size());
  return local_renumbering[i];
}


inline
unsigned int
BlockInfo::base_element(unsigned int i) const
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



DEAL_II_NAMESPACE_CLOSE

#endif

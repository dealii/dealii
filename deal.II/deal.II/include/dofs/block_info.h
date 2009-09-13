//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
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
 * In addition to grouping all the block indices into a single
 * instance, this class also provides functions to fill them
 * automatically.
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
				      * structure of the MGDoFHandler.
				      *
				      * This function will also clear
				      * the local() indices.
				      */
    template <int dim, int spacedim>
    void initialize(const MGDoFHandler<dim, spacedim>&);
    
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
    std::vector<unsigned int> base_element;
    
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


DEAL_II_NAMESPACE_CLOSE

#endif

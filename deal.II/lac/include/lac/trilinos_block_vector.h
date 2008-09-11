//---------------------------------------------------------------------------
//    $Id: trilinos_block_vector.h 14783 2007-06-18 14:52:01Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_block_vector_h
#define __deal2__trilinos_block_vector_h


#include <base/config.h>
#include <lac/trilinos_vector.h>
#include <lac/block_indices.h>
#include <lac/block_vector_base.h>
#include <lac/block_vector.h>
#include <lac/exceptions.h>

#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

                                   // forward declaration
template <typename Number> class BlockVector;

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{
                                   // forward declaration
  class BlockVector;


/**
 * An implementation of block vectors based on the vector class
 * implemented in TrilinosWrappers. While the base class provides for
 * most of the interface, this class handles the actual allocation of
 * vectors and provides functions that are specific to the underlying
 * vector type.
 *
 * The model of distribution of data is such that each of the blocks
 * is distributed across all MPI processes named in the MPI
 * communicator. I.e. we don't just distribute the whole vector, but
 * each component. In the constructors and reinit() functions, one
 * therefore not only has to specify the sizes of the individual
 * blocks, but also the number of elements of each of these blocks to
 * be stored on the local process.
 *
 * @ingroup Vectors
 * @ingroup TrilinosWrappers
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
 */
    class BlockVector : public BlockVectorBase<Vector>
    {
      public:
                                       /**
					* Typedef the base class for simpler
					* access to its own typedefs.
					*/
        typedef BlockVectorBase<Vector> BaseClass;
    
                                       /**
					* Typedef the type of the underlying
					* vector.
					*/
        typedef BaseClass::BlockType  BlockType;

                                       /**
					* Import the typedefs from the base
					* class.
					*/
        typedef BaseClass::value_type      value_type;
        typedef BaseClass::pointer         pointer;
        typedef BaseClass::const_pointer   const_pointer;
        typedef BaseClass::reference       reference;
        typedef BaseClass::const_reference const_reference;
        typedef BaseClass::size_type       size_type;
        typedef BaseClass::iterator        iterator;
        typedef BaseClass::const_iterator  const_iterator;

                                       /**
					* Default constructor. Generate an
					* empty vector without any blocks.
					*/
        BlockVector ();
        
                                       /**
					* Constructor. Generate a block
					* vector with as many blocks as
					* there are entries in Input_Maps.
					* Each Epetra_Map already knows
					* the distribution of data among
					* the MPI processes.
					*/
        explicit BlockVector (const std::vector<Epetra_Map> &InputMaps);
    
                                       /**
					* Copy-Constructor. Set all the
					* properties of the parallel vector
					* to those of the given argument and
					* copy the elements.
					*/
        BlockVector (const BlockVector  &V);
    
                                       /**
					* Creates a block vector
					* consisting of
					* <tt>num_blocks</tt>
					* components, but there is no
					* content in the individual
					* components and the user has to
					* fill appropriate data using a
					* reinit of the blocks.
					*/
        BlockVector (const unsigned int num_blocks);
    
                                       /**
					* Destructor. Clears memory
					*/
        ~BlockVector ();

                                       /**
					* Copy operator: fill all
					* components of the vector that
					* are locally stored with the
					* given scalar value.
					*/
        BlockVector &
	operator = (const value_type s);

                                       /**
					* Copy operator for arguments of
					* the same type.
					*/
        BlockVector &
        operator = (const BlockVector &V);

                                       /**
					* Another copy function. This
					* one takes a deal.II block
					* vector and copies it into a
					* TrilinosWrappers block
					* vector. Note that the number
					* of blocks has to be the same
					* in the vector as in the input
					* vector. Use the reinit()
					* command for resizing the
					* BlockVector or for changing
					* the internal structure of the
					* block components.
					*
					* Since Trilinos only works on
					* doubles, this function is
					* limited to accept only one
					* possible number type in the
					* deal.II vector.
					*/
      template <typename Number>
      BlockVector & 
      operator = (const ::dealii::BlockVector<Number> &V);

                                         /**
                                          * Reinitialize the BlockVector to
                                          * contain as many blocks as there 
					  * are Epetra_Maps given in the input
					  * argument, according to the
					  * parallel distribution of the 
					  * individual components described
					  * in the maps.
                                          *
                                          * If <tt>fast==false</tt>, the vector
                                          * is filled with zeros.
                                          */
        void reinit (const std::vector<Epetra_Map> &input_maps);

                                         /**
                                          * Change the dimension to that
                                          * of the vector <tt>V</tt>. The same
                                          * applies as for the other
                                          * reinit() function.
                                          *
                                          * The elements of <tt>V</tt> are not
                                          * copied, i.e.  this function is
                                          * the same as calling <tt>reinit
                                          * (V.size(), fast)</tt>.
                                          *
                                          * Note that you must call this
                                          * (or the other reinit()
                                          * functions) function, rather
                                          * than calling the reinit()
                                          * functions of an individual
                                          * block, to allow the block
                                          * vector to update its caches of
                                          * vector sizes. If you call
                                          * reinit() on one of the
                                          * blocks, then subsequent
                                          * actions on this object may
                                          * yield unpredictable results
                                          * since they may be routed to
                                          * the wrong block.
                                          */
        void reinit (const BlockVector &V,
		     const bool fast=false);

                                         /**
                                          * Change the number of blocks to
					  * <tt>num_blocks</tt>. The individual
					  * blocks will get initialized with
					  * zero size, so it is assumed that
					  * the user resizes the
					  * individual blocks by herself 
					  * in an appropriate way, and
					  * calls <tt>collect_sizes</tt> 
					  * afterwards.
                                          */
        void reinit (const unsigned int num_blocks);

					 /** 
					  * Compresses all the components
					  * after assembling together all
					  * elements.
					  */
	void compress ();

                                         /**
                                          * Swap the contents of this
                                          * vector and the other vector
                                          * <tt>v</tt>. One could do this
                                          * operation with a temporary
                                          * variable and copying over the
                                          * data elements, but this
                                          * function is significantly more
                                          * efficient since it only swaps
                                          * the pointers to the data of
                                          * the two vectors and therefore
                                          * does not need to allocate
                                          * temporary storage and move
                                          * data around.
                                          *
                                          * Limitation: right now this
                                          * function only works if both
                                          * vectors have the same number
                                          * of blocks. If needed, the
                                          * numbers of blocks should be
                                          * exchanged, too.
                                          *
                                          * This function is analog to the
                                          * the swap() function of all C++
                                          * standard containers. Also,
                                          * there is a global function
                                          * swap(u,v) that simply calls
                                          * <tt>u.swap(v)</tt>, again in analogy
                                          * to standard functions.
                                          */
        void swap (BlockVector &v);

                                         /**
                                          * Exception
                                          */
        DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);

                                         /**
                                          * Exception
                                          */
        DeclException0 (ExcNonMatchingBlockVectors);
    };

/*----------------------- Inline functions ----------------------------------*/


    inline
    BlockVector::BlockVector ()
    {}



    inline
    BlockVector::BlockVector (const std::vector<Epetra_Map> &InputMaps)
    {
      reinit (InputMaps);
    }



    inline
    BlockVector::BlockVector (const unsigned int num_blocks)
    {
      reinit (num_blocks);
    }



    inline
    BlockVector::BlockVector (const BlockVector& v)
                    :
                    BlockVectorBase<Vector > ()
    {
      this->components.resize (v.n_blocks());
      this->block_indices = v.block_indices;
    
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }



    inline
    BlockVector &
    BlockVector::operator = (const value_type s)
    {
      BaseClass::operator = (s);
      return *this;
    }



    inline
    BlockVector &
    BlockVector::operator = (const BlockVector &v)
    {
      if (this->n_blocks() != v.n_blocks())
	reinit(v.n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
	this->components[i] = v.block(i);

      collect_sizes();
	
      return *this;
    }



    template <typename Number>
    inline
    BlockVector &
    BlockVector::operator = (const ::dealii::BlockVector<Number> &v)
    {
      Assert (n_blocks() == v.n_blocks(),
	      ExcDimensionMismatch(n_blocks(),v.n_blocks()));

      for (unsigned int i=0; i<this->n_blocks(); ++i)
	this->components[i] = v.block(i);

      return *this;
    }



    inline
    BlockVector::~BlockVector ()
    {}



    inline
    void
    BlockVector::reinit (const std::vector<Epetra_Map> &input_maps)
    {
      unsigned int no_blocks = input_maps.size();
      std::vector<unsigned int> block_sizes (no_blocks);

      for (unsigned int i=0; i<no_blocks; ++i)
	{
	  block_sizes[i] = input_maps[i].NumGlobalElements();
	}

      this->block_indices.reinit (block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->components[i].reinit(input_maps[i]);

      collect_sizes();
    }



    inline
    void
    BlockVector::reinit (const BlockVector& v,
			 const bool fast)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());
  
      for (unsigned int i=0;i<this->n_blocks();++i)
        block(i).reinit(v.block(i), fast);
      
      collect_sizes();
    }



    inline
    void
    BlockVector::reinit (const unsigned int num_blocks)
    {
      std::vector<unsigned int> block_sizes (num_blocks, 0);
      this->block_indices.reinit (block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());
  
      for (unsigned int i=0;i<this->n_blocks();++i)
        block(i).clear();
    }



    inline
    void
    BlockVector::compress ()
    {
      for (unsigned int i=0; i<n_blocks(); ++i)
	block(i).compress();
    }



    inline
    void
    BlockVector::swap (BlockVector &v)
    {
      Assert (this->n_blocks() == v.n_blocks(),
              ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));
  
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->components[i].swap (v.components[i]);
      ::dealii::swap (this->block_indices, v.block_indices);
    }
  


/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors. 
 *
 * @ingroup TrilinosWrappers
 * @relates TrilinosWrappers::BlockVector
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
 */
    inline
    void swap (BlockVector &u,
               BlockVector &v)
    {
      u.swap (v);
    }

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif  // DEAL_II_USE_TRILINOS

#endif

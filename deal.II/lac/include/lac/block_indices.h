//----------------------------  block_indices.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_indices.h  ---------------------------
#ifndef __deal2__block_indices_h
#define __deal2__block_indices_h


#include <vector>


/**
 * Class that manages the conversion of global indices into a block
 * vector or matrix to the local indices within this block. This is
 * required when you address a global element in a block vector and
 * want to know which element within which block this is. It is also
 * useful if a matrix is composed of several blocks, where you have to
 * translate global row and column indices to local ones.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 2000
 */
class BlockIndices
{
  public:

				     /**
				      * Default
				      * constructor. Initialize for
				      * @p{n_blocks} blocks and set
				      * all block sizes to zero.
				      */
    BlockIndices (const unsigned int n_blocks);

				     /**
				      * Constructor. Initialize the
				      * number of entries in each
				      * block @p{i} as @p{n[i]}. The
				      * number of blocks will be the
				      * size of the vector
				      */
    BlockIndices (const vector<unsigned int> &n);
    
				     /**
				      * Reinitialize the number of
				      * indices within each block from
				      * the given argument. The number
				      * of blocks will be adjusted to
				      * the size of @p{n} and the size
				      * of block @p{i} is set to
				      * @p{n[i]}.
				      */
    void reinit (const vector<unsigned int> &n);
    
				     /**
				      * Return the block and the
				      * index within that block
				      * for the global index @p{i}. The
				      * first element of the pair is
				      * the block, the second the
				      * index within it.
				      */
    pair<unsigned int,unsigned int>
    global_to_local (const unsigned int i) const;

				     /**
				      * Return the global index of
				      * @p{index} in block @p{block}.
				      */
    unsigned int local_to_global (const unsigned int block,
				  const unsigned int index) const;

				     /**
				      * Number of blocks in index field.
				      */
    unsigned int size () const;
  
				     /**
				      * Return the total number of
				      * indices accumulated over all
				      * blocks, that is, the dimension
				      * of the vector space of the
				      * block vector.
				      */
    unsigned int total_size () const;

				     /**
				      * Copy operator.
				      */
    BlockIndices & operator = (const BlockIndices &b);

				     /**
				      * Compare whether two objects
				      * are the same, i.e. whether the
				      * starting indices of all blocks
				      * are equal.
				      */
    bool operator == (const BlockIndices &b) const;
    
				     /**
				      * Swap the contents of these two
				      * objects.
				      */
    void swap (BlockIndices &b);
    
  private:
				     /**
				      * Number of blocks. This is made
				      * constant to avoid accidental
				      * changes during lifetime.
				      */
    unsigned int n_blocks;

                                     /**
				      * Global starting index of each
				      * vector. The last and redundant
				      * value is the total number of
				      * entries.
				      */
    vector<unsigned int> start_indices;
};



/* ---------------------- template and inline functions ------------------- */

inline
BlockIndices::BlockIndices (unsigned int n_blocks)
  : n_blocks(n_blocks),
    start_indices(n_blocks)
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    start_indices[i] = 0;
};



inline
BlockIndices::BlockIndices (const vector<unsigned int> &n)
{
  reinit (n);
};



inline
void
BlockIndices::reinit (const vector<unsigned int> &n)
{
  if (start_indices.size() != n.size()+1)
    {
      n_blocks = n.size();
      start_indices.resize(n_blocks+1);
    }
  start_indices[0] = 0;
  for (unsigned int i=1; i<=n_blocks; ++i)
    start_indices[i] = start_indices[i-1] + n[i-1];
};



inline
pair<unsigned int,unsigned int>
BlockIndices::global_to_local (const unsigned int i) const
{
  Assert (i<total_size(), ExcIndexRange(i, 0, total_size()));

  int block = n_blocks-1;
  while (i < start_indices[block])
    --block;

  return make_pair<unsigned int>(block, i-start_indices[block]);
};


inline
unsigned int
BlockIndices::local_to_global (const unsigned int block,
					 const unsigned int index) const
{
  Assert (block < n_blocks, ExcIndexRange(block, 0, n_blocks));
  Assert (index < start_indices[block+1]-start_indices[block],
	  ExcIndexRange (index, 0, start_indices[block+1]-start_indices[block]));

  return start_indices[block]+index;
};


inline
unsigned int
BlockIndices::size () const
{
  return n_blocks;
}



inline
unsigned int
BlockIndices::total_size () const
{
  return start_indices[n_blocks];
};



inline
BlockIndices &
BlockIndices::operator = (const BlockIndices &b)
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    start_indices[i] = b.start_indices[i];
  return *this;
};



inline
bool
BlockIndices::operator == (const BlockIndices &b) const
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    if (start_indices[i] != b.start_indices[i])
      return false;
  
  return true;
};



inline
void
BlockIndices::swap (BlockIndices &b)
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    std::swap (start_indices[i], b.start_indices[i]);
};



/* ----------------- global functions ---------------------------- */


/**
 * Global function @p{swap} which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two objects.
 *
 * @author Wolfgang Bangerth, 2000
 */
inline
void swap (BlockIndices &u, BlockIndices &v)
{
  u.swap (v);
};




#endif

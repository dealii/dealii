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
 * @author Wolfgang Bangerth, 2000
 */
template <int n_blocks>
class BlockIndices
{
  public:

				     /**
				      * Defaulut constructor. Set all
				      * indices denoting the start of
				      * the blocks to zero.
				      */
    BlockIndices ();

				     /**
				      * Constructor. Initialize the
				      * number of indices within each
				      * block from the given
				      * argument. The size of the
				      * vector shall be equal to
				      * #n_blocks#.
				      */
    BlockIndices (const vector<unsigned int> &n);
    
				     /**
				      * Reset the number of indices
				      * within each block from the
				      * given argument. The size of
				      * the vector shall be equal to
				      * #n_blocks#.
				      */
    void reinit (const vector<unsigned int> &n);
    
				     /**
				      * Return the block and the
				      * index within that block
				      * for the global index #i#. The
				      * first element of the pair is
				      * the block, the second the
				      * index within that.
				      */
    pair<unsigned int,unsigned int> global_to_local (const unsigned int i) const;

				     /**
				      * Return the total number of
				      * indices accumulated over all
				      * blocks.
				      */
    unsigned int total_size () const;

				     /**
				      * Copy operator.
				      */
    BlockIndices<n_blocks> & operator = (const BlockIndices<n_blocks> &b);

				     /**
				      * Swap the contents of these two
				      * objects.
				      */
    void swap (BlockIndices<n_blocks> &b);
    
  private:
                                     /**
				      * Global starting index of each
				      * vector. The last and redundant
				      * value is the total number of
				      * entries.
				      */
    unsigned int start_indices[n_blocks+1];
};



/* ---------------------- template and inline functions ------------------- */

template <int n_blocks>
inline
BlockIndices<n_blocks>::BlockIndices ()
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    start_indices[i] = 0;
};



template <int n_blocks>
inline
BlockIndices<n_blocks>::BlockIndices (const vector<unsigned int> &n)
{
  reinit (n);
};



template <int n_blocks>
inline
void
BlockIndices<n_blocks>::reinit (const vector<unsigned int> &n)
{
  Assert(n.size()==n_blocks,
	 ExcDimensionMismatch(n.size(), n_blocks));

  start_indices[0] = 0;
  for (unsigned int i=1; i<=n_blocks; ++i)
    start_indices[i] = start_indices[i-1] + n[i];
};



template <int n_blocks>
inline
pair<unsigned int,unsigned int>
BlockIndices<n_blocks>::find (const unsigned int i) const
{
  Assert (i<total_size(), ExcIndexRange(i, 0, total_size()));

  int block = n_blocks-1;
  while (i < start_indices[block])
    --block;

  return make_pair<unsigned int>(block, i-start_indices[block]);
};



template <int n_blocks>
inline
unsigned int
BlockIndices<n_blocks>::total_size () const
{
  return start_indices[n_blocks];
};



template <int n_blocks>
inline
BlockIndices<n_blocks> &
BlockIndices<n_blocks>::operator = (const BlockIndices<n_blocks> &b)
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    start_indices[i] = b.start_indices[i];
  return *this;
};



template <int n_blocks>
inline
void
BlockIndices<n_blocks>::swap (BlockIndices<n_blocks> &b)
{
  for (unsigned int i=0; i<=n_blocks; ++i)
    std::swap (start_indices[i], b.start_indices[i]);
};



/* ----------------- global functions ---------------------------- */


/**
 * Global function #swap# which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two objects.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int n_blocks>
inline
void swap (BlockIndices<n_blocks> &u,
	   BlockIndices<n_blocks> &v)
{
  u.swap (v);
};




#endif
